import io
import numpy as np
import time
import pandas as pd
import astropy.units as u

from astroquery.jplhorizons import Horizons, conf
from astropy.table import vstack
from sbpy.photometry import HG


conf.horizons_server = 'https://ssd.jpl.nasa.gov/api/horizons.api'

def load_all_multifile(folder, extension, id_string=None, time_col = 0, data_col = 1):
    """"""
    file_list = io.create_file_list(folder, extension, id_string)
    if not file_list:
        print(f'No files with extension {extension} found in folder: {folder}\n')
        return None, None, -1
    
    data, metadata, exit_status = io.multi_lightcurve_reader(file_list, time_col, data_col)

    return data, metadata, exit_status

def load_all_single(file_path, time_col = 0, data_col = 1):
    return io.lightcurve_reader(file_path, time_col, data_col)

def fetch_ephemeris(time_jds, asteroid_id, coordinates = None, chunk_size=50):
    """
    Query JPL Horizons to get RA, DEC, phase angel, etc.
    chunk_size = 50 seems to be a good default, but it's worth exploring other parameters
    statistics for query times should be in the wiki
    """
    if coordinates is not None:
        location = {
            'lon': coordinates[1],
            'lat': coordinates [0],
            'elevation' : 0 #TODO: Add a method for calcualte elevation from coordinates, not a big deal right now
        }
    else: 
        location ='500' # Falback to the center of the Earth

    chunks = [time_jds[i:i+chunk_size] for i in range(0, len(time_jds), chunk_size)]

    results=[]

    for i, chunk in enumerate(chunks):
        for attempt in range(3):
            try:
                obj = Horizons(id=str(asteroid_id), location=location, epochs=chunk)
                eph = obj.ephemerides(quantities='1,2,19,20,23,24')
                results.append(eph)
                time.sleep(0.3)
                break
            except Exception as e:
                print(f'Chunk {i} attempt {attempt+1} failed: {e}')
                time.sleep(1)
        else:
            print(f'Chunk {i} failed after 3 attempts, skipping.')

    return vstack(results)


def apply_photo_corrections(df, ephemeris, G=0.15):
    """Correct for both phase angle variations (from Earth) in magnitude and distance variations (from Sun)"""
    df = df.copy()

    r = ephemeris['r'].data.astype(float)
    delta = ephemeris['delta'].data.astype(float)
    alpha = ephemeris['alpha'].data.astype(float)

    distance_correction = 5 * np.log10(r * delta)

    hg = HG(H=0*u.mag, G=G*u.dimensionless_unscaled)
    phase_correction = hg(alpha * u.deg).value

    df[1] = df[1] - distance_correction - phase_correction

    return df

def convert_times(df, ephemeris):
    """Converts times from time recieved at Earth to time emitted from asteroid"""
    df = df.copy()
    c_au_per_day = 173.145
    delta = ephemeris['delta'].data.astype(float)

    df['jd_corrected'] = df[0].values - (delta / c_au_per_day)

    return df


def sigma_filter(df, data_col=1, sigma=3, window=20):
    """
    Apply running sigma clipping using a local window.
    Removes points where |y - local_mean| > sigma * local_std
    """
    df = df.copy()

    while True:
        # Compute rolling mean and std
        rolling_mean = df[data_col].rolling(window=window, center=True, min_periods=3).mean()
        rolling_std = df[data_col].rolling(window=window, center=True, min_periods=3).std()

        mask = np.abs(df[data_col] - rolling_mean) <= sigma * rolling_std

        new_df = df[mask]

        if len(new_df) == len(df):
            break

        df = new_df.reset_index(drop=True)

    return df.reset_index(drop=True)


def check_corrections(metadata):
    """
    Check the metadata to see if the correction have been applied, default to False
    TODO: This method is worth revisiting in the future to handle edge cases
    """
    time_corrected, phase_corrected = False, False
    try:
        time_corrected = metadata['LTCAPP'] != "NONE"
    except KeyError:
        pass

    try:
        phase_corrected = metadata['REDUCEDMAGS'] != "NONE"
    except KeyError:
        pass

    return time_corrected, phase_corrected

def check_relative_mags(metadata):
    """"""
    relative = False

    try:
        relative = metadata['DIFFERMAGS'] == "TRUE"
    except KeyError:
        pass

    return relative
        
     
def load_coordinates(metadata):
    try:
        long = float(metadata['OBSLONGITUDE'])
        lat = float(metadata['OBSLATITUDE'])
    except KeyError:
        print('Coordiantes not in the metadata, please supply them as a flag.')
        return None

    return (lat, long)

def undo_ltcapp(jds, ltcdays):
    """Undo light travel time correction using the value stored in metadata"""
    return jds - float(ltcdays) # The subtraction might be confusing since the correction is always a negative value, but we're undoing the correction

def fetch_asteroid_id(metadata):
    return metadata.get('OBJECTNUMBER') or metadata.get('OBJECTNAME', None) # Try both options, return None otherwise

def process_lightcurve(data, metadata, asteroid_id = None):
    """Apply corrections one by one"""

    processed = []

    if asteroid_id is None:
        asteroid_id = fetch_asteroid_id(metadata[0]) # This assumes all data is for the same asteroid

    # First, we check what has already been corrected, since our ephemeris calls rely on uncorrected data
    # We need to undo any light time travel corrections that have already been applied
    for observation, obs_meta in zip(data, metadata):
        df = observation.copy()
        time_corrected, _ = check_corrections(obs_meta)

        if time_corrected:
            ltcdays = float(obs_meta.get('LTCDAYS', 0))
            df[0] = undo_ltcapp(df[0].values, ltcdays)

        processed.append(df)
    
    all_jds = np.concatenate([df[0].values for df in processed])
    ephemeris = fetch_ephemeris(all_jds, asteroid_id)

    sizes = [len(df) for df in processed]
    split_indices = np.cumsum(sizes)[:-1]
    eph_splits = np.split(np.arange(len(ephemeris)), split_indices)

    final = []

    for df, obs_meta, eph_idx in zip(processed, metadata, eph_splits):
        obs_eph = ephemeris[eph_idx]
        _, phase_corrected = check_corrections(obs_meta)
        coordinates = load_coordinates(obs_meta)

        if not phase_corrected:
            df = apply_photo_corrections(df, obs_eph)

        df = convert_times(df, obs_eph)

        #TODO: Find a better way to homogenize the zero point of each session
        if check_relative_mags(obs_meta):
            pass
        df[1] = df[1] - df[1].mean() # Remove the variable 0 point offset from each session, decent approximation
        df = sigma_filter(df)

        final.append(df)

    combined = pd.concat(final, ignore_index = True)
    
    return combined




