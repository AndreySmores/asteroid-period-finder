import io
import numpy as np
from astroquery.jplhorizons import Horizons, conf
from astropy.table import vstack
import time

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
    """Query JPL Horizons to get RA, DEC, phase angel, etc."""
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


def apply_photo_corrections(df, ephemeris):
    pass

def covert_times(df, ra, dec, observatory):
    pass

def sigma_filter(df, col = 'corected_mag', sigma = 3):
    pass

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
     
def load_coordinates(metadata):
    try:
        long = float(metadata['OBSLONGITUDE'])
        lat = float(metadata['OBSLATITUDE'])
    except KeyError:
        print('Coordiantes not in the metadata, please supply them as a flag.')
        return None

    return (lat, long)

def process_lightcurve(data, metadata, asteroid_id = None):
    """Apply corrections one by one"""

    all_jds = np.concatenate([df[0].values for df in data])
    ephemeris = fetch_ephemeris(all_jds, asteroid_id)

    # First, we check what has already been corrected, then we apply the necessary methods
    for observation, obs_meta in zip(data, metadata):
        # TODO Find a way to reduce the number of JPL queries
        # There is always a chance that multiple locations are mixed it, so we can't just
        # Create one query since they are location dependent

        time_corrected, phase_corrected = check_corrections(obs_meta)
        coordinates = load_coordinates(obs_meta)
    return ephemeris




