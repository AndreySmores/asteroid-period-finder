import io

def load_all_multifile(folder, extension, id_string=None, time_col = 0, data_col = 1):
    """"""
    file_list = io.create_file_list(folder, extension, id_string)
    if file_list == None:
        print(f'No files with extension {extension} found in folder: {folder}\n')
        return None, None, -1
    
    data, metadata, exit_status = io.multi_lightcurve_reader(file_list, time_col, data_col)

    return data, metadata, exit_status

def load_all_single(file_path, time_col = 0, data_col = 1):
    return io.lightcurve_reader(file_path, time_col, data_col)

def fetch_ephemeris(time_jds, asteroid_id, coordinates):
    pass

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
        print(metadata['LTCAPP'])
    except KeyError:
        pass

    try:
        phase_corrected = metadata['REDUCEDMAGS'] != "NONE"
    except KeyError:
        pass

    return time_corrected, phase_corrected
    




def process_lightcurve(data, metadata, asteroid_id = None):
    """Apply corrections one by one"""

    time_corrected_array =  []
    phase_corrected_array = []

    #First, we check what has already been corrected, then we apply the necessary methods
    for observation, obs_meta in zip(data, metadata):
        time_corrected, phase_corrected = check_corrections(obs_meta)
        time_corrected_array.append(time_corrected)
        phase_corrected_array.append(phase_corrected)
    
    return time_corrected_array, phase_corrected_array
        





