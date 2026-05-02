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

def apply_photo_corrections(df, ephemeris):

def covert_times(df, ra, dec, observatory):

def sigma_filter(df, col = 'corected_mag', sigma = 3):



def process_lightcurve(data, metadata, asteroid_id = None):
    """Apply corrections one by one"""
    
    time_corrected, phase_corrected, time_corrected = 






