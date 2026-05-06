from pathlib import Path
import pandas as pd


def lightcurve_reader(filepath, time_col = 0, data_col = 1):
    """"""
    try:
        path = Path(filepath)
        extension = path.suffix
    except Exception as e:
        print(f"Invalid filepath: {filepath}! \nPlease retry with a valid path.\nException: {e}\n")
        return None, -1
    
    match str(extension):
        case '.xml': 
            data, metadata, exit_status = xml_handler(path, time_col, data_col)
        case '.json':
            data, metadata, exit_status = json_handler(path, time_col, data_col)
        case '.csv':
            data, metadata, exit_status = csv_handler(path, time_col, data_col)
        case ".tab":
            data, metadata, exit_status = tab_handler(path, time_col, data_col)
        case ".txt":
            data, metadata, exit_status = txt_handler(path, time_col, data_col)
        case _:
            # In normal operations, this should never run. Wrapper methods check for this as well, but it is a good double check I suppose.
            print(f'Extension: {extension} of file {filepath} is not valid! \nPlease enter a valid filetype.\n')
            return None, None, -1
        
    if exit_status < 0:
        print(f'Reading of file: {filepath} failed.')
        return None, None, -1
    
    if not isinstance(data, list):
        data = [data]
        metadata = [metadata]

    return data, metadata, 0

def xml_handler(file_path, time_col = 0, data_col = 1):
    return 100, None, -1

def json_handler(file_path, time_col = 0, data_col = 1):
    return 200, None, 0

def csv_handler(file_path, time_col = 0, data_col = 1):
    return 300, None, 0

def tab_handler(file_path, time_col=0, data_col=1):
    """Read .tab files"""
    try:
        df = pd.read_csv(file_path, sep=r"\s+", header=None)
    except Exception as e:
        print(f'Error reading file: {str(file_path)}\n Exception: {e}\n')
        return None, None, -1

    # Tab files have no embedded metadata, so we stub it out with known defaults.
    # Tab files from this pipeline are assumed to be pre-calibrated differential magnitudes
    # with no light travel time correction applied.
    metadata = {
        'LTCAPP': 'NONE',
        'REDUCEDMAGS': 'NONE',
        'DIFFERMAGS': 'FALSE',
        'TIME_FORMAT': 'MJD',
        'source_file': str(file_path),
        'origin': 'tab'
    }

    return df, metadata, 0

def txt_handler(file_path, time_col = 0, data_col = 1):
    """Read txt files that have metadata built in"""
    try:
        text = file_path.read_text()
    except Exception as e:
        print(f"[txt_handler] Failed to read file: {e}")
        return [], [], 1
 
    raw_blocks = text.split("STARTMETADATA")
    raw_blocks = [b.strip() for b in raw_blocks if b.strip()]
 
    if not raw_blocks:
        print("[txt_handler] No STARTMETADATA blocks found.")
        return [], [], 1
 
    all_metadata = []
    all_data = []
 
    for block in raw_blocks:
        meta_raw, _, data_raw = block.partition("ENDMETADATA")
 
        meta_dict = {}
        for line in meta_raw.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                key, _, value = line.partition("=")
                meta_dict[key.strip()] = value.strip()
 
        rows = []
        in_data = False
        for line in data_raw.splitlines():
            line = line.strip()
            if line == "ENDDATA":
                break
            if line.startswith("DATA="):
                in_data = True
                payload = line[len("DATA="):]
                parts = [p.strip() for p in payload.split("|")]
                if len(parts) >= 2:
                    try:
                        jd   = float(parts[0])
                        mag  = float(parts[1])
                        unc  = float(parts[2]) if len(parts) >= 3 else float("nan")
                        rows.append({time_col: jd, data_col: mag, "uncertainty": unc})
                    except ValueError:
                        continue
 
        df = pd.DataFrame(rows, columns=[time_col, data_col, "uncertainty"])
        all_metadata.append(meta_dict)
        all_data.append(df)
 
    return all_data, all_metadata, 0


def create_file_list(folder_path, extension, id_string = None):
    """Create a file list from a given folder and extension"""
    if ('.' not in extension):
        extension = '.' + extension
    
    file_list = list((Path(folder_path).glob(f"*{extension}")))
    if id_string is not None:
        file_list = [f for f in file_list if id_string in f.name]
    
    return file_list

def multi_lightcurve_reader(file_list, time_col=0, data_col=1):
    """Process multipole lightcurves"""
    data = []
    metadata = []
    exit_status = 0

    for file in file_list:
        file_data, file_metadata, exit_status = lightcurve_reader(file, time_col, data_col)
        if exit_status == 0:
            data.extend(file_data)
            metadata.extend(file_metadata)
            continue

        exit_status = 1
    
    #combined_df = pd.concat(data, ignore_index=True) deprecated, delete later?
    #combined_df = combined_df.sort_values(by=time_col, ignore_index=True)

    return data, metadata, exit_status
