from pathlib import Path


def lightcurve_reader(filepath):
    """"""
    try:
        path = Path(filepath)
        extension = path.suffix
    except Exception as e:
        print(f"Invalid filepath: {filepath}! \nPlease retry with a valid path.\nException: {e}\n")
        return None, -1
    
    match str(extension):
        case '.xml': 
            data, exit_status = xml_handler(path)
        case '.json':
            data, exit_status = json_handler(path)
        case '.csv':
            data, exit_status = csv_handler(path)
        case ".tab":
            data, exit_status = tab_handler(path)
        case _:
            print(f'Extension: {extension} of file {filepath} is not valid! \nPlease enter a valid filetype.\n')
            return None, -1
        
    if exit_status < 0:
        print(f'Reading of file: {filepath} failed.')
        return None, -1
    
    return data, 0

def xml_handler(path):
    return 100, -1

def json_handler(path):
    return 200, 0

def csv_handler(path):
    return 300, 0

def tab_handler(path):
    return 400, 0

def create_file_list(folder_path, extension):
    """Create a file list from a given folder and extension"""
    if ('.' not in extension):
        extension = '.' + extension
    
    return list((Path(folder_path).glob(f"*{extension}")))