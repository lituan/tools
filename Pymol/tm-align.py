"""
use TM-align to align structures
"""
def get_names(directory):

    def files_in_dir(directory):
        for root, dirs, files in os.walk(directory):
            for f in files:
                if '.pdb' in f or '.PDB' in f:
                    yield os.path.join(root, f)

    names = []
    for f in files_in_dir(directory):
        f_path,f_name = os.path.split(f)
        f_name,f_extention = os.path.splitext(f_name)
        names.append(f_name)
    return names

