import os
# from .dysts import *

data_dirpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

data_path = os.path.join(data_dirpath, 'chaotic_attractors.json')

data_path2 = os.path.join(data_dirpath, 'discrete_maps.json')


from pathlib import Path
PACKAGEDIR = Path(__file__).parent.absolute()



# my_file = pkg_resources.resource_filename('my_data_pack', 'my_data/data_file.txt')

# with open(my_file2) as fin:
#     my_data_object = fin.readlines()