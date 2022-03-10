# list of groups or datasets in the root of the input HDF5 file

import os
import re
import h5py
import sys
from datetime import datetime

if len(sys.argv)<2:
    print(("Usage: {0:s} input_file").format(sys.argv[0]))
    sys.exit()
          
fileName=sys.argv[1]

with h5py.File(fileName, "r") as f:
    for key in list(f.keys()):
        print(key)
