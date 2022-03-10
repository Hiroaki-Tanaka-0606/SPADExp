# SPADExp make database of atomic potentials

import os
import re
import h5py
from datetime import datetime

with h5py.File("Atomic_potentials.hdf5", "w") as f:
    f.attrs.create("datetime", datetime.now().isoformat(" "))
    for Z in range(1, 87):
        Z_string=("{0:03d}").format(Z)
        for h5Name in os.listdir("./Output_files"):
            if(h5Name[0:3]==Z_string):
                print("Adding "+h5Name)
                g=f.create_group(Z_string)
                with h5py.File(os.path.join("./Output_files", h5Name)) as f_input:
                    for key in f_input[Z_string].keys():
                        print(" adding dataset "+key)
                        if(key=="Orbitals"):
                            orbitals=[]
                            for i in range(0, f_input[Z_string][key].shape[0]):
                                orbitals.append(f_input[Z_string][key][i].decode("utf_8"))
                            g.create_dataset(key, data=orbitals)
                        else:
                            g.create_dataset(key, data=f_input[Z_string][key])
                    for key in f_input[Z_string].attrs.keys():
                        print(" adding attribute "+key)
                        g.attrs.create(key, f_input[Z_string].attrs[key])
                break
        

    
