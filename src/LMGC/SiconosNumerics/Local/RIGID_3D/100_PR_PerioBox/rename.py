import os
import shutil as shutil
from glob import glob 
base = '100_PR_PerioBox'
os.mkdir(base)
counter =0
for filename in glob('*.hdf5'):
    print( filename)
    #statinfo = os.stat(os.path.join(dirname, filename))
    #print statinfo.st_size,  statinfo_previous.st_size
    new_filename = filename.replace('FC3D',base)
    if counter%50==0 :
        print "copy", filename, "in ", os.path.join(base, new_filename) 
        shutil.copy(filename, os.path.join(base, new_filename))
    counter =counter +1

