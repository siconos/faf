import os
import shutil as shutil
from glob import glob 
base = './TolvaChute_selected'
os.mkdir(base)
counter =0
max_size = 0

def attributes_split(filename):
    attributes=filename.split('-')
    print attributes
    id = int(attributes[-1].split('.')[0])
    print id
    size = int(attributes[-2])
    print size
    iteration = int(attributes[-3][1:])
    print iteration
    return id, size, iteration


for filename in glob('./TolvaChute_size/*.hdf5'):
    print( filename)
    id, size, iteration = attributes_split(filename)
    max_size=max(size,max_size)
    print "max_size", max_size



n_packets=20
n_files =10
dnp = max_size/n_packets
print "dnp",dnp

list_filename= []
for i in range(n_packets):
    list_filename.append([])

    
for filename in glob('./TolvaChute_size/*.hdf5'):
    #print( filename)
    id, size, iteration = attributes_split(filename)
    print 'size/dnp',size/dnp
    list_filename[size/dnp-1].append((iteration,filename))

for i in range(n_packets):
    print('len(list_filename[',i,'])',len(list_filename[i]))
    list_filename[i]= sorted(list_filename[i], key=lambda data: data[0])
    print (list_filename[i][-(n_files+1):-1])
    list_filename[i]=list_filename[i][-(n_files+1):-1]
    print('len(list_filename[',i,'])',len(list_filename[i]))

    for size_f,f in  list_filename[i]:
        print size_f,f
        new_filename = f.split('/')[-1]
        print "copy", f, "in ", os.path.join(base, new_filename)
        shutil.copy(f, os.path.join(base, new_filename))
    #statinfo = os.stat(os.path.join(dirname, filename))
    #print statinfo.st_size,  statinfo_previous.st_size
    #new_filename = filename.replace('FC3D',base)
    # if counter%1==0 :
    #     print "copy", filename, "in ", os.path.join(base, new_filename) 
    #     shutil.copy(filename, os.path.join(base, new_filename))
    # counter =counter +1
