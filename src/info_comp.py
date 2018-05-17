#!/usr/bin/env python
import h5py

def long_substr(data):
    substr = ''
    #print("data=",data)
    if len(data) > 0 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and is_substr(data[0][i:i+j], data):
                    substr = data[0][i:i+j]
    return substr

def is_substr(find, data):
    if len(data) < 1 and len(find) < 1:
        return False
    for i in range(len(data)):
        if find not in data[i]:
            return False
    return True

with h5py.File('comp.hdf5','r+') as comp_file:
    data = comp_file.get('data')
    precision=data['comp'].attrs.get('precision')
    print("info_comp -- precision :", precision)
    timeout=data['comp'].attrs.get('timeout')
    print("info_comp -- timeout :", timeout)

    hostname=data['comp'].attrs.get('hostname')
    print("info_comp -- hostname :", hostname)
    measure_name=data['comp'].attrs.get('measure_name')
    print("info_comp -- measure_name :", measure_name)


    comp_data = data['comp']
    test_names= []
    for solvername in comp_data:
        #print("  ",solvername)
        filenames =[filename  for filename in comp_data[solvername]]
        #print(filenames)
        test_names.append(long_substr(filenames).partition('-')[0])
    test_name=long_substr(test_names)
        
    print("info_comp -- test_name :", test_name)
    with open('info_comp.txt', "w") as report_file:
        print("precision", file=report_file)
        print(precision, file=report_file)
        print("timeout",file=report_file)
        print(timeout, file=report_file)
        print("test_name", file=report_file)
        print(test_name, file=report_file)

    
