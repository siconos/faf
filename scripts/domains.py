#!/usr/bin/env python
import sys
import h5py
import getopt
import os
def drepr(x, sort = True, indent = 0):
        """ Provides nice print format for a dictionnary """
        if isinstance(x, dict):
                r = '{\n'
                for (key, value) in (sorted(x.items()) if sort else x.iteritems()):
                        r += (' ' * (indent + 4)) + repr(key) + ': '
                        r += drepr(value, sort, indent + 4) + ',\n'
                r = r.rstrip(',\n') + '\n'
                r += (' ' * indent) + '}'
        # elif hasattr(x, '__iter__'):
        #       r = '[\n'
        #       for value in (sorted(x) if sort else x):
        #               r += (' ' * (indent + 4)) + drepr(value, sort, indent + 4) + ',\n'
        #       r = r.rstrip(',\n') + '\n'
        #       r += (' ' * indent) + ']'
        else:
                r = repr(x)
        return r

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
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                   ['help','target=', 'test_name', 'domain', 'precision', 'timeout', 'step'])


    
except getopt.GetoptError as err:
        sys.stderr.write('{0}\n'.format(str(err)))
        usage()
        exit(2)
with h5py.File('comp.hdf5','r+') as comp_file:
    data = comp_file.get('data')
    precision=data['comp'].attrs.get('precision')
    #print("info_comp : precision", precision)
    timeout=data['comp'].attrs.get('timeout')
    #print("info_comp :timeout", timeout)
    comp_data = data['comp']
    test_names= []
    for solvername in comp_data:
        #print("  ",solvername)
        filenames =[filename  for filename in comp_data[solvername]]
        #print(filenames)
        test_names.append(long_substr(filenames).partition('-')[0])
    test_name=long_substr(test_names)
    if test_name.endswith('_'):
        test_name  = test_name[:-1]


# print("domain : test_name", test_name)
# with open('domain.txt', "w") as report_file:
#     print("precision", file=report_file)
#     print(precision, file=report_file)
#     print("timeout",file=report_file)
#     print(timeout, file=report_file)
#     print("test_name", file=report_file)
#     print(test_name, file=report_file)

#print([x for x in next(os.walk('/Users/acary/Work/fclib-library/'))[1]])
tests= ['LMGC_100_PR_PerioBox', 'LMGC_945_SP_Box_PL', 'LMGC_AqueducPR', 'BoxesStack1', 'LMGC_Bridge_PR', 'Capsules', 'Chain', 'Chute_1000', 'Chute_4000', 'Chute_local_problems', 'LMGC_Cubes_H8', 'Global', 'KaplasTower', 'LMGC_LowWall_FEM']
targets= 'vi nsgs_localtol_ac_gp nsgs_localtol_vi  nsgs_localsolver nsgs_localsolver_hybrid nsgs_shuffled psor_solvers nsn_solvers prox_solvers prox_series regul_series opti_solvers comp_solvers comp_solvers_large'
list_target=targets.split(' ')

#print('len of targets', len(list_target))
#default Values
default_values = [5,2,2,4,4,5,5,5,15,100,100,100,10,20,100]
#print('len of default_values', len(default_values))

data=dict()
for t in tests:
    data[t]=dict()
    i=0
    for tt in list_target:
        data[t][tt]=dict()
        data[t][tt]['domain'] = default_values[i]
        i=i+1
        


        
#specific values
test= 'Chain'
data[test]['nsgs_localsolver']['domain']=10
data[test]['nsn_solvers']['domain']=35


test= 'LMGC_AqueducPR'
data[test]['nsgs_localtol_ac_gp']['domain']=4
data[test]['opti_solvers']['domain']=10
data[test]['comp_solvers']['domain']=5
data[test]['comp_solvers_large']['domain']=30

test= 'LMGC_Bridge_PR'
data[test]['nsgs_localtol_ac_gp']['domain']=4
data[test]['nsgs_localsolver']['domain']=5
data[test]['nsgs_shuffled']['domain']=3
data[test]['nsn_solvers']['domain']=100
data[test]['opti_solvers']['domain']=5
data[test]['comp_solvers']['domain']=5

test = 'LMGC_Cubes_H8'
data[test]['vi']['domain']=5
data[test]['nsgs_localtol_ac_gp']['domain']=3
data[test]['nsgs_localsolver']['domain']=5
data[test]['nsgs_shuffled']['domain']=10
data[test]['psor_solvers']['domain']=20
data[test]['nsn_solvers']['domain']=50
data[test]['opti_solvers']['domain']=5
data[test]['comp_solvers']['domain']=5

test= 'LMGC_945_SP_Box_PL'
data[test]['nsgs_localtol_ac_gp']['domain']=4
data[test]['nsgs_localsolver']['domain']=10
data[test]['nsgs_localsolver_hybrid']['domain']=10
data[test]['nsgs_shuffled']['domain']=3
data[test]['psor_solvers']['domain']=10
data[test]['opti_solvers']['domain']=5
data[test]['comp_solvers']['domain']=15


test = 'LMGC_100_PR_PerioBox'
data[test]['vi']['domain']=4
data[test]['psor_solvers']['domain']=15
data[test]['nsn_solvers']['domain']=20
data[test]['prox_solvers']['domain']=2
data[test]['prox_series']['domain']=20
data[test]['nsgs_localsolver']['domain']=15
data[test]['nsgs_localtol_ac_gp']['domain']=8

data[test]['nsgs_localsolver_hybrid']['domain']=25
data[test]['prox_series']['domain']=200
data[test]['regul_series']['domain']=50



test = 'LMGC_LowWall_FEM'
data[test]['vi']['domain']=4
data[test]['nsgs_localtol_ac_gp']['domain']=4
data[test]['nsgs_localsolver']['domain']=5
data[test]['nsgs_localtol_ac_gp']['domain']=3
data[test]['nsgs_localsolver']['domain']=5
data[test]['nsgs_shuffled']['domain']=3
data[test]['psor_solvers']['domain']=15
data[test]['nsn_solvers']['domain']=20
data[test]['prox_solvers']['domain']=20
data[test]['prox_series']['domain']=20
data[test]['comp_solvers']['domain']=2
data[test]['comp_solvers_large']['domain']=10


test = 'Chute_local_problems'
data[test]['prox_solvers']['domain']=3
data[test]['prox_series']['domain']=3
data[test]['comp_solvers']['domain']=10
data[test]['comp_solvers_large']['domain']=50
data[test]['regul_series']['domain']=15
data[test]['nsn_solvers']['domain']=4

test = 'Chute_1000'
#data[test]['nsgs_shuffled']['domain']=
data[test]['psor_solvers']['domain']=10

for t in tests:
    for tt in list_target:
        
        #print(t,tt)
        #print(data[t][tt])
        #print(data[t][tt]['domain'])
        data[t][tt]['step'] =(data[t][tt]['domain'])/100.0
        #print(data[t][tt]['precision'])



#print(drepr(data))

for o, a in opts:
    if (o == '--test_name'):
        print(test_name)
    if (o == '--target'):
        target=str(a)
        
for o, a in opts:
      if (o == '--domain'):
            print(data[test_name][target]['domain'])
      if (o == '--step'):
            print("{0:.3f}".format(data[test_name][target]['step']))
      if (o == '--precision'):
            print("{0:.1e}".format(precision))
      if (o == '--timeout'):
            print("{0:d}".format(int(timeout)))
