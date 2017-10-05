from glob import glob
import os 

# list_f= glob('../figure/VI/UpdateRule/**/*.pdf', recursive=True)
# list_f= glob('../figure/NSGS/LocalSolver/**/*.pdf', recursive=True)


# list_f= glob('../figure/NSGS/LocalTol/VI/**/*.pdf', recursive=True)
# list_f= glob('../figure/NSGS/LocalTol/**/*.pdf', recursive=True)

# list_f_copy=[]
# for f in list_f:
#     if 'VI' in f:
#         continue
#     else:
#         list_f_copy.append(f)
# list_f=list_f_copy

#list_f= glob('../figure/NSGS/rho/**/*.pdf', recursive=True)
#list_f= glob('../figure/NSGS/LocalSolverHybrid/**/*.pdf', recursive=True)
list_f= glob('../figure/NSGS/Shuffled/**/*.pdf', recursive=True)

#list_f= glob('../figure/PSOR/**/*.pdf', recursive=True)
#list_f= glob('../figure/NSN/**/*.pdf', recursive=True)
list_f= glob('../figure/PROX/NSN/InternalSolvers/**/*.pdf', recursive=True)
#list_f= glob('../figure/PROX/NSGS/InternalSolvers/**/*.pdf', recursive=True)
#list_f= glob('../figure/PROX/Parameters/nu05/**/*.pdf', recursive=True)
#list_f= glob('../figure/OPTI/**/*.pdf', recursive=True)
#list_f= glob('../figure/COMP/**/*.pdf', recursive=True)
#print('list_f=',list_f)
dict_f={}



for f in list_f:
    print('f=',f)
            
    if "legend" in f:
        legendfile = f
    else:

        f_s=f.split("/")
        #print("f_s",f_s)
        #print(f_s[-3])
        #print(f_s[-4])        

        precision = f_s[-4].replace('1.0e','10^{{')+ '}}'
        example_name = (os.path.basename(f).partition('-')[2]).partition('.')[0]

        if example_name == 'BoxesStack1' :
            example_name = 'BoxesStack'
            
        #print("precision=",float(f_s[-4]))
        #print("example_name",example_name)
        
        if example_name == 'LMGC_Bridge_PR' and float(f_s[-4]) == 1e-04 :
            example_name = example_name + ' II'
        if example_name == 'LMGC_Cubes_H8' and float(f_s[-4]) == 1e-04 :
            example_name = example_name + ' II'
        if example_name == 'Chute_local_problems' and float(f_s[-4]) == 1e-04 :
            example_name = example_name + ' II'
        if example_name == 'KaplasTower' and float(f_s[-4]) == 1e-04 :
            example_name = example_name + ' II'
        if example_name == 'LMGC_100_PR_PerioBox' and float(f_s[-4]) == 1e-04 :
            example_name = example_name + ' II'

        
            
        key=example_name
        #print("key=", key)
        dict_f[key] = dict()
        dict_f[key]['timeout']=f_s[-3]
        dict_f[key]['precision']=precision
        dict_f[key]['full_name']=f


        if 'LMGC_' in example_name:
            example_name=example_name[5:]
        example_name=example_name.replace('_','\_')  
        #print(example_name)
        dict_f[key]['example_name']=example_name

   
    
print(dict_f)
count =0
for k in sorted(dict_f, reverse=True):
        count = count +1
        if count%3 ==0:
                extension = '\\\\'
        else:
                extension = ''
        #print('\\subfloat[\\scriptsize {0}  precision ${1}$ timeout ${2}$]{{\\includegraphics[width=\\figwidth]{{{3}}}}} '.format(dict_f[k]['example_name'],dict_f[k]['precision'],dict_f[k]['timeout'],dict_f[k]['full_name'])+extension)
        print('\\subfloat[\\scriptsize {0}]\n   {{\\includegraphics[width=\\figwidth]{{{3}}}}} '.format(dict_f[k]['example_name'],dict_f[k]['precision'],dict_f[k]['timeout'],dict_f[k]['full_name'])+extension)

print('{{\\includegraphics[height=\\legendheight]{{{0}}}}} '.format(legendfile))
