import os

ls_dir=os.listdir('.')

for d in ls_dir:
    if os.path.isdir(d):
        ls_file=os.listdir(d)
        #print(ls_file)
        for f in ls_file:
            if f[0:8]=='run-admm':
                print('./'+d+'/'+f)
            
