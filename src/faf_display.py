from __future__ import print_function
import h5py
from faf_tools import *
from faf_matrix_tools import *
class Faf_display():
    def __init__(self, result_file,
                 solvers,
                 time,
                 domain,
                 filenames,
                 filename,
                 rhos= None,
                 gnp_output=False,
                 gnp_with_color=False,
                 gnp_separate_keys=False,
                 no_matplot=False,
                 logscale=False):
        self._result_file = result_file
        self._solvers = solvers
        self._rhos=rhos
        self._domain= domain
        self._time=time
        self._filenames=filenames
        self._filename=filename

        self._gnuplot_output=gnp_output
        self._gnuplot_with_color=gnp_with_color
        self._gnuplot_separate_keys= gnp_separate_keys
        self._gnuplot_add_title = True
        self._no_matplot=no_matplot
        self._logscale=logscale
        
        
    def default_display_task(self, solver_r):
        print("4 -- Running default display tasks ")

        with h5py.File(self._result_file, 'r') as comp_file:
            data = comp_file['data']
            comp_data = data['comp']
            date_str =  self._time.ctime(creation_date(self._result_file))
            
            #print('###########', time.ctime(creation_date('comp.hdf5')))
        
            if (self._gnuplot_output and (self._filename != None)) :
                def write_report(r, filename):
                    with open(filename, "w") as input_file:
                        for k, v in r.items():
                            line = '{}, {}'.format(k, v)
                            print(line, file=input_file)

                out_data=np.empty([len(self._domain),len(comp_data)+1])
                write_report(self._rhos,'rhos.txt')
                write_report(solver_r,'solver_r.txt')
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


                with open('profile.gp','w') as gp:
                    #  all_rhos = [ domain ] + [ rhos[solver_name] for solver_name in comp_data ]
                    all_rhos = [ self._domain ] + [ self._rhos[solver.name()] for solver in filter(lambda s: s._name in comp_data, self._solvers) ]
                    np.savetxt('profile.dat', np.matrix(all_rhos).transpose())
                    gp.write('resultfile = "profile.dat"\n')
                    #print("self._filenames=",self._filenames)
                    #print("long_substr(self._filenames)=", long_substr(self._filenames))
                    test_name = long_substr(self._filenames).partition('-')[0]
                    print("test_name=",test_name)
                    test_name_gnuplot = test_name.replace('_',' ')
                    print("test_name_gnuplot=",test_name_gnuplot)
                    if test_name.endswith('_'):
                        test_name  = test_name[:-1]
                    gp.write('basename="profile-{0}"\n'.format(test_name))
                    #print filename.partition('-')[0]
                    print("test_name=",test_name)
                    gp.write('\n')
                    gp.write('term_choice_tikz=1\n')
                    gp.write('if (term_choice_tikz == 1) \\\n')
                    if (self._gnuplot_with_color):
                        gp.write('set term tikz standalone size 5in,3in font \'\\small\\sf\';  \\\n')
                    else:
                        gp.write('set term tikz standalone monochrome  size 5in,3in font \'\\small\\sf\';  \\\n')
                    gp.write('extension = \'.tex\'; \\\n')
                    gp.write('extension_legend = \'_legend.tex\'; \\\n')
                    gp.write('set output basename.extension; \\\n')
                    gp.write('print "output = ", basename.extension; \\\n')

                    gp.write('else \\\n')
                    gp.write('set term aqua;\\\n')
                    gp.write('\n')

                    if self._gnuplot_add_title:
                        gp.write('set title\'{0} - precision: {1} - timeout: {2} -  {3}\';; \n'.format(test_name_gnuplot,comp_data.attrs['precision'],comp_data.attrs['timeout'], date_str))


                    
                    gp.write('set xrange [{0}:{1}]\n'.format(self._domain[0]-0.01, self._domain[len(self._domain)-1]))
                    gp.write('set yrange [-0.01:1.01]\n')
                    gp.write('set ylabel \'$\\rho(\\tau)$ \' \n')
                    maxrows=len(self._solvers)/2+1
                    gp.write('set key below right vertical maxrows {0}\n'.format(maxrows))


                    x_label=False
                    if x_label:
                        if logscale:
                            gp.write('set logscale x\n')
                            gp.write('set xlabel \'$\\tau$ ({0}) (logscale)\' \n'.format(measure_name))
                        else:
                            gp.write('set xlabel \'$\\tau$ ({0})\' \n'.format(measure_name))
                            
                    #gp.write('set title \'{0}\'\n'.format(filename.partition('-')[0]));
                    gp.write('plot ')
                    if self._gnuplot_separate_keys:
                        if (self._gnuplot_with_color):
                            gp.write(
                                ','.join(['resultfile using 1:{0} notitle w l  dashtype {1} linecolor {2} lw 3'.format(index + 2,index+1,index%6+1)
                                          for index, solver in enumerate(filter(lambda s: s._name in comp_data, self._solvers)) ]))

                            gp.write('\n set output basename.extension_legend; \n')
                            gp.write('print "output = ", basename.extension_legend; \n \n')
                            gp.write('unset border; \n \n')
                            gp.write('unset title; \n \n')
                            gp.write('unset tics; \n \n')
                            gp.write('unset xlabel; \n \n')
                            gp.write('unset ylabel; \n \n')
                            gp.write('set term tikz standalone  size 5in,1.5in font \'\\small\\sf\';  \\\n')
                            gp.write('set key right inside vertical maxrows {0}\n'.format(maxrows))
                            gp.write('\n plot [0:1] [0:1]')
                            gp.write(
                                ','.join([' NaN t "{1}" w l dashtype {2} linecolor {3} lw 3'.format(index + 2, solver.gnuplot_name(),index+1,index%6+1)
                                          for index, solver in enumerate(filter(lambda s: s._name in comp_data, self._solvers)) ]))



                        else:
                            gp.write(
                                ','.join(['resultfile using 1:{0} notitle w l dashtype {1} linecolor {2} lw 3'.format(index + 2,index+1,8)
                                          for index, solver in enumerate(filter(lambda s: s._name in comp_data, self._solvers)) ]))

                            gp.write('\n set output basename.extension_legend; \n')
                            gp.write('print "output = ", basename.extension_legend; \n \n')
                            gp.write('unset border; \n \n')
                            gp.write('unset title; \n \n')
                            gp.write('unset tics; \n \n')
                            gp.write('unset xlabel; \n \n')
                            gp.write('unset ylabel; \n \n')
                            gp.write('set term tikz standalone monochrome  size 5in,1.5in font \'\\small\\sf\';  \\\n')
                            gp.write('set key right inside vertical maxrows {0}\n'.format(maxrows))
                            gp.write('\n plot [0:1] [0:1]')
                            gp.write(
                                ','.join([' NaN t "{1}" w l dashtype {2} linecolor {3} lw 3'.format(index + 2, solver.gnuplot_name(),index+1,8)
                                          for index, solver in enumerate(filter(lambda s: s._name in comp_data, self._solvers)) ]))


                    else:
                        if (self._gnuplot_with_color):
                            gp.write(
                                ','.join(['resultfile using 1:{0} t "{1}" w l dashtype {2} linecolor {3} lw 3'.format(index + 2, solver.gnuplot_name(),index+1,index%6+1)
                                          for index, solver in enumerate(filter(lambda s: s._name in comp_data, self._solvers)) ]))
                        else:
                            gp.write(
                                ','.join(['resultfile using 1:{0} t "{1}" w l dashtype {2} linecolor {3} lw 3'.format(index + 2, solver.gnuplot_name(),index+1,8)
                                          for index, solver in enumerate(filter(lambda s: s._name in comp_data, self._solvers)) ]))

                # all_rhos = [ rhos[solver_name] for solver_name in comp_data ]
                # g.plot(*all_rhos)

            if (self._gnuplot_output and (self._filename == None)) :
                print("Warning: no problem corresponding to the required solver")
                if (os.path.isfile('profile.gp')):
                    os.remove('profile.gp')

            if not self._no_matplot:
                # 5 plot
                from matplotlib.pyplot import subplot, title, plot, grid, show, get_fignums, legend, figure, xlim, ylim, xscale

                #for solver_name in comp_data:
                for solver in self._solvers:
                    solver_name=solver.name()
                    if self._logscale:
                        xscale('log')

                    if solver_name in comp_data :
                        plot(self._domain, self._rhos[solver_name], label=solver_name)
                        ylim(0, 1.0001)
                        xlim(self._domain[0], self._domain[-1])
                        legend(loc=4)
                    grid()
        



    def display_convergence(self,user_filenames, random_sample_proba, max_problems, cond_nc):
        import matplotlib.pyplot as plt
        from matplotlib.pyplot import show
        with h5py.File(self._result_file, 'r') as comp_file:

            data = comp_file['data']
            comp_data = data['comp']

            if user_filenames == []:
                for solver_name in comp_data:

                    filenames = subsample_problems(comp_data[solver_name],
                                                   random_sample_proba,
                                                   max_problems, cond_nc)
            else:
                filenames = user_filenames
            
            for filename in filenames:

                fig, axs = plt.subplots(1, 1)
                ax = axs

                ax.set_title('Convergence on {0}'.format(filename))
                ax.grid(True, which="both")

                ax.set_yscale('symlog', linthreshy=0.001)
                
                for solver_name in comp_data:

                    try:
                        pfilename = os.path.splitext(filename)[0]
                        solver_problem_data = comp_data[solver_name][pfilename]
                    
                        ax.plot(np.arange(len(solver_problem_data['errors'][:])),
                                np.log(solver_problem_data['errors']),
                                label='{0}'.format(solver_name))
                        ax.legend(loc='lower left')
                        
                    except:
                        pass

    def display_distribution(self,problem_filenames):
        from matplotlib.pyplot import title, subplot, grid, show, get_fignums, legend, figure, hist, xlim, ylim, xscale

        nc = []
        nds = []
        cond_nc = []
        cond_W = []
        cond_W_lsmr = []
        rank_dense_W=[]
        rank_estimate_W=[]
        rank_ratio =[]
        rank_estimate_ratio =[]
        for problem_filename in problem_filenames:
            try:
                 nc.append(numberOfDegreeofFreedomContacts(problem_filename)/3)
            except:
                 pass
            try:
                 nds.append(numberOfDegreeofFreedom(problem_filename))
            except:
                 pass
            try:
                 cond_nc.append(cond_problem(problem_filename))
            except:
                 pass
            try:
                 cond_W.append(cond(problem_filename))
            except:
                 pass
            try:
                 cond_W_lsmr.append(cond_lsmr(problem_filename))
            except:
                 pass
            try:
                 rank_dense_W.append(rank_dense(problem_filename))
            except:
                 pass
            try:
                 rank_estimate_W.append(rank_estimate(problem_filename))
            except:
                 pass
            try:
                 rank_ratio.append(numberOfDegreeofFreedomContacts(problem_filename)/float(rank_dense(problem_filename)))
                 rank_estimate_ratio.append(numberOfDegreeofFreedomContacts(problem_filename)/float(rank_estimate(problem_filename)))
            except:
                 pass
        print("number of problems", len(nds))
        print("nds", nds)
        print("max ndof", max(nds))
        print("min ndof", min(nds))
        print("max nc", max(nc))
        print("min nc", min(nc))


        print("rank_dense_W",  rank_dense_W)
        print("cond_nc", cond_nc)

        print("cond_W", cond_W)
        if (len(cond_W) >0):
            print("max cond_W", max(cond_W))
            print("min cond_W", min(cond_W))

            
        print("cond_W_lsmr", cond_W_lsmr)
        print("max cond_W_lsmr", max(cond_W_lsmr))
        print("min cond_W_lsmr", min(cond_W_lsmr))
        print("max cond_nc", max(cond_nc))
        print("min cond_nc", min(cond_nc))
        
        print("max rank_dense_W", max(rank_dense_W))
        print("min rank_dense_W", min(rank_dense_W))
       
        print("max rank_estimate_W", max(rank_estimate_W))
        print("min rank_estimate_W", min(rank_estimate_W))

        print("max rank_ratio", max(rank_ratio))
        print("min rank_ratio", min(rank_ratio))

        print("max rank_estimate_ratio", max(rank_estimate_ratio))
        print("min rank_estimate_ratio", min(rank_estimate_ratio))


        if (len(cond_W) == 0):
            cond_W = cond_W_lsmr



        print(nc)
        figure()
        subplot(311)
        hist(nc, 100, label='nc', histtype='stepfilled')
        grid()
        legend()
        subplot(312)
        import math
        if not math.isnan(min(nds)):
            hist(nds, 100, label='nds', histtype='stepfilled')
        grid()
        legend()
        subplot(313)

        if not math.isnan(min(cond_nc)):
            hist(cond_nc, 100, label='cond_nc', histtype='stepfilled')
        grid()
        legend()

        figure()
        subplot(311)
        if not math.isnan(min(cond_W)):
            hist(cond_W, 100, label='cond(W)', histtype='stepfilled')
        grid()
        legend()
        subplot(312)
        if not math.isnan(min(rank_dense_W)):
            hist(rank_dense_W, 100, label='rank(W)', histtype='stepfilled')
        grid()
        legend()
        subplot(313)

        if not math.isnan(min(rank_ratio)):
            hist(rank_ratio, 100, label='rank_ratio(W)', histtype='stepfilled')
        grid()
        legend()

        if self._gnuplot_output :

            with open('distrib.gp','w') as gp:
                 #  all_rhos = [ domain ] + [ rhos[solver_name] for solver_name in comp_data ]
                 all_distrib = [ nc] + [nds] + [cond_nc]
                 np.savetxt('distrib.dat', np.matrix(all_distrib).transpose())
                 gp.write('resultfile = "distrib.dat"\n')
                 gp.write('basename="distrib-{0}"\n'.format(problem_filename.partition('-')[0]))
                 gp.write('\n')
                 gp.write('ismin(x) = (x<min)?min=x:0\n')
                 gp.write('ismax(x) = (x>max)?max=x:0\n')
                 gp.write('\n')
                 gp.write('max=-1e38;min=1e38;\n')
                 gp.write('plot resultfile u 1:(ismin($1)*ismax($1))\n')
                 gp.write('min_nc = min; max_nc = max\n')
                 gp.write('max=-1e38;min=1e38;\n')
                 gp.write('plot resultfile u 1:(ismin($2)*ismax($2))\n')
                 gp.write('min_ndof = min; max_ndof = max\n')
                 gp.write('max=-1e38;min=1e38;\n')
                 gp.write('plot resultfile u 1:(ismin($3)*ismax($3))\n')
                 gp.write('min_ncond = min; max_ncond = max\n')
                 gp.write('\n')
                 gp.write('\n')

                 gp.write('term_choice_tikz=1\n')
                 gp.write('if (term_choice_tikz == 1) \\\n')
                 gp.write('set term tikz standalone monochrome  size 5in,3in font \'\\small\\sf\';  \\\n')
                 gp.write('extension = \'.tex\'; \\\n')
                 gp.write('set output basename.extension; \\\n')
                 gp.write('print "output = ", basename.extension; \\\n')
                 gp.write('else \\\n')
                 gp.write('set term aqua;\\\n')
                 gp.write(' \n')
                 gp.write('set xtics offset 0,0.5 \n')
                 gp.write('set key left top\n')

                 gp.write('basheight = 0.36; heightoff = 0.0; winratio = 1.0; winheight = basheight*winratio ; trans = 0.9\n')
                 gp.write('set multiplot \n')
                 gp.write('set size winratio,winheight \n')
                 gp.write('\n')
                 # gp.write('set xrange [{0}:{1}]\n'.format(domain[0]-0.01, domain[len(domain)-1]))
                 # gp.write('set yrange [-0.01:1.01]\n')
                 gp.write('set ylabel \'\\shortstack{Number of \\\ problems} \' \n')
                 gp.write('bin(x, width) = width*floor(x/width) + binwidth/2.0\n')
                 #gp.write('set title \'{0}\'\n'.format(problem_filename.partition('-')[0]));
                 gp.write('\n')
                 gp.write('set origin 0.0,winheight*2.0*trans+heightoff\n')
                 gp.write('numberofbox=50\n')
                 gp.write('print \'max_nc =\', max_nc,\' min_nc =\', min_nc \n')
                 gp.write('binwidth = (max_nc-min_nc)/numberofbox\n')
                 gp.write('set boxwidth binwidth\n')

                 gp.write('print \'binwidth =\', binwidth \n')

                 gp.write('set xlabel \'number of contacts\' offset 0,1.2 \n')
                 #gp.write('plot resultfile u (bin($1, binwidth)):(1.0) smooth freq w boxes title \'number of contacts\'  \n')
                 gp.write('plot resultfile u (bin($1, binwidth)):(1.0) smooth freq w boxes notitle  \n')
                 gp.write('\n')
                 gp.write('print \'max_ndof =\', max_ndof,\' min_ndof =\', min_ndof \n')
                 gp.write('if ( (max_ndof-min_ndof) < numberofbox) {binwidth =1 ;} else {binwidth = (max_ndof-min_ndof)/numberofbox}\n');

                 gp.write('set boxwidth binwidth\n')
                 gp.write('set origin 0.0,winheight*1.0*trans+heightoff\n')

                 gp.write('set xlabel \'number of degrees of freedom \' offset 0,1.2 \n')

                 #gp.write('plot resultfile u (bin($2, binwidth)):(1.0) smooth freq w boxes title  \'number of degrees of freedom \' \n')
                 gp.write('plot resultfile u (bin($2, binwidth)):(1.0) smooth freq w boxes notitle \n')
                 gp.write('\n')

                 gp.write('set origin 0.0,winheight*0.0+heightoff\n')
                 gp.write('binwidth = (max_ncond-min_ncond)/numberofbox\n')
                 gp.write('print \'binwidth =\', binwidth \n')

                 gp.write('set boxwidth binwidth\n')
                 gp.write('set xlabel \'ratio number of contacts unknowns/number of degrees of freedom\' offset 0,1.2 \n')
                 #gp.write('plot resultfile u (bin($3, binwidth)):(1.0) smooth freq w boxes title \'ratio number of contacts unknowns/number of degrees of freedom\' \n')
                 gp.write('plot resultfile u (bin($3, binwidth)):(1.0) smooth freq w boxes notitle \n')



                 #gp.write(','.join(['resultfile using 1:{0} t "{1}" w l'.format(index + 2, solver.name())
            #                    for index, solver in enumerate(filter(lambda s: s._name in comp_data, solvers)) ]))
            # all_rhos = [ rhos[solver_name] for solver_name in comp_data ]
            # g.plot(*all_rhos)


    def display_speedup(self):
        from matplotlib.pyplot import subplot, title, plot, grid, show, get_fignums, legend, figure, hist, bar, xlabel, ylabel, boxplot
        print('\n display speedup is starting ...')
        with h5py.File(self._result_file, 'r') as comp_file:

            data = comp_file['data']
            comp_data = data['comp']
            result_utimeout=comp_file['data']['comp'].attrs.get('timeout')
            print('solver in comp_data =',[s for s in comp_data] )
            solvers=[]
            if user_solvers != []:
                #print "user_solvers", user_solvers
                solvers.extend( list(filter(lambda s: any(us in s for us in user_solvers), comp_data)))

                if solvers == []:
                    raise RuntimeError ("Cannot find any matching solver")

            elif user_solvers_exact != []:
                #print "user_solvers_exact", user_solvers_exact
                solvers.extend(list(filter(lambda s: any(us ==  s  for us in user_solvers_exact), comp_data)))

                if solvers == []:
                    raise RuntimeError("Cannot find any solvers in specified list")

            else:
                solvers= comp_data

            print(' filtered solver in comp_data =',[s for s in solvers])

            solvers= list(filter(lambda s: ('OPENMP' in s), solvers))


            print(' solver for speedup-display =',[s for s in solvers])
            if solvers == []:
                raise RuntimeError("Cannot find any solvers in specified list")
            #---#
            # collect results by filename
            #---#
            results_filename={}
            nthread_set= set()
            n_filename=[]

            for solver_name in solvers:
                if user_filenames == []:
                    filenames = subsample_problems(comp_data[solver_name],
                                                   random_sample_proba,
                                                   max_problems, cond_nc)
                else:
                    filenames = user_filenames
                n_filename.append(len(filenames))
                for filename in filenames:
                    pfilename = os.path.splitext(filename)[0]
                    nthread=int(solver_name.split('-')[-1])
                    nthread_set.add(nthread)
                    measure_data =comp_data[solver_name][pfilename].attrs[measure_name]
                    nc =comp_data[solver_name][pfilename].attrs['nc']
                    n_iter =comp_data[solver_name][pfilename].attrs['iter']
                    if filename in results_filename.keys():
                        results_filename[filename].append([solver_name,nthread,measure_data,nc,n_iter])
                    else:
                        results_filename[filename] = [[solver_name,nthread,measure_data,nc,n_iter]]


            #print("n_filename", n_filename)

            results_filename_fails = {}
            for filename,solver in results_filename.items():
                for s in solver:
                    if np.isnan(s[2]):
                        if filename in results_filename.keys():
                            results_filename_fails[filename]=results_filename[filename]
                            print("\nremove failed instance for filename:",filename)
                            print(results_filename.pop(filename))




            nthread_list=list(nthread_set)
            nthread_list.sort()

            # collect results by thread

            measure_by_nthread=[]
            measure_penalized_by_nthread=[]
            measure_mean_by_nthread=[]
            measure_mean_penalized_by_nthread=[]
            index_non_failed=[]
            index_failed=[]

            for n in nthread_list:
                measure_by_nthread.append([])
                measure_mean_by_nthread.append(0.0)

            speedup_list =[]
            speedup_size_list =[]

            for n in nthread_list:
                speedup_list.append([])
                speedup_size_list.append([])

            for filename,solver in results_filename.items():
                for s in solver:
                        nthread= s[1]
                        thread_index=nthread_list.index(nthread)
                        measure_by_nthread[thread_index].append(s[2])
                        measure_mean_by_nthread[thread_index] += s[2]


            for n in nthread_list:
                thread_index    = nthread_list.index(n)
                measure_mean_by_nthread[thread_index] /= len(measure_by_nthread[thread_index])

            #raw_input()
            #print('measure_mean_by_nthread', measure_mean_penalized_by_nthread)

            speedup_list =[]
            speedup_size_list =[]
            iter_size_list=[]
            speedup_avg =[]

            for n in nthread_list:
                speedup_list.append([])
                speedup_size_list.append([])
                iter_size_list.append([])
                speedup_avg.append(0.0)
                thread_index=nthread_list.index(n)
                for i in range(len(measure_by_nthread[thread_index])):
                    speedup_list[thread_index].append(0.0)

            for n in nthread_list:
                thread_index_ref= nthread_list.index(1)
                thread_index    = nthread_list.index(n)
                for i  in range(len(measure_by_nthread[thread_index])):
                    speedup_list[thread_index][i] =  measure_by_nthread[thread_index_ref][i]/measure_by_nthread[thread_index][i]
                    speedup_avg[thread_index] += speedup_list[thread_index][i]
                speedup_avg[thread_index] /=len(measure_by_nthread[thread_index])

            print('speedup_avg', speedup_avg)
            #print('speedup_list', speedup_list)

            _cmp=0
            for filename,solver in results_filename.items():
                for s in solver:
                    thread_index=nthread_list.index(s[1])
                    #print _cmp, speedup_list[thread_index], speedup_list[thread_index][_cmp]
                    speedup_size_list[thread_index].append([s[3] , speedup_list[thread_index][_cmp]])
                    iter_size_list[thread_index].append([s[3] , s[4]])
                _cmp+=1

            count_failed=[]
            for n in nthread_list:
                count_failed.append(0)
            for filename,solver in results_filename_fails.items():
                for s in solver:
                    nthread=  s[1]
                    thread_index=nthread_list.index(nthread)
                    if np.isnan(s[2]):
                        count_failed[thread_index] +=1

            # ---------------------------- #
            # figure #
            # ---------------------------- #



            figure(figsize=(8,14))

            subplot('211')
            xlabel('problem size')
            ylabel('speedup for each thread')
            for n in nthread_list:
                index=nthread_list.index(n)
                list_size_speedup= speedup_size_list[index]
                list_size_speedup= sorted(list_size_speedup, key=lambda data: data[0])
                plot(np.array([t[0]  for t in list_size_speedup]), np.array([t[1]  for t in list_size_speedup]), label='n='+str(n))
            legend()
            subplot('212')
            xlabel('problem size')
            ylabel('iter for each thread')
            for n in nthread_list:
                index=nthread_list.index(n)
                iter_size_speedup= iter_size_list[index]
                iter_size_speedup= sorted(iter_size_speedup, key=lambda data: data[0])
                plot(np.array([t[0]  for t in iter_size_speedup]), np.array([t[1]  for t in iter_size_speedup]), label='n='+str(n))
            legend()

            figure(figsize=(8,14))

            subplot('311')
            boxplot(speedup_list,positions=nthread_list)
            ylabel('speedup distribution')
            legend()


            subplot('312')
            plot(nthread_list,speedup_avg)
            ylabel('avg. speed up')
            legend()

            subplot('313')
            bar(nthread_list,count_failed)
            ylabel('# fails')
            xlabel('number of threads')
            legend()

            figure(figsize=(16,14))

            for filename,solver in results_filename.items():
                data_tuples = []
                for s in solver:
                    data_tuples.append((s[1],s[2],s[3],s[4]))
                data_tuples=sorted(data_tuples, key=lambda data: data[0])
                try:
                    subplot('211')
                    #plot(np.array([data[0] for data in data_tuples])[:],np.array([data[1] for data in data_tuples])[:], label =filename)
                    plot(np.array([data[0] for data in data_tuples])[:],np.array([data[1] for data in data_tuples])[:])
                    ylabel('cpu time for each problems')
                    xlabel('number of threads')
                    legend()


                    subplot('212')
                    #plot(np.array([data[0] for data in data_tuples])[:],np.array([data[1] for data in data_tuples])[:], label =filename)
                    plot(np.array([data[0] for data in data_tuples])[:],np.array([data_tuples[1][1]/data[1] for data in data_tuples])[:])
                    ylabel('speedup for each problems')
                    xlabel('number of threads')
                    legend()
                except:
                    pass
