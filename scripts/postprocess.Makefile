#measure_name=time
#make -f ~/Work/faf/scripts/postprocess.Makefile -k  measure_name=flpops
#make -f ~/Work/faf/scripts/postprocess.Makefile -k  measure_name=time


faf_dir=$(HOME)/Work/faf
domain=python ~/Work/faf/scripts/domains.py
comp=python ~/Work/faf/src/comp.py

all: vi nsgs_localtol_ac_gp nsgs_localtol_vi  nsgs_localsolver nsgs_localsolver_rho nsgs_localsolver_limited nsgs_localsolver_hybrid nsgs_shuffled psor_solvers nsn_solvers prox_solvers prox_series_nu05 prox_series_nu10 prox_series_nu20 regul_series opti_solvers comp_solvers comp_solvers_large prox_solvers_nsgs
#all:  nsgs_localsolver nsgs_shuffled psor_solvers nsn_solvers prox_solvers opti_solvers comp_solvers comp_solvers_large


# VI solvers Update rule
vi : 	domain_max=$(shell $(domain) --target=vi --domain)
vi : 	precision=$(shell $(domain) --target=vi --precision)
vi : 	step=$(shell $(domain) --target=vi --step)
vi : 	timeout=$(shell $(domain) --target=vi --timeout)
vi : 	test_name=$(shell $(domain) --test_name)
vi : 	save_dir=figure/VI/UpdateRule/${precision}/${timeout}/${measure_name}
vi :	
	echo "VI solvers Update rule"
	$(comp) --measure=${measure_name} --display --no-matplot --solvers=FixedPoint-DeSaxce,FixedPoint-VI-,ExtraGrad-VI --domain=1.0:$(step):${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-$(test_name).tex;
	pdflatex -interaction=batchmode  profile-$(test_name)_legend.tex;
	echo ${save_dir}
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS Local tolerances
nsgs_localtol_ac_gp : 	domain_max=$(shell $(domain) --target=nsgs_localtol_ac_gp --domain)
nsgs_localtol_ac_gp : 	precision=$(shell $(domain) --target=nsgs_localtol_ac_gp --precision)
nsgs_localtol_ac_gp : 	timeout=$(shell $(domain) --target=nsgs_localtol_ac_gp --timeout)
nsgs_localtol_ac_gp : 	step=$(shell $(domain) --target=nsgs_localtol_ac_gp --step)
nsgs_localtol_ac_gp : 	test_name=$(shell $(domain) --test_name)
nsgs_localtol_ac_gp :	save_dir=figure/NSGS/LocalTol/${precision}/${timeout}/${measure_name}
nsgs_localtol_ac_gp :
	echo "NSGS Local tolerances"
	$(comp) --measure=${measure_name} --display --no-matplot --solvers='NSGS-AC-GP-1e','NSGS-AC-GP-ADAPTIVE1' --domain=1.0:$(step):${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS Local tolerances
nsgs_localtol_vi : 	domain_max=$(shell $(domain) --target=nsgs_localtol_vi --domain)
nsgs_localtol_vi : 	step=$(shell $(domain) --target=nsgs_localtol_vi --step)
nsgs_localtol_vi : 	precision=$(shell $(domain) --target=nsgs_localtol_vi --precision)
nsgs_localtol_vi : 	timeout=$(shell $(domain) --target=nsgs_localtol_vi --timeout)
nsgs_localtol_vi : 	test_name=$(shell $(domain) --test_name)
nsgs_localtol_vi :	save_dir=figure/NSGS/LocalTol/VI/${precision}/${timeout}/${measure_name}
nsgs_localtol_vi :
	echo "NSGS Local tolerances"
	$(comp) --measure=${measure_name} --display --no-matplot --solvers='NSGS-PLI-1e','NSGS-PLI-ADAPTIVE1' --domain=1.0:$(step):${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS Local Solvers
nsgs_localsolver : 	domain_max=$(shell $(domain) --target=nsgs_localsolver --domain)
nsgs_localsolver : 	step=$(shell $(domain) --target=nsgs_localsolver --step)
nsgs_localsolver : 	precision=$(shell $(domain) --target=nsgs_localsolver --precision)
nsgs_localsolver : 	timeout=$(shell $(domain) --target=nsgs_localsolver --timeout)
nsgs_localsolver : 	test_name=$(shell $(domain) --test_name)
nsgs_localsolver : save_dir=figure/NSGS/LocalSolver/${precision}/${timeout}/${measure_name}
nsgs_localsolver :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact=\
	'NSGS-AC','NSGS-AC-GP',\
	'NSGS-JM','NSGS-JM-GP',\
	'NSGS-PLI-1e-14','NSGS-PLI-1e-06',\
	'NSGS-P','NSGS-Quartic' --domain=1.0:$(step):${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS Local Solvers rho studies
nsgs_localsolver_rho : 	domain_max=$(shell $(domain) --target=nsgs_localsolver_rho --domain)
nsgs_localsolver_rho : 	step=$(shell $(domain) --target=nsgs_localsolver_rho --step)
nsgs_localsolver_rho : 	precision=$(shell $(domain) --target=nsgs_localsolver_rho --precision)
nsgs_localsolver_rho : 	timeout=$(shell $(domain) --target=nsgs_localsolver_rho --timeout)
nsgs_localsolver_rho : 	test_name=$(shell $(domain) --test_name)
nsgs_localsolver_rho :  save_dir=figure/NSGS/rho/${precision}/${timeout}/${measure_name}
nsgs_localsolver_rho :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC',\
	'NSGS-AC-RHO-GIVEN','NSGS-AC-RHO-SPECTRAL-NORM','NSGS-AC-RHO-SPLIT-SPECTRAL-NORM','NSGS-AC-RHO-SPLIT-SPECTRAL-NORM-COND' \
	--domain=1.0:$(step):${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS Local Solvers
nsgs_localsolver_limited : 	domain_max=$(shell $(domain) --target=nsgs_localsolver --domain)
nsgs_localsolver_limited : 	step=$(shell $(domain) --target=nsgs_localsolver --step)
nsgs_localsolver_limited : 	precision=$(shell $(domain) --target=nsgs_localsolver --precision)
nsgs_localsolver_limited : 	timeout=$(shell $(domain) --target=nsgs_localsolver --timeout)
nsgs_localsolver_limited : 	test_name=$(shell $(domain) --test_name)
nsgs_localsolver_limited :      save_dir=figure/NSGS/limited/${precision}/${timeout}/${measure_name}
nsgs_localsolver_limited :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC',\
	'NSGS-AC-100','NSGS-AC-GP-100',\
	'NSGS-PLI-1e-14','NSGS-PLI-1e-06','NSGS-PLI-100','NSGS-PLI-10' \
	 --domain=1.0:$(step):${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


# NSGS Local Solvers
nsgs_localsolver_hybrid : 	domain_max=$(shell $(domain) --target=nsgs_localsolver --domain)
nsgs_localsolver_hybrid : 	step=$(shell $(domain) --target=nsgs_localsolver --step)
nsgs_localsolver_hybrid : 	precision=$(shell $(domain) --target=nsgs_localsolver --precision)
nsgs_localsolver_hybrid : 	timeout=$(shell $(domain) --target=nsgs_localsolver --timeout)
nsgs_localsolver_hybrid : 	test_name=$(shell $(domain) --test_name)
nsgs_localsolver_hybrid : save_dir=figure/NSGS/LocalSolverHybrid/${precision}/${timeout}/${measure_name}
nsgs_localsolver_hybrid :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC-GP','NSGS-PLI-100','NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-10-1','NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-10-10','NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-100-1','NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-100-10','NSGS-AC-GP-HYBRID-NSN-PLI-NSN-LOOP-10-1','NSGS-AC-GP-HYBRID-NSN-PLI-NSN-LOOP-10-10','NSGS-AC-GP-HYBRID-NSN-PLI-NSN-LOOP-100-1','NSGS-AC-GP-HYBRID-NSN-PLI-NSN-LOOP-100-10' --domain=1.0:$(step):${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS  Shuffled
nsgs_shuffled : 	domain_max=$(shell $(domain) --target=nsgs_shuffled --domain)
nsgs_shuffled : 	step=$(shell $(domain) --target=nsgs_shuffled --step)
nsgs_shuffled : 	precision=$(shell $(domain) --target=nsgs_shuffled --precision)
nsgs_shuffled : 	timeout=$(shell $(domain) --target=nsgs_shuffled --timeout)
nsgs_shuffled : 	test_name=$(shell $(domain) --test_name)
nsgs_shuffled: save_dir=figure/NSGS/Shuffled/${precision}/${timeout}/${measure_name}
nsgs_shuffled :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC-GP','NSGS-AC-GP-Shuffled-full','NSGS-AC-GP-Shuffled' --domain=1.0:$(step):${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}



#PSOR solvers
psor_solvers : 	domain_max=$(shell $(domain) --target=psor_solvers --domain)
psor_solvers : 	step=$(shell $(domain) --target=psor_solvers --step)
psor_solvers : 	precision=$(shell $(domain) --target=psor_solvers --precision)
psor_solvers : 	timeout=$(shell $(domain) --target=psor_solvers --timeout)
psor_solvers : 	test_name=$(shell $(domain) --test_name)
psor_solvers: save_dir=figure/PSOR/${precision}/${timeout}/${measure_name}
psor_solvers :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers='PSOR-'   --domain=1.0:${step}:${domain_max}  --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


#NSN solvers
nsn_solvers_all : 	domain_max=$(shell $(domain) --target=nsn_solvers --domain)
nsn_solvers_all : 	step=$(shell $(domain) --target=nsn_solvers --step)
nsn_solvers_all : 	precision=$(shell $(domain) --target=nsn_solvers --precision)
nsn_solvers_all : 	timeout=$(shell $(domain) --target=nsn_solvers --timeout)
nsn_solvers_all : 	test_name=$(shell $(domain) --test_name)
nsn_solvers_all : save_dir=figure/NSN/${precision}/${timeout}/${measure_name}
nsn_solvers_all :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSN-AlartCurnier','NSN-AlartCurnier-NLS','NSN-AlartCurnier-NLS-RHO-GIVEN','NSN-AlartCurnier-Generated','NSN-AlartCurnier-Generated-NLS','NSN-AlartCurnier-WRAP-FPP','NSN-AlartCurnier-WRAP-EG','NSN-AlartCurnier-FBLSA','NSN-JeanMoreau','NSN-JeanMoreau-NLS','NSN-JeanMoreau-Generated','NSN-JeanMoreau-Generated-NLS','NSN-JeanMoreau-FBLSA','NSN-FischerBurmeister-GP','NSN-FischerBurmeister-NLS','NSN-FischerBurmeister-FBLSA','NSN-NaturalMap-GP','NSN-NaturalMap-NLS','NSN-NaturalMap-FBLSA','NSN-AlartCurnier-NLS-HYBRID'  --domain=1.0:$(step):${domain_max}  --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

#NSN solvers
nsn_solvers : 	domain_max=$(shell $(domain) --target=nsn_solvers --domain)
nsn_solvers : 	step=$(shell $(domain) --target=nsn_solvers --step)
nsn_solvers : 	precision=$(shell $(domain) --target=nsn_solvers --precision)
nsn_solvers : 	timeout=$(shell $(domain) --target=nsn_solvers --timeout)
nsn_solvers : 	test_name=$(shell $(domain) --test_name)
nsn_solvers: save_dir=figure/NSN/${precision}/${timeout}/${measure_name}
nsn_solvers :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact=\
	'NSN-AlartCurnier','NSN-AlartCurnier-NLS','NSN-AlartCurnier-NLS-RHO-GIVEN','NSN-AlartCurnier-FBLSA',\
	'NSN-JeanMoreau','NSN-JeanMoreau-NLS','NSN-JeanMoreau-FBLSA',\
	'NSN-FischerBurmeister-GP','NSN-FischerBurmeister-NLS','NSN-FischerBurmeister-FBLSA',\
	'NSN-NaturalMap-GP','NSN-NaturalMap-NLS','NSN-NaturalMap-FBLSA',\
	'NSN-AlartCurnier-NLS-HYBRID'  --domain=1.0:$(step):${domain_max}  --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


#PROX solvers
prox_solvers : 	step=$(shell $(domain) --target=prox_solvers --step)
prox_solvers : 	domain_max=$(shell $(domain) --target=prox_solvers --domain)
prox_solvers : 	precision=$(shell $(domain) --target=prox_solvers --precision)
prox_solvers : 	timeout=$(shell $(domain) --target=prox_solvers --timeout)
prox_solvers : 	test_name=$(shell $(domain) --test_name)
prox_solvers: save_dir=figure/PROX/NSN/InternalSolvers/${precision}/${timeout}/${measure_name}
prox_solvers :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact=\
	'NSN-AlartCurnier','PROX-NSN-AC','PROX-NSN-AC-NLS',\
	'PROX-NSN-FB-GP','PROX-NSN-FB-NLS','PROX-NSN-FB-FBLSA',\
	'PROX-NSN-AC-regulVar-1e+03','PROX-NSN-AC-nu1.0-sigma0.5' \
	--domain=1.0:$(step):${domain_max}  --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

#PROX solvers NSGS
prox_solvers_nsgs : 	step=$(shell $(domain) --target=prox_solvers --step)
prox_solvers_nsgs : 	domain_max=$(shell $(domain) --target=prox_solvers --domain)
prox_solvers_nsgs : 	precision=$(shell $(domain) --target=prox_solvers --precision)
prox_solvers_nsgs : 	timeout=$(shell $(domain) --target=prox_solvers --timeout)
prox_solvers_nsgs : 	test_name=$(shell $(domain) --test_name)
prox_solvers_nsgs: save_dir=figure/PROX/NSGS/InternalSolvers/${precision}/${timeout}/${measure_name}
prox_solvers_nsgs :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC','PROX-NSGS-NSN-AC' \
	--domain=1.0:$(step):${domain_max}  --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


#PROX series
prox_series_nu05 : 	step=$(shell $(domain) --target=prox_series --step)
prox_series_nu05 : 	domain_max=$(shell $(domain) --target=prox_series --domain)
prox_series_nu05 : 	precision=$(shell $(domain) --target=prox_series --precision)
prox_series_nu05 : 	timeout=$(shell $(domain) --target=prox_series --timeout)
prox_series_nu05 : 	test_name=$(shell $(domain) --test_name)
prox_series_nu05: save_dir=figure/PROX/Parameters/nu05/${precision}/${timeout}/${measure_name}

prox_series_nu05 :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers='PROX-NSN-AC-nu0.5'  --domain=1.0:$(step):${domain_max}  --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

#PROX series
prox_series_nu10 : 	step=$(shell $(domain) --target=prox_series --step)
prox_series_nu10 : 	domain_max=$(shell $(domain) --target=prox_series --domain)
prox_series_nu10 : 	precision=$(shell $(domain) --target=prox_series --precision)
prox_series_nu10 : 	timeout=$(shell $(domain) --target=prox_series --timeout)
prox_series_nu10 : 	test_name=$(shell $(domain) --test_name)
prox_series_nu10: save_dir=figure/PROX/Parameters/nu10/${precision}/${timeout}/${measure_name}
prox_series_nu10 :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers='PROX-NSN-AC-nu1.0'  --domain=1.0:$(step):${domain_max}  --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}
#PROX series
prox_series_nu20 :      step=$(shell $(domain) --target=prox_series --step)
prox_series_nu20 : 	domain_max=$(shell $(domain) --target=prox_series --domain)
prox_series_nu20 : 	precision=$(shell $(domain) --target=prox_series --precision)
prox_series_nu20 : 	timeout=$(shell $(domain) --target=prox_series --timeout)
prox_series_nu20 : 	test_name=$(shell $(domain) --test_name)
prox_series_nu20: save_dir=figure/PROX/Parameters/nu20/${precision}/${timeout}/${measure_name}
prox_series_nu20 :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers='PROX-NSN-AC-nu2.0'  --domain=1.0:$(step):${domain_max}  --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

#Regul series
regul_series : 	domain_max=$(shell $(domain) --target=regul_series --domain)
regul_series : 	step=$(shell $(domain) --target=regul_series --step)
regul_series : 	precision=$(shell $(domain) --target=regul_series --precision)
regul_series : 	timeout=$(shell $(domain) --target=regul_series --timeout)
regul_series : 	test_name=$(shell $(domain) --test_name)
regul_series: save_dir=figure/REGUL/${precision}/${timeout}/${measure_name}
regul_series :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact=\
	'PROX-NSN-AC-regul-1e+04','PROX-NSN-AC-regul-1e+04','PROX-NSN-AC-regul-1e+06','PROX-NSN-AC-regul-1e+08','PROX-NSN-AC-regul-1e+10',\
	'PROX-NSN-AC-NLS',\
	'PROX-NSN-AC-regulVar-1e+03','PROX-NSN-AC-regulVar-1e+04','PROX-NSN-AC-regulVar-1e+06' \
	 --domain=1.0:$(step):${domain_max}  --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}




#OPTI
opti_solvers : step=$(shell $(domain) --target=opti_solvers --step)
opti_solvers : 	domain_max=$(shell $(domain) --target=opti_solvers --domain)
opti_solvers : 	precision=$(shell $(domain) --target=opti_solvers --precision)
opti_solvers : 	timeout=$(shell $(domain) --target=opti_solvers --timeout)
opti_solvers : 	test_name=$(shell $(domain) --test_name)
opti_solvers: save_dir=figure/OPTI/${precision}/${timeout}/${measure_name}

opti_solvers :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact=\
	'SOCLCP-NSGS-PLI',\
	'TRESCA-NSGS-PLI','TRESCA-VI-FP','TRESCA-VI-EG','TRESCA-PG',\
	'ACLMFixedPoint-SOCLCP-NSGS-PLI','ACLMFixedPoint-SOCLCP-VI-EG','ACLMFixedPoint-SOCLCP-VI-FPP',\
	'PANA-PGS-VI-FPP','PANA-PGS-VI-EG','PANA-PGS-CONVEXQP-PG','PANA-CONVEXQP-PG' \
	--domain=1.0:${step}:${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


#COMP_solvers (best solvers in each families)

comp_solvers : 	domain_max=$(shell $(domain) --target=comp_solvers --domain)
comp_solvers : step=$(shell $(domain) --target=comp_solvers --step)
comp_solvers : 	precision=$(shell $(domain) --target=comp_solvers --precision)
comp_solvers : 	timeout=$(shell $(domain) --target=comp_solvers --timeout)
comp_solvers : 	test_name=$(shell $(domain) --test_name)
comp_solvers: save_dir=figure/COMP/zoom/${precision}/${timeout}/${measure_name}
comp_solvers :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact=\
	'NSGS-AC','NSN-AlartCurnier','NSN-AlartCurnier-NLS',\
	'PROX-NSN-AC','PROX-NSN-AC-nu2.0-sigma5.0',\
	'TRESCA-NSGS-PLI','ACLMFixedPoint-SOCLCP-NSGS-PLI','FixedPoint-VI','ExtraGradient-VI' \
	--domain=1.0:${step}:${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

comp_solvers_large : 	domain_max=$(shell $(domain) --target=comp_solvers_large --domain)
comp_solvers_large : 	step=$(shell $(domain) --target=comp_solvers_large --step)
comp_solvers_large : 	precision=$(shell $(domain) --target=comp_solvers_large --precision)
comp_solvers_large : 	timeout=$(shell $(domain) --target=comp_solvers_large --timeout)
comp_solvers_large : 	test_name=$(shell $(domain) --test_name)
comp_solvers_large: save_dir=figure/COMP/large/${precision}/${timeout}/${measure_name}
comp_solvers_large :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact=\
	'NSGS-AC','NSGS-PLI-1e-14','NSGS-PLI-1e-06',\
	'NSN-AlartCurnier','NSN-AlartCurnier-NLS',\
	'PROX-NSN-AC','PROX-NSGS-NSN-AC','PROX-NSN-AC-regulVar-1e+03','PROX-NSN-AC-nu1.0-sigma0.5',\
	'TRESCA-NSGS-PLI','ACLMFixedPoint-SOCLCP-NSGS-PLI','ACLMFixedPoint-SOCLCP-VI-EG',\
	'ExtraGradient-VI' \
	--domain=1.0:$(step):${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}



# ADMM solvers Update rule
admm : 	domain_max=$(shell $(domain) --target=admm --domain)
admm : 	precision=$(shell $(domain) --target=admm --precision)
admm : 	step=$(shell $(domain) --target=admm --step)
admm : 	timeout=$(shell $(domain) --target=admm --timeout)
admm : 	test_name=$(shell $(domain) --test_name)
admm : 	save_dir=figure/ADMM/${precision}/${timeout}/${measure_name}
admm :
	echo "ADMM"
	$(comp) --measure=${measure_name} --display --no-matplot --solvers=ADMM-NORM-INF,ADMM-CST,ADMM-BR,ADMM-SBR,ADMM-BR-NO,ADMM-ASYM-CST,ADMM-ASYM-BR,ADMM-ASYM-SBR \
	--domain=1.0:$(step):${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-$(test_name).tex;
	pdflatex -interaction=batchmode  profile-$(test_name)_legend.tex;
	echo ${save_dir}
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# ADMM solvers Update rule
admm_vs_nsgs : 	domain_max=$(shell $(domain) --target=admm --domain)
admm_vs_nsgs : 	precision=$(shell $(domain) --target=admm --precision)
admm_vs_nsgs : 	step=$(shell $(domain) --target=admm --step)
admm_vs_nsgs : 	timeout=$(shell $(domain) --target=admm --timeout)
admm_vs_nsgs : 	test_name=$(shell $(domain) --test_name)
admm_vs_nsgs : 	save_dir=figure/ADMM_NSGS/${precision}/${timeout}/${measure_name}
admm_vs_nsgs :
	echo "ADMM VS. NSGS"
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact=NSGS-AC,ADMM-BR,ADMM-SBR \
	--domain=1.0:$(step):${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-$(test_name).tex;
	pdflatex -interaction=batchmode  profile-$(test_name)_legend.tex;
	echo ${save_dir}
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# ADMM global solvers Update rule
admm_global_convex : 	domain_max=$(shell $(domain) --target=admm_global_convex --domain)
admm_global_convex : 	precision=$(shell $(domain) --target=admm_global_convex --precision)
admm_global_convex : 	step=$(shell $(domain) --target=admm_global_convex --step)
admm_global_convex : 	timeout=$(shell $(domain) --target=admm_global_convex --timeout)
admm_global_convex : 	test_name=$(shell $(domain) --test_name)
admm_global_convex : 	save_dir=figure/ADMM/Global/convex/${precision}/${timeout}/${measure_name}
admm_global_convex :
	echo "ADMM GLOBAL CONVEX"
	$(comp) --global --measure=${measure_name} --display --no-matplot --solvers-exact='ADMM-CST','ADMM-NORM-INF','ADMM-EIG','ADMM-CST-NO','ADMM-NORM-INF-NO','ADMM-EIG-NO','ADMM-BR','ADMM-SBR','ADMM-BR-SCALED','ADMM-BR-NO' \
	--domain=1.0:$(step):${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-$(test_name).tex;
	pdflatex -interaction=batchmode  profile-$(test_name)_legend.tex;
	echo ${save_dir}
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# ADMM global solvers Update rule
admm_global : 	domain_max=$(shell $(domain) --target=admm_global --domain)
admm_global : 	precision=$(shell $(domain) --target=admm_global --precision)
admm_global : 	step=$(shell $(domain) --target=admm_global --step)
admm_global : 	timeout=$(shell $(domain) --target=admm_global --timeout)
admm_global : 	test_name=$(shell $(domain) --test_name)
admm_global : 	save_dir=figure/ADMM/Global/${precision}/${timeout}/${measure_name}
admm_global :
	echo "ADMM GLOBAL"
	$(comp) --global --measure=${measure_name} --display --no-matplot --solvers-exact='ADMM-CST','ADMM-NORM-INF','ADMM-EIG','ADMM-CST-NO','ADMM-NORM-INF-NO','ADMM-EIG-NO','ADMM-BR','ADMM-SBR','ADMM-BR-SCALED','ADMM-BR-NO','ADMM-BR-FH' \
	--domain=1.0:$(step):${domain_max} --gnuplot-output --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-$(test_name).tex;
	pdflatex -interaction=batchmode  profile-$(test_name)_legend.tex;
	echo ${save_dir}
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

clean :
	rm -rf ./figure
	rm -f *.tex *.gp *.aux *.dat *.txt *.log

publish:
	cp -r figure $(faf_dir)/TeX
publish_admm:
	cp -r figure $(faf_dir)/TeX/ADMM/
