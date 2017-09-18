#measure_name=time
#make -f ~/Work/faf/scripts/postprocess.Makefile -k  measure_name=flpops
#make -f ~/Work/faf/scripts/postprocess.Makefile -k  measure_name=time


faf_dir=$(HOME)/Work/faf
domain=python ~/Work/faf/scripts/domains.py
comp=python ~/Work/faf/src/comp.py

all: vi nsgs_localtol_ac_gp nsgs_localtol_vi  nsgs_localsolver nsgs_localsolver_hybrid nsgs_shuffled psor_solvers nsn_solvers prox_solvers prox_series regul_series opti_solvers comp_solvers comp_solvers_large
#all:  nsgs_localsolver nsgs_shuffled psor_solvers nsn_solvers prox_solvers opti_solvers comp_solvers comp_solvers_large


# VI solvers Update rule
vi : 	save_dir=figure/VI/UpdateRule/${measure_name}
vi : 	domain_max=$(shell $(domain) --target=vi --domain)
vi : 	precision=$(shell $(domain) --target=vi --precision)
vi : 	test_name=$(shell $(domain) --test_name)
vi :	
	echo "VI solvers Update rule"
	$(comp) --measure=${measure_name} --display --no-matplot --solvers=FixedPoint-VI-,ExtraGrad-VI --domain=1.0:$(precision):${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-$(test_name).tex;
	pdflatex -interaction=batchmode  profile-$(test_name)_legend.tex;
	echo ${save_dir}
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS Local tolerances
nsgs_localtol_ac_gp :	save_dir=figure/NSGS/LocalTol/${measure_name}
nsgs_localtol_ac_gp : 	domain_max=$(shell $(domain) --target=nsgs_localtol_ac_gp --domain)
nsgs_localtol_ac_gp : 	precision=$(shell $(domain) --target=nsgs_localtol_ac_gp --precision)
nsgs_localtol_ac_gp : 	test_name=$(shell $(domain) --test_name)
nsgs_localtol_ac_gp :
	echo "NSGS Local tolerances"
	$(comp) --measure=${measure_name} --display --no-matplot --solvers=NSGS-AC-GP-1e --domain=1.0:$(precision):${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS Local tolerances
nsgs_localtol_vi :	save_dir=figure/NSGS/LocalTol/VI/${measure_name}
nsgs_localtol_vi : 	domain_max=$(shell $(domain) --target=nsgs_localtol_vi --domain)
nsgs_localtol_vi : 	precision=$(shell $(domain) --target=nsgs_localtol_vi --precision)
nsgs_localtol_vi : 	test_name=$(shell $(domain) --test_name)
nsgs_localtol_vi :
	echo "NSGS Local tolerances"
	$(comp) --measure=${measure_name} --display --no-matplot --solvers=NSGS-PLI-1e --domain=1.0:$(precision):${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS Local Solvers
nsgs_localsolver: save_dir=figure/NSGS/LocalSolver/${measure_name}
nsgs_localsolver : 	domain_max=$(shell $(domain) --target=nsgs_localsolver --domain)
nsgs_localsolver : 	precision=$(shell $(domain) --target=nsgs_localsolver --precision)
nsgs_localsolver : 	test_name=$(shell $(domain) --test_name)
nsgs_localsolver :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC','NSGS-AC-GP','NSGS-AC-100','NSGS-AC-GP-100','NSGS-JM','NSGS-JM-GP','NSGS-PLI-1e-14','NSGS-PLI-1e-06','NSGS-PLI-100','NSGS-PLI-10','NSGS-P','NSGS-Quartic' --domain=1.0:$(precision):${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS Local Solvers
nsgs_localsolver_hybrid : save_dir=figure/NSGS/LocalSolverHybrid/${measure_name}
nsgs_localsolver_hybrid : 	domain_max=$(shell $(domain) --target=nsgs_localsolver --domain)
nsgs_localsolver_hybrid : 	precision=$(shell $(domain) --target=nsgs_localsolver --precision)
nsgs_localsolver_hybrid : 	test_name=$(shell $(domain) --test_name)
nsgs_localsolver_hybrid :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC-GP','NSGS-PLI-100','NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-10-1','NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-10-10','NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-100-1','NSGS-AC-GP-HYBRID-PLI-NSN-LOOP-100-10','NSGS-AC-GP-HYBRID-NSN-PLI-NSN-LOOP-10-1','NSGS-AC-GP-HYBRID-NSN-PLI-NSN-LOOP-10-10','NSGS-AC-GP-HYBRID-NSN-PLI-NSN-LOOP-100-1','NSGS-AC-GP-HYBRID-NSN-PLI-NSN-LOOP-100-10' --domain=1.0:$(precision):${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS  Shuffled
nsgs_shuffled: save_dir=figure/NSGS/Shuffled/${measure_name}
nsgs_shuffled : 	domain_max=$(shell $(domain) --target=nsgs_shuffled --domain)
nsgs_shuffled : 	precision=$(shell $(domain) --target=nsgs_shuffled --precision)
nsgs_shuffled : 	test_name=$(shell $(domain) --test_name)
nsgs_shuffled :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC-GP','NSGS-AC-GP-Shuffled-full','NSGS-AC-GP-Shuffled' --domain=1.0:$(precision):${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}



#PSOR solvers
psor_solvers: save_dir=figure/PSOR/${measure_name}
psor_solvers : 	domain_max=$(shell $(domain) --target=psor_solvers --domain)
psor_solvers : 	precision=$(shell $(domain) --target=psor_solvers --precision)
psor_solvers : 	test_name=$(shell $(domain) --test_name)
psor_solvers :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers='PSOR-'   --domain=1.0:0.1:${domain_max}  --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


#NSN solvers
nsn_solvers: save_dir=figure/NSN/${measure_name}
nsn_solvers : 	domain_max=$(shell $(domain) --target=nsn_solvers --domain)
nsn_solvers : 	precision=$(shell $(domain) --target=nsn_solvers --precision)
nsn_solvers : 	test_name=$(shell $(domain) --test_name)
nsn_solvers :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSN-AlartCurnier','NSN-AlartCurnier-NLS','NSN-AlartCurnier-Generated','NSN-AlartCurnier-Generated-NLS','NSN-AlartCurnier-WRAP-FPP','NSN-AlartCurnier-WRAP-EG','NSN-AlartCurnier-FBLSA','NSN-JeanMoreau','NSN-JeanMoreau-NLS','NSN-JeanMoreau-Generated','NSN-JeanMoreau-Generated-NLS','NSN-JeanMoreau-FBLSA','NSN-FischerBurmeister-GP','NSN-FischerBurmeister-NLS','NSN-FischerBurmeister-FBLSA','NSN-NaturalMap-GP','NSN-NaturalMap-NLS','NSN-NaturalMap-FBLSA','NSN-AlartCurnier-NLS-HYBRID'  --domain=1.0:$(precision):${domain_max}  --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


#PROX solvers
prox_solvers: save_dir=figure/PROX/InternalSolvers/${measure_name}
prox_solvers : 	domain_max=$(shell $(domain) --target=prox_solvers --domain)
prox_solvers : 	precision=$(shell $(domain) --target=prox_solvers --precision)
prox_solvers : 	test_name=$(shell $(domain) --test_name)
prox_solvers :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSN-AlartCurnier','PROX-NSN-AC','PROX-NSN-AC-NLS','PROX-NSN-FB-GP','PROX-NSN-FB-NLS','PROX-NSGS-NSN-AC','PROX-NSN-FB-FBLSA','NSGS-AC'   --domain=1.0:$(precision):${domain_max}  --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

#PROX series
prox_series: save_dir=figure/PROX/Parameters/${measure_name}
prox_series : 	domain_max=$(shell $(domain) --target=prox_series --domain)
prox_series : 	precision=$(shell $(domain) --target=prox_series --precision)
prox_series : 	test_name=$(shell $(domain) --test_name)
prox_series :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers='PROX-NSN-AC-'  --domain=1.0:$(precision):${domain_max}  --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

#Regul series
regul_series: save_dir=figure/REGUL/${measure_name}
regul_series : 	domain_max=$(shell $(domain) --target=regul_series --domain)
regul_series : 	precision=$(shell $(domain) --target=regul_series --precision)
regul_series : 	test_name=$(shell $(domain) --test_name)
regul_series :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='PROX-NSN-AC-regulVar-1e+03','PROX-NSN-AC-regul-1e+04','PROX-NSN-AC-regul-1e+04','PROX-NSN-AC-regul-1e+06','PROX-NSN-AC-regul-1e+08','PROX-NSN-AC-regul-1e+10','PROX-NSN-AC-NLS'  --domain=1.0:$(precision):1000 --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}




#OPTI
opti_solvers: save_dir=figure/OPTI/${measure_name}
opti_solvers : 	domain_max=$(shell $(domain) --target=opti_solvers --domain)
opti_solvers : 	precision=$(shell $(domain) --target=opti_solvers --precision)
opti_solvers : 	test_name=$(shell $(domain) --test_name)
opti_solvers :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-PLI-1e-14','NSGS-PLI-10','TrescaFixedPoint-NSGS-PLI','SOCLCP-NSGS-PLI','ACLMFixedPoint-SOCLCP-NSGS-PLI','PANA-PGS-VI-FPP' --domain=1.0:0.01:${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


#COMP_solvers (best solvers in each families)

comp_solvers: save_dir=figure/COMP/zoom/${measure_name}
comp_solvers : 	domain_max=$(shell $(domain) --target=comp_solvers --domain)
comp_solvers : 	precision=$(shell $(domain) --target=comp_solvers --precision)
comp_solvers : 	test_name=$(shell $(domain) --test_name)
comp_solvers :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC','NSN-AlartCurnier','NSN-AlartCurnier-NLS','PROX-NSN-AC','TrescaFixedPoint-NSGS-PLI','SOCLCP-NSGS-PLI','ACLMFixedPoint-SOCLCP-NSGS-PLI',FixedPoint-VI,ExtraGradient-VI --domain=1.0:0.01:${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

comp_solvers_large: save_dir=figure/COMP/large/${measure_name}
comp_solvers_large : 	domain_max=$(shell $(domain) --target=comp_solvers_large --domain)
comp_solvers_large : 	precision=$(shell $(domain) --target=comp_solvers_large --precision)
comp_solvers_large : 	test_name=$(shell $(domain) --test_name)
comp_solvers_large :
	$(comp) --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC','NSN-AlartCurnier','NSN-AlartCurnier-NLS','PROX-NSN-AC','TrescaFixedPoint-NSGS-PLI','SOCLCP-NSGS-PLI','ACLMFixedPoint-SOCLCP-NSGS-PLI',FixedPoint-VI,ExtraGradient-VI --domain=1.0:$(precision):${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

clean :
	rm -rf ./figure
	rm -f *.tex *.gp *.aux *.dat *.txt *.log

publish:
	cp -r figure $(faf_dir)/TeX
