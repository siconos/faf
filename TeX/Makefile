#measure_name=time
#test_name=LMGC_AqueducPR
#make -f ~/Work/faf/TeX/Makefile -k test_name=LMGC_945_SP_Box_PL measure_name=flpops

all: vi nsgs_localtol nsgs_localsolver nsgs_shuffled nsn_solvers prox_solvers opti_solvers comp_solvers

domain_value :
ifeq ($(test_name),Capsules)
	echo "10" | cat > domain.txt ;
	echo "2" | cat >> domain.txt ;
	echo "4" | cat >> domain.txt
	echo "5" | cat >> domain.txt
	echo "15" | cat >> domain.txt
	echo "100" | cat >> domain.txt
	echo "100" | cat >> domain.txt
	echo "10" | cat >> domain.txt
	echo "20" | cat >> domain.txt
endif
ifeq ($(test_name),LMGC_AqueducPR)
	echo "10" | cat > domain.txt ;
	echo "4" | cat >> domain.txt ;
	echo "5" | cat >> domain.txt
	echo "3" | cat >> domain.txt
	echo "100" | cat >> domain.txt
	echo "100" | cat >> domain.txt
	echo "100" | cat >> domain.txt
	echo "5" | cat >> domain.txt
	echo "5" | cat >> domain.txt
endif
ifeq ($(test_name),LMGC_Cubes_H8_5)
	echo "10" | cat > domain.txt ;
	echo "4" | cat >> domain.txt ;
	echo "5" | cat >> domain.txt
	echo "3" | cat >> domain.txt
	echo "100" | cat >> domain.txt
	echo "100" | cat >> domain.txt
	echo "100" | cat >> domain.txt
	echo "5" | cat >> domain.txt
	echo "5" | cat >> domain.txt
endif
ifeq ($(test_name),LMGC_945_SP_Box_PL)
	echo "10" | cat > domain.txt ;
	echo "4" | cat >> domain.txt ;
	echo "5" | cat >> domain.txt
	echo "3" | cat >> domain.txt
	echo "100" | cat >> domain.txt
	echo "100" | cat >> domain.txt
	echo "100" | cat >> domain.txt
	echo "5" | cat >> domain.txt
	echo "5" | cat >> domain.txt
endif
ifeq ($(test_name),)
echo "enter test_name=[name of the test]"
endif

#export list="10:4:3" 
#export sep=":"
#vi :    domain_max=$(shell IFS=":"; set - $$list; shift 1; echo $1)



# VI solvers Update rule
vi : 	save_dir=figure/VI/UpdateRule/${measure_name}
vi : 	domain_max=$(shell sed '1q;d' domain.txt)
vi :	domain_value
	echo "VI solvers Update rule"
	comp.py --measure=${measure_name} --display --no-matplot --solvers=FixedPoint-VI-,ExtraGrad-VI --domain=1.0:.01:${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-$(test_name).tex;
	pdflatex -interaction=batchmode  profile-$(test_name)_legend.tex;
	echo ${save_dir}
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS Local tolerances
nsgs_localtol:	save_dir=figure/NSGS/LocalTol/${measure_name}
nsgs_localtol: 	domain_max=$(shell sed '2q;d' domain.txt)
nsgs_localtol:	domain_value
	echo "NSGS Local tolerances"
	comp.py --measure=${measure_name} --display --no-matplot --solvers=NSGS-AC-GP-1e,NSGS-AC-GP-0,NSGS-PLI-1e,NSGS-PLI-0 --domain=1.0:.01:${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


# NSGS Local Solvers
nsgs_localsolver: save_dir=figure/NSGS/LocalSolver/${measure_name}
nsgs_localsolver: domain_max=$(shell sed '3q;d' domain.txt)
nsgs_localsolver: domain_value
	comp.py --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC','NSGS-AC-GP','NSGS-JM','NSGS-JM-GP','NSGS-PLI','NSGS-PLI-10','NSGS-P' --domain=1.0:.01:${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS Local Shuffled
nsgs_shuffled: save_dir=figure/NSGS/Shuffled/${measure_name}
nsgs_shuffled: domain_max=$(shell sed '4q;d' domain.txt)
nsgs_shuffled: domain_value
	comp.py --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC-GP','NSGS-AC-GP-Shuffled-full','NSGS-AC-GP-Shuffled' --domain=1.0:.01:${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


#NSN solvers
nsn_solvers: save_dir=figure/NSN/${measure_name}
nsn_solvers: domain_max=$(shell sed '5q;d' domain.txt)
nsn_solvers: domain_value
	comp.py --measure=${measure_name} --display --no-matplot --solvers-exact='NSN-AlartCurnier','NSN-AlartCurnier-NLS','NSN-AlartCurnier-Generated','NSN-AlartCurnier-Generated-NLS','NSN-AlartCurnier-R','NSN-JeanMoreau','NSN-JeanMoreau-NLS','NSN-JeanMoreau-Generated','NSN-JeanMoreau-Generated-NLS','NSN-FischerBurmeister-GP','NSN-FischerBurmeister-NLS'  --domain=1.0:1.0:${domain_max}  --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


#PROX solvers
prox_solvers: save_dir=figure/PROX/InternalSolvers/${measure_name}
prox_solvers: domain_max=$(shell sed '6q;d' domain.txt)
prox_solvers: domain_value
	comp.py --measure=${measure_name} --display --no-matplot --solvers-exact='NSN-AlartCurnier','NSN-AlartCurnier-NLS','NSN-AlartCurnier-Generated','NSN-AlartCurnier-Generated-NLS','PROX-NSN-AC','PROX-NSN-FB-GP','PROX-NSN-FB-FBLSA'  --domain=1.0:1.0:${domain_max}  --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

#PROX series
prox_series: save_dir=figure/PROX/Parameters/${measure_name}
prox_series: domain_max=$(shell sed '6q;d' domain.txt)
prox_series: domain_value
	comp.py --measure=${measure_name} --display --no-matplot --solvers='PROX-NSN-AC-'  --domain=1.0:1.0:${domain_max}  --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}




#OPTI
opti_solvers: save_dir=figure/OPTI/${measure_name}
opti_solvers: domain_max=$(shell sed '7q;d' domain.txt)
opti_solvers: domain_value
	comp.py --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-PLI','NSGS-PLI-10','TrescaFixedPoint-NSGS-PLI','SOCLCP-NSGS-PLI','ACLMFixedPoint-SOCLCP-NSGS-PLI' --domain=1.0:0.01:${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


#COMP_solvers (best solvers in each families)

comp_solvers: save_dir=figure/COMP/${measure_name}
comp_solvers: domain_max=$(shell sed '8q;d' domain.txt)

comp_solvers: domain_value
	comp.py --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC','NSN-AlartCurnier','NSN-AlartCurnier-NLS','PROX-NSN-AC','TrescaFixedPoint-NSGS-PLI','SOCLCP-NSGS-PLI','ACLMFixedPoint-SOCLCP-NSGS-PLI',FixedPoint-VI,ExtraGradient-VI --domain=1.0:0.01:${domain_max} --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

clean :
	rm -rf ./figure
	rm -f *.tex *.gp *.aux *.dat *.txt *.log

