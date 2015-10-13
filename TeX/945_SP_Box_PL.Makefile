#test_name=LMGC_Cubes_H8_5
test_name=LMGC_945_SP_Box_PL

#measure_name=time
measure_name=flpops

#export domain="10:4:3" 
#sep=":"

all: vi nsgs_localtol nsgs_localsolver nsgs_shuffled nsn_solvers prox_solvers opti_solvers

# VI solvers Update rule
vi : 	save_dir=figure/VI/UpdateRule/${measure_name}
#vi : 	element=$(shell `echo $$domain ; IFS=":"; set - $$domain; shift 1; echo  $1`)
vi :
	echo "VI solvers Update rule"
#	echo "element=", ${element}	
	comp.py --measure=${measure_name} --display --no-matplot --solvers=FixedPoint-VI-,ExtraGrad-VI --domain=1.0:.01:10 --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-$(test_name).tex;
	pdflatex -interaction=batchmode  profile-$(test_name)_legend.tex;
	echo ${save_dir}
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS Local tolerances
nsgs_localtol: save_dir=figure/NSGS/LocalTol/${measure_name}
nsgs_localtol:
	echo "NSGS Local tolerances"
	comp.py --measure=${measure_name} --display --no-matplot --solvers=NSGS-AC-GP-1e,NSGS-AC-GP-0,NSGS-PLI-1e,NSGS-PLI-0 --domain=1.0:.01:4 --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


# NSGS Local Solvers
nsgs_localsolver: save_dir=figure/NSGS/LocalSolver/${measure_name}
nsgs_localsolver:
	comp.py --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC','NSGS-AC-GP','NSGS-JM','NSGS-JM-GP','NSGS-PLI','NSGS-PLI-10','NSGS-P' --domain=1.0:.01:5 --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}

# NSGS Local Shuffled
nsgs_shuffled: save_dir=figure/NSGS/Shuffled/${measure_name}
nsgs_shuffled:
	comp.py --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC-GP','NSGS-AC-GP-Shuffled-full','NSGS-AC-GP-Shuffled' --domain=1.0:.01:3 --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


#NSN solvers
nsn_solvers: save_dir=figure/NSN/${measure_name}
nsn_solvers:
	comp.py --measure=${measure_name} --display --no-matplot --solvers-exact='NSN-AlartCurnier','NSN-AlartCurnier-NLS','NSN-AlartCurnier-Generated','NSN-AlartCurnier-Generated-NLS','NSN-AlartCurnier-R','NSN-JeanMoreau','NSN-JeanMoreau-NLS','NSN-JeanMoreau-Generated','NSN-JeanMoreau-Generated-NLS','NSN-FischerBurmeister-GP','NSN-FischerBurmeister-NLS'  --domain=1.0:1.0:100 --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


#PROX solvers
prox_solvers: save_dir=figure/PROX/${measure_name}
prox_solvers:
	comp.py --measure=${measure_name} --display --no-matplot --solvers=PROX,'NSN-AlartCurnier'  --domain=1.0:1.0:100 --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}




#OPTI
opti_solvers: save_dir=figure/OPTI/${measure_name}
opti_solvers:
	comp.py --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-PLI','TrescaFixedPoint-NSGS-PLI','SOCLCP-NSGS-PLI','ACLMFixedPoint-SOCLCP-NSGS-PLI' --domain=1.0:0.01:5 --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}


#COMP_solvers (best solvers in each families)

comp_solvers: save_dir=figure/COMP/${measure_name}
comp_solvers:
	comp.py --measure=${measure_name} --display --no-matplot --solvers-exact='NSGS-AC','TrescaFixedPoint-NSGS-PLI','SOCLCP-NSGS-PLI','ACLMFixedPoint-SOCLCP-NSGS-PLI' --domain=1.0:0.01:5 --gnuplot-profile --gnuplot-separate-keys
	gnuplot profile.gp ;
	pdflatex -interaction=batchmode  profile-${test_name}.tex;
	pdflatex -interaction=batchmode  profile-${test_name}_legend.tex;
	mkdir -p ${save_dir}
	mv  profile-${test_name}.pdf ${save_dir}
	mv  profile-${test_name}_legend.pdf ${save_dir}
