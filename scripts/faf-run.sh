#!/bin/sh -fx

faf_dir=$HOME/faf
faf_src_dir=$faf_dir/src
faf_scripts_dir=$faf_dir/scripts

fclib_library_dir=$HOME/fclib-library

comp=$faf_src_dir/comp.py
.

example_name=${OAR_JOB_ID}_${example}_${precision}

rundir=/nfs_scratch/$USER/faf/$example_name
mkdir -p $rundir
cd $rundir
. $rundir

# 

rsync -av $fclib_library_dir/$example .
for d in $example; do
    cd $d
    $comp --max-problems=$max_problems --no-compute --no-collect # output problems.txt
    cat problems.txt | $preload parallel $comp --timeout=$timeout --precision=$precision $solvers --no-collect $with_mumps --maxiterls=6 '--files={}'
    $comp --just-collect --timeout=$timeout --precision=$precision --with-mumps --maxiterls=6
    cd ..
done
#cat $HOME/faf/$examples/$0 > command
cp $faf_scripts_dir/$example/$0 $rundir/$0
tar zcvf comps-$example_name.tar.gz `find . -name comp.hdf5` $rundir/$0

mkdir -p $faf_scripts_dir/results
mv comps-$example_name.tar.gz $faf_scripts_dir/results
