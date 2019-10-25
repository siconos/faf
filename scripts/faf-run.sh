#!/bin/sh -fx

faf_dir=$HOME/faf
faf_src_dir=$faf_dir/src
faf_scripts_dir=$faf_dir/scripts

fclib_library_dir=$HOME/fclib-library
fclib_library_dir=/scratch/vincent/fclib-library

comp=$faf_src_dir/comp.py

echo `pwd`


example_prefix_without_slash=`echo "$example_prefix" | sed -r 's/[\/]+/_/g'`
echo $example_prefix_without_slash

if test -z ${OAR_JOB_ID}
then
example_name=$(date +%F--%T)_${example_prefix_without_slash}_${example}_${precision}_${timeout}
else
example_name=${OAR_JOB_ID}_${example_prefix_without_slash}_${example}_${precision}_${timeout}
fi
#rundir=/bettik/$USER/faf/$example_name
rundir_base='/scratch/vincent/faf/'
rundir=$rundir_base$example_name
mkdir -p $rundir



cd $rundir
echo `pwd`
echo `python --version`
echo $HOSTNAME

#
cp -r $fclib_library_dir/$example_prefix/$example .
for d in $example; do
    cd $d
    $comp $global --max-problems=$max_problems --no-compute --no-collect # output problems.txt
    # if a file comp.hdf5 is present, we will complete it with the new comparisons.
    if [ -f $faf_scripts_dir/$example/comp.hdf5 ]
    then
       cp  $faf_scripts_dir/$example/comp.hdf5 .
    fi
    #cat problems.txt | $comp --timeout=$timeout --precision=$precision $solvers --no-compute --no-collect $with_mumps --maxiterls=6 '--files={}'# dry run
    cat problems.txt | $preload parallel $comp $global --timeout=$timeout --precision=$precision $solvers --no-collect $with_mumps --maxiterls=6 '--files={}'
    $comp $global --just-collect --timeout=$timeout --precision=$precision --with-mumps --maxiterls=6
    cd ..
done
#cat $HOME/faf/$examples/$0 > command
echo `pwd`
cp $faf_scripts_dir/$example_prefix/$example/$0 $rundir/$0
cd $rundir_base
tar zcvf comps-$example_name.tar.gz `find ${example_name} -name comp.hdf5`  ${example_name}/$0 --force-local
mkdir -p $faf_dir/results
mv comps-$example_name.tar.gz $faf_dir/results
