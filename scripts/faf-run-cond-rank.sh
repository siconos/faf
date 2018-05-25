
#!/bin/sh -fx

faf_dir=$HOME/faf
faf_src_dir=$faf_dir/src
faf_scripts_dir=$faf_dir/scripts

fclib_library_dir=$HOME/fclib-library
#fclib_library_dir=/scratch/Vincent/fclib-library

comp=$faf_src_dir/comp.py
.
if test -z ${OAR_JOB_ID}
then
example_name=$(date +%F--%T)_${example}_${precision}_${timeout}
else
example_name=${OAR_JOB_ID}_${example}_${precision}_${timeout}
fi
rundir=/nfs_scratch/$USER/faf/$example_name
#rundir=/scratch/Vincent/faf/$example_name
mkdir -p $rundir



cd $rundir
. $rundir


echo $HOSTNAME

# 
rsync -av $fclib_library_dir/$example .
for d in $example; do
    cd $d
    $comp --max-problems=$max_problems --no-compute --no-collect # output problems.txt
    cat problems.txt | $preload parallel $comp --compute-cond-rank --forced '--files={}'
    cd ..
done
#cat $HOME/faf/$examples/$0 > command

cp $faf_scripts_dir/$example/$0 $rundir/$0
cd ..
tar zcvf comps-$example_name.tar.gz `find ${example_name} -name comp.hdf5`  ${example_name}/$0 --force-local
mkdir -p $faf_dir/results
mv comps-$example_name.tar.gz $faf_dir/results
