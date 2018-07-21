#! /bin/bash

for model in "0" "2" "4"; do
    for h in "0.0" "0.5" "s" ; do

    cd /u/home/b/bkim331/project-klohmueldata/bernard_data/admixture_simulations_fullsims/burnin_model${model}/h_${h}
    cwd=$( pwd )

cat << EOF > submitparserjob.sh
#! /bin/bash
#$ -cwd
#$ -A bkim331
#$ -l h_rt=24:00:00,h_data=4G,arch=intel*
#$ -N par.m${model}.h${h}
#$ -o parser.output
#$ -j y
#$ -M bkim331@gmail.com
#$ -m a
#$ -t 1-100:1

cd ${cwd}
. /u/home/b/bkim331/miniconda2/etc/profile.d/conda.sh
conda activate msprime
python /u/home/b/bkim331/project-klohmueldata/bernard_data/admixture_simulations_fullsims/parseancestrytrees.py h${h}_full_\${SGE_TASK_ID}_p1.trees

EOF

qsub submitparserjob.sh

    done
done
