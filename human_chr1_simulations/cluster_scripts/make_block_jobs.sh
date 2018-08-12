#! /bin/bash

for model in "0" "2" "4"; do
    for h in "0.0" "0.5" "s" ; do

    cd /u/home/b/bkim331/project-klohmueldata/bernard_data/admixture_simulations_fullsims/sims_model${model}/h_${h}
    cwd=$( pwd )

cat << EOF > blockjob.sh
#! /bin/bash
#$ -cwd
#$ -A bkim331
#$ -l h_rt=24:00:00,h_data=4G,arch=intel*
#$ -N b.m${model}.h${h}
#$ -o block.output
#$ -j y
#$ -M bkim331@gmail.com
#$ -m a

cd ${cwd}
. /u/home/b/bkim331/miniconda2/etc/profile.d/conda.sh
conda activate msprime
for i in {1..100}; do
    python /u/home/b/bkim331/project-klohmueldata/bernard_data/admixture_simulations_fullsims/blocklengthdistribution.py h${h}_full_\${i}_p1.ancestry
done
EOF

qsub blockjob.sh

    done
done
