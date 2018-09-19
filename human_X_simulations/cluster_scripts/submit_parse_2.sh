#! /bin/bash

for model in "model0" "model2" "model4"; do
    for h in "0.5" "0.0" "s"; do

    cd /u/home/b/bkim331/project-klohmueldata/bernard_data/admixture_simulations_x/burnin_${model}/h_${h}
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
python /u/home/b/bkim331/project-klohmueldata/bernard_data/admixture_simulations_x/compilestats_treeseq1.py h${h}_full_\${SGE_TASK_ID}_p1.ancestry

EOF

qsub submitparserjob.sh

    done
done
