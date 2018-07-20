#! /bin/bash
#$ -cwd
#$ -A bkim331
#$ -l h_rt=24:00:00,h_data=4G
#$ -N slim_hselfing_full
#$ -o slim_hselfing_full.output
#$ -j y
#$ -M bkim331@gmail.com
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load python
cd /u/home/b/bkim331/project-klohmueldata/bernard_data/admixture_simulations_arabidopsis/burnin_model0/h_selfing

for i in {1..100}; do
    for j in "p1" "p2"; do
        ./parse_admixture.py hselfing_full_${i}_${j}.vcf > hselfing_aprop_${i}_${j}.txt
    done
done