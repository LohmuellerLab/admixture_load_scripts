#! /bin/bash
#$ -cwd
#$ -A bkim331
#$ -l highp,h_rt=100:00:00,h_data=32G
#$ -N self_75
#$ -o selfing_75.output
#$ -j y
#$ -M bkim331@gmail.com
#$ -m a
#$ -t 1-100:1
. /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3
cd /u/home/b/bkim331/project-klohmueldata/bernard_data/admixture_simulations_arabidopsis/burnin_model0/h_selfing
burnseed=${SGE_TASK_ID}
/u/home/b/bkim331/project-klohmueldata/bernard_data/SLiM/bin/slim -seed ${burnseed} slim_hselfing_full_split1000_self75.job
