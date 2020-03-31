#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=240:0:0 # Request 24 hour runtime
#$ -l h_vmem=1G   # Request 12GB RAM
#$ -m ea
#$ -M m.j.haughey@qmul.ac.uk
#$ -t 1-1
#$ -N dataProcessing

#bash run_codes.sh $SGE_TASK_ID $1 $SGE_TASK_LAST

bash concatenateFullPopulationData.sh $1

bash concatenateLineageData.sh $1
