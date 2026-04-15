#!/bin/bash

# ------ SLURM Parameters ------

#SBATCH -p general
#SBATCH --cpus-per-task=3
#SBATCH --mem=350000


# ------ Make sure it's run as a job ------

if [ -z "${SLURM_JOBID}" ]; then
    echo "Error: Must be run as a job"
    exit 1
fi

eval "$(conda shell.bash hook)"
conda init
conda activate r-env

eval $(spack load --sh fftw)

export LD_LIBRARY_PATH=/ref/vslab/software/conda/envs/r-env/lib:/ref/vslab/software/R/4.2/$LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
export LD_PRELOAD=/ref/vslab/software/conda/envs/r-env/lib/libstdc++.so.6
echo "LD_PRELOAD: $LD_PRELOAD"

path=$1
project=$2
input_1=$3
input_2=$4
merged=$5
output=$6

Rscript merge_two_seurat_obj.R $path $project $input_1 $input_2 $merged

Rscript data_analysis.R $path $project $merged $output 
