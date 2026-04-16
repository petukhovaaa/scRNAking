#!/bin/bash

# ------ SLURM Parameters ------

#SBATCH -p general
#SBATCH --cpus-per-task=3
#SBATCH --mem=500000


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
seurat_objects_dir=$3
input=$4
output=$5

Rscript data_analysis_with_merge.R $path $project $seurat_objects_dir $input $output
