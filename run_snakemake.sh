#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -o cluster_logs/slurm-%x-%j-%N.out

umask 002

WORKFLOW_DIR="workflow"

# execute snakemake
snakemake --reason \
    --rerun-incomplete \
    --nolock \
    --keep-going \
    --printshellcmds \
    --config workflow_dir=${WORKFLOW_DIR} \
    --configfile ${WORKFLOW_DIR}/config.yaml \
    --local-cores 4 \
    --jobs 500 \
    --max-jobs-per-second 1 \
    --use-conda --conda-frontend mamba \
    --use-singularity --singularity-args '--nv ' \
    --latency-wait 120 \
    --cluster-config ${WORKFLOW_DIR}/cluster.yaml \
    --cluster "sbatch --partition={cluster.partition} \
                      --cpus-per-task={threads} \
                      --output={cluster.out} {cluster.extra} " \
    --snakefile ${WORKFLOW_DIR}/Snakefile 
