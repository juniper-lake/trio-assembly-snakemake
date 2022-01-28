#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -o cluster_logs/slurm-%x-%j-%N.out

umask 002

TRIO_ID=$1
WORKFLOW_DIR="workflow_trio_assembly"
OUT_DIR="output_trio_assembly"
mkdir -p ${OUT_DIR}/${TRIO_ID}/
LOCKFILE=${OUT_DIR}/${TRIO_ID}/snakemake.lock

# add lockfile to directory to prevent multiple simultaneous jobs
lockfile -r 0 ${LOCKFILE} || exit 1
trap "rm -f ${LOCKFILE}; exit" SIGINT SIGTERM ERR EXIT


# execute snakemake
snakemake --reason \
    --rerun-incomplete \
    --nolock \
    --keep-going \
    --printshellcmds \
    --config workflow_dir=${WORKFLOW_DIR} output_dir=${OUT_DIR} trio_id=${TRIO_ID} \
    --configfile config.yaml \
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
