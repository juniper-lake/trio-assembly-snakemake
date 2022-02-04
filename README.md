# Trio Assembly Snakemake Workflow

## Input

- Child reads
  - HiFi in `.bam` or `.fastq.gz`
- Parent reads
  - HiFi in `.bam` or `.fastq.gz`
  - Paired end reads in `.bam`, `.cram`, or (`.R1.fasq.gz` + `.R2.fasq.gz`)

### Directory structure within basedir

```text
.
├── cluster_logs # slurm stderr/stdout logs
├── config.yaml # customized copy of example_config.yaml
├── output_trio_assembly # output directory
│   └── <trio_id>
│       ├── fasta
│       ├── hifiasm
│       ├── logs
│       └── yak
└── workflow_trio_assembly # clone of this repo
    ├── calN50
    └── rules
        └── envs
```

## Run the pipeline

```bash
# clone workflow from github into directory named `workflow_trio_assembly`
git clone --recursive git@github.com:juniper-lake/trio-assembly-snakemake.git workflow_trio_assembly

# create directory for cluster logs to be stored.
mkdir cluster_logs

# create `config.yaml` based on `example_config.yaml`.
cp workflow_trio_assembly/example_config.yaml config.yaml

# adjust paths in Snakefile and sbatch script as necessary

# create conda environment to run snakemake workflow.
conda create --prefix ./conda_env --channel bioconda --channel conda-forge lockfile==0.12.2 python=3 snakemake mamba

# activate conda environment. **This conda env must be activated each time you run the workflow.**
conda activate ./conda_env

# run workflow by submitting sbatch script with <trio_id>.
sbatch workflow_trio_assembly/run_snakemake.sh <trio_id>
```

## Configuration

See `example_config.yaml` for an example configuration file.

Trio relationships and input files are defined in `config.yaml` as follows:

```bash
<trio_id>:
  child: 
    id: <child_id>
    reads: 
      - <reads_file>
      - <reads_file>
      - <reads_file>
  father:
    id: <father_id>
    reads:
      - <reads_file>
      - <reads_file>
  mother:
    id: <mother_id>
    reads:
      - <reads_file>
```

## To Do

- find out if yak requires R1 and R2 separate for paired end
  - if not: remove sort step and remove step separating reads into separate files
