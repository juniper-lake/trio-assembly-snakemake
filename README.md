# Trio Assembly Snakemake Workflow

## Input

- Child reads
  - HiFi in `.bam` or `.fastq.gz`
- Parent reads
  - HiFi in `.bam` or `.fastq.gz`
  - Paired end reads in `.bam`, `.cram`, or (`.R1.fasq.gz` + `.R2.fasq.gz`)

## Getting started

Clone workflow from github into directory named `workflow_trio_assembly` and move into directory.

```git
git clone --recursive git@github.com:juniper-lake/trio-assembly-snakemake.git workflow_trio_assembly
cd workflow_trio_assembly
```

Create directory for cluster logs to be stored.

```bash
mkdir cluster_logs
```

Create `config.yaml` based on `example_config.yaml`.

```bash
cp workflow_trio_assembly/example_config.yaml config.yaml
```

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

Create conda environment to run snakemake workflow.

```bash
conda create -n snakemake -c bioconda -c conda-forge lockfile==0.12.2 python=3 snakemake mamba
```

Activate conda environment. **This conda env must be activated each time you run the workflow.**

```bash
conda activate snakemake
```

Run workflow by submitting sbatch script with <trio_id>.

```bash
sbatch workflow_trio_assembly/run_snakemake.sh <trio_id>
```

### Directory structure within basedir

The workflow creates the directory `output_trio_assembly` where results are deposited. 

```bash
.
├── cluster_logs
├── output_trio_assembly
│   └── HG002_sr
│       ├── fasta
│       ├── hifiasm
│       ├── logs
│       └── yak
└── workflow_trio_assembly
    ├── calN50
    └── rules
        └── envs
```

## To Do

- find out if yak requires R1 and R2 separate for paired end
  - if not: remove sort step and remove step separating reads into separate files