# Trio Assembly Snakemake Workflow

## Input

- HiFi reads for child in bam format
- Either HiFi or short reads for parents in either bam or cram format 

## Directory structure within basedir
```
.
├── cluster_logs
├── reads
│   ├── HG002
│   │   └── hifi
│   ├── HG003
│   │   ├── hifi
│   │   └── paired_end
│   ├── HG004
│   │   ├── hifi
│   │   └── paired_end
├── results
│   ├── HG002
│   │   ├── fasta
│   │   ├── hifiasm
│   │   ├── logs
│   │   └── yak
│   └── HG002_sr
│       ├── fasta
│       ├── hifiasm
│       ├── logs
│       └── yak
└── workflow_trio_assembly
    ├── calN50
    └── rules
        └── envs
```