import re
import yaml
from pathlib import Path
from collections import defaultdict
import os


shell.prefix(f"set -o pipefail; umask 002; export TMPDIR={config['tmpdir']}; export SINGULARITY_TMPDIR={config['tmpdir']}; ")  # set g+w


# configuration
# child will be provided command line with `--config cohort=$COHORT`
if "trio_id" in config: trio_id = config['trio_id']
else: raise ValueError("Error - trio_id not specified in command line")
# father and mother will be provided in config.yaml
child = config[trio_id]['child']['id']
ref = config['ref']['shortname']
print(f"Processing child {child} with reference {ref}.")
# yak/hifiasm parameters
yak_param = '' if 'include_singleton_kmers' in config[trio_id] and config[trio_id]['include_singleton_kmers'] else '-b37'
hifiasm_param = config[trio_id]['hifiasm_kmer'] if 'hifiasm_kmer' in config[trio_id] else ''
yak_rule = 'rules/yak_hifi.smk' if 'parent_reads' in config[trio_id] and config[trio_id]['parent_reads'] == 'hifi' else 'rules/yak_sr.smk'
# set directory names 
workflow_dir = 'workflow_trio_assembly' if 'workflow_dir' not in config else config['workflow_dir']
output_dir = 'output_trio_assembly' if 'output_dir' not in config else config['output_dir']


# make dict of input file paths and types
hifi_dict = defaultdict(lambda: defaultdict(dict))
hifi_prefixes = defaultdict(list)
sr_dict = defaultdict(lambda: defaultdict(dict))
sr_prefixes = defaultdict(list)


for reads_path in config[trio_id]['child']['reads']:
        for extension in ['.fastq.gz', '.bam']:
            if reads_path.endswith(extension):
                hifi_dict[child][extension][os.path.basename(reads_path)[:-len(extension)]] = reads_path
                hifi_prefixes[child].append(os.path.basename(reads_path)[:-len(extension)])
for parent in ['mother', 'father']:
    parent_id = config[trio_id][parent]['id']
    for reads_path in config[trio_id][parent]['reads']:
        if 'parent_reads' in config[trio_id] and config[trio_id]['parent_reads'] == 'hifi':
            for extension in ['.fastq.gz', '.bam']:
                if reads_path.endswith(extension):
                    hifi_dict[parent_id][extension][os.path.basename(reads_path)[:-len(extension)]] = reads_path
                    hifi_prefixes[parent_id].append(os.path.basename(reads_path)[:-len(extension)])
        else:
            for extension in ['.fastq.gz', '.bam', '.cram']:
                if reads_path.endswith(extension):
                    sr_dict[parent_id][extension][os.path.basename(reads_path)[:-len(extension)]] = reads_path
                    prefix = os.path.basename(reads_path)[:-(len(extension)+3)] if extension == '.fastq.gz' else os.path.basename(reads_path)[:-len(extension)]
                    sr_prefixes[parent_id].append(prefix)


for i in hifi_prefixes.keys():
    print("HiFi reads:")
    print(f"Sample: {i}")
    print(f"Filename prefixes {hifi_prefixes[i]}")

for i in sr_prefixes.keys():
    print("Paired end reads:")
    print(f"Sample: {i}")
    print(f"Filename prefixes {sr_prefixes[i]}")


# rules to include
include: 'rules/trio_hifiasm.smk'
include: yak_rule


# build a list of targets
targets = []
# keep yak dbs?
if 'kmers' in config['targets']:
    targets.extend([f"{output_dir}/{trio_id}/yak/{parent}.yak"
        for parent in [config[trio_id]['father']['id'], config[trio_id]['mother']['id']]])
# assembly and stats
if 'assembly' in config['targets']:
    targets.extend([f"{output_dir}/{trio_id}/hifiasm/{child}.asm.dip.{infix}.{suffix}"
        for suffix in ['fasta.gz', 'fasta.stats.txt', 'trioeval.txt']
        for infix in ['hap1.p_ctg', 'hap2.p_ctg']])
# assembly alignments
if 'alignment' in config['targets']:
    targets.extend([f"{output_dir}/{trio_id}/hifiasm/{child}.{hap}.{ref}.{suffix}"
        for hap in ['hap1', 'hap2']
        for suffix in ['bam', 'bam.bai']])


localrules: all, md5sum


rule all:
    input: targets + [f"{x}.md5" for x in targets]


rule md5sum:
    input: "{prefix}"
    output: "{prefix}.md5"
    message: "Creating md5 checksum for {input}."
    shell: "md5sum {input} > {output}"
