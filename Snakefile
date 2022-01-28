import re
import yaml
from pathlib import Path
from collections import defaultdict


shell.prefix(f"set -o pipefail; umask 002; export TMPDIR={config['tmpdir']}; export SINGULARITY_TMPDIR={config['tmpdir']}; ")  # set g+w


# child will be provided command line with `--config cohort=$COHORT`
if "trio_id" in config:
    trio_id = config['trio_id']
else: 
    raise ValueError("Error - trio_id not specified in command line")


# father and mother will be provided in config.yaml
child = config[trio_id]['child']
father = config[trio_id]['father']
mother = config[trio_id]['mother']
samples = [child, father, mother]
ref = config['ref']['shortname']
print(f"Processing child {child} with reference {ref}.")


# set directory names 
workflow_dir = 'workflow' if 'workflow_dir' not in config else config['workflow_dir']
output_dir = 'results' if 'output_dir' not in config else config['output_dir']
input_dir = 'reads' if 'input_dir' not in config else config['input_dir']


# yak/hifiasm parameters
yak_param = '' if 'include_singleton_kmers' in config[trio_id] and config[trio_id]['include_singleton_kmers'] else '-b37'
hifiasm_param = config[trio_id]['hifiasm_kmer'] if 'hifiasm_kmer' in config[trio_id] else ''


# scan smrtcells/ready directory for uBAMs or FASTQs that are ready to process
# uBAMs have priority over FASTQs in downstream processes if both are available
ubam_pattern = re.compile(fr'{input_dir}/(?P<sample>[A-Za-z0-9_-]+)/hifi/(?P<movie>.+).bam')
ubam_dict = defaultdict(dict)
fastq_pattern = re.compile(fr'{input_dir}/(?P<sample>[A-Za-z0-9_-]+)/hifi/(?P<movie>.+).fastq.gz')
fastq_dict = defaultdict(dict)

fastq_sr_pattern = re.compile(fr'{input_dir}/(?P<sample>[A-Za-z0-9_-]+)/paired_end/(?P<movie>.+).R(?P<read>.).fastq.gz')
bamcram_sr_pattern = re.compile(fr'{input_dir}/(?P<sample>[A-Za-z0-9_-]+)/paired_end/(?P<movie>.+).(bam|cram)')
sr_dict = defaultdict(list)

for infile in Path(input_dir).glob('**/hifi/*'):
    ubam_match = ubam_pattern.search(str(infile))
    if ubam_match:
        # create a dict-of-dict to link samples to movie context to uBAM filenames
        ubam_dict[ubam_match.group('sample')][ubam_match.group('movie')] = str(infile)
    fastq_match = fastq_pattern.search(str(infile))
    if fastq_match:
        # create a dict-of-dict to link samples to movie context to FASTQ filenames
        fastq_dict[fastq_match.group('sample')][fastq_match.group('movie')] = str(infile)
for infile in Path(input_dir).glob('**/*'):
    sr_fastq_match = fastq_sr_pattern.search(str(infile))
    if sr_fastq_match:
        sr_dict[sr_fastq_match.group('sample')].append(sr_fastq_match.group('movie'))
    sr_bam_match = bamcram_sr_pattern.search(str(infile))
    if sr_bam_match:
        sr_dict[sr_bam_match.group('sample')].append(sr_bam_match.group('movie'))
   
# get list of unique movie names for each sample, ignoring redundancy between ubam and fastq
hifi_dict = {sample:list(set(list(ubam_dict[sample].keys()) + list(fastq_dict[sample].keys()))) for sample in list(ubam_dict.keys()) + list(fastq_dict.keys())}


# rules to include
include: 'rules/trio_hifiasm.smk'
yak_rule = 'rules/yak_hifi.smk' if 'parent_reads' in config[trio_id] and config[trio_id]['parent_reads'] == 'hifi' else 'rules/yak_sr.smk'
include: yak_rule


# build a list of targets
targets = []
# keep yak dbs?
if 'kmers' in config['targets']:
    targets.extend([f"{output_dir}/{trio_id}/yak/{parent}.yak"
        for parent in [father, mother]])
# assembly and stats
if 'assembly' in config['targets']:
    targets.extend([f"{output_dir}/{trio_id}/hifiasm/{child}.asm.dip.{infix}.{suffix}"
        for suffix in ['fasta.gz', 'fasta.stats.txt', 'trioeval.txt']
        for infix in ['hap1.p_ctg', 'hap2.p_ctg']])
# assembly alignments
if 'alignment' in config['targets']:
    targets.extend([f"{output_dir}/{trio_id}/hifiasm/{child}.asm.{ref}.{suffix}"
        for suffix in ['bam', 'bam.bai']])


localrules: all, md5sum


rule all:
    input: targets + [f"{x}.md5" for x in targets]


rule md5sum:
    input: "{prefix}"
    output: "{prefix}.md5"
    message: "Creating md5 checksum for {input}."
    shell: "md5sum {input} > {output}"
