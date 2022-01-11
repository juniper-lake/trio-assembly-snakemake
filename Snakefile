import re
import yaml
from pathlib import Path
from collections import defaultdict


shell.prefix(f"set -o pipefail; umask 002; export TMPDIR={config['tmpdir']}; export SINGULARITY_TMPDIR={config['tmpdir']}; ")  # set g+w


# cohort will be provided at command line with `--config cohort=$COHORT`
proband = config['proband']
parents = [config['parent1'], config['parent2']]
samples = [config['proband'], config['parent1'], config['parent2']]

ref = config['ref']['shortname']
print(f"Processing proband {proband} with reference {ref}.")


# find all trios in cohort
if parents:
    print(f"Parents IDs for {proband}: {parents}")
else:
    print(f"No parents found for {proband}.")

# scan smrtcells/ready directory for uBAMs or FASTQs that are ready to process
# uBAMs have priority over FASTQs in downstream processes if both are available
ubam_pattern = re.compile(r'smrtcells/ready/(?P<sample>[A-Za-z0-9_-]+)/(?P<movie>m\d{5}[Ue]?_\d{6}_\d{6}).(ccs|hifi_reads).bam')
ubam_dict = defaultdict(dict)
for infile in Path('smrtcells/ready').glob('**/*.bam'):
    ubam_match = ubam_pattern.search(str(infile))
    if ubam_match:
        # create a dict-of-dict to link samples to movie context to uBAM filenames
        ubam_dict[ubam_match.group('sample')][ubam_match.group('movie')] = str(infile)


# build a list of targets
targets = []

# assemble with hifiasm
include: 'rules/trio_hifiasm.smk'

# assembly and stats
targets.extend([f"results/{proband}/hifiasm/{proband}.asm.dip.{infix}.{suffix}"
    for suffix in ['fasta.gz', 'fasta.stats.txt']
    for infix in ['hap1.p_ctg', 'hap2.p_ctg']])
# assembly alignments
targets.extend([f"results/{proband}/hifiasm/{proband}.asm.{ref}.{suffix}"
    for suffix in ['bam', 'bam.bai']])


print(targets)


localrules: all, md5sum


rule all:
    input: targets + [f"{x}.md5" for x in targets]


rule md5sum:
    input: "{prefix}"
    output: "{prefix}.md5"
    message: "Creating md5 checksum for {input}."
    shell: "md5sum {input} > {output}"
