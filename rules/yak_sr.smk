ruleorder: seqtk_fastq_to_fasta_sr > samtools_fasta_sr > samtools_sort_sr_bam > samtools_sort_sr_cram

    
rule samtools_sort_sr_cram:
    input: 
        cram = lambda wildcards: sr_dict[wildcards.sample][".cram"][wildcards.prefix],
        ref = config['cram_ref']
    output: f"{output_dir}/{trio_id}/fasta/{{sample}}/{{prefix}}.sorted.bam"
    log: f"{output_dir}/{trio_id}/logs/yak_sr/{{sample}}.{{prefix}}.samtools_sort_sr_cram.log"
    threads: 8
    conda: "envs/samtools.yaml"
    params: 
        view_threads = 3,
        sort_threads = 3
    shell: 
        """
        (samtools view -@ {params.view_threads} -b -T {input.ref} {input.cram} \
            | samtools sort -n -@ {params.sort_threads} -o {output} ) > {log} 2>&1
        """


rule samtools_sort_sr_bam:
    input: lambda wildcards: sr_dict[wildcards.sample][".bam"][wildcards.prefix],
    output: f"{output_dir}/{trio_id}/fasta/{{sample}}/{{prefix}}.sorted.bam"
    log: f"{output_dir}/{trio_id}/logs/yak_sr/{{sample}}.{{prefix}}.samtools_sort_sr_bam.log"
    threads: 4
    conda: "envs/samtools.yaml"
    shell: "(samtools sort -n -@ 3 -o {output} {input}) > {log} 2>&1"


rule seqtk_fastq_to_fasta_sr:
    input: lambda wildcards: sr_dict[wildcards.sample][".fastq.gz"][wildcards.prefix],
    output: f"{output_dir}/{trio_id}/fasta/{{sample}}/{{prefix}}.fasta"
    log: f"{output_dir}/{trio_id}/logs/yak_sr/{{sample}}.{{prefix}}.seqtk_fastq_to_fasta_sr.log"
    conda: "envs/seqtk.yaml"
    shell: "(seqtk seq -A {input} > {output}) > {log} 2>&1"


rule samtools_fasta_sr:
    input: f"{output_dir}/{trio_id}/fasta/{{sample}}/{{prefix}}.sorted.bam"
    output: 
        r1=f"{output_dir}/{trio_id}/fasta/{{sample}}/{{prefix}}_R1.fasta",
        r2=f"{output_dir}/{trio_id}/fasta/{{sample}}/{{prefix}}_R2.fasta"
    log: f"{output_dir}/{trio_id}/logs/yak_sr/{{sample}}.{{prefix}}.samtools_fasta_sr.log"
    threads: 4
    conda: "envs/samtools.yaml"
    shell: "(samtools fasta {input} -@ 3 -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n) > {log} 2>&1"


rule yak_count_sr:
    input: 
        r1=lambda wildcards: expand(f"{output_dir}/{trio_id}/fasta/{wildcards.sample}/{{prefix}}_R1.fasta", prefix=sr_prefixes[wildcards.sample]),
        r2=lambda wildcards: expand(f"{output_dir}/{trio_id}/fasta/{wildcards.sample}/{{prefix}}_R2.fasta", prefix=sr_prefixes[wildcards.sample])
    output: f"{output_dir}/{trio_id}/yak/{{sample}}.yak"
    log: f"{output_dir}/{trio_id}/logs/yak_sr/{{sample}}.yak_count_sr.log"
    conda: "envs/yak.yaml"
    params: "-b37"
    threads: 32
    shell: "(yak count -t {threads} {params} -o {output} <(cat {input.r1} {input.r2}) <(cat {input.r1} {input.r2})) > {log} 2>&1"
