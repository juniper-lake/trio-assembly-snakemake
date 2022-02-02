ruleorder: samtools_fasta > seqtk_fastq_to_fasta


rule samtools_fasta:
    input: lambda wildcards: ubam_dict[wildcards.sample][wildcards.prefix]
    output: temp(f"{output_dir}/{trio_id}/fasta/{{sample}}/{{prefix}}.fasta")
    log: f"{output_dir}/{trio_id}/logs/trio_hifiasm/{{sample}}.{{prefix}}.samtools_fasta.log"
    threads: 4
    conda: "envs/samtools.yaml"
    shell: "(samtools fasta -@ 3 {input} > {output}) > {log} 2>&1"


rule seqtk_fastq_to_fasta:
    input: lambda wildcards: fastq_dict[wildcards.sample][wildcards.prefix]
    output: temp(f"{output_dir}/{trio_id}/fasta/{{sample}}/{{prefix}}.fasta")
    log: f"{output_dir}/{trio_id}/logs/trio_hifiasm/{{sample}}.{{prefix}}.seqtk_fastq_to_fasta.log"
    conda: "envs/seqtk.yaml"
    shell: "(seqtk seq -A {input} > {output}) > {log} 2>&1"


rule hifiasm:
    input: 
        fasta = lambda wildcards: expand(f"{output_dir}/{trio_id}/fasta/{child}/{{movie}}.fasta", movie=hifi_dict[child]),
        pat_yak = f"{output_dir}/{trio_id}/yak/{father}.yak",
        mat_yak = f"{output_dir}/{trio_id}/yak/{mother}.yak"
    output:
        temp(f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.hap1.p_ctg.gfa"),
        f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.hap1.p_ctg.lowQ.bed",
        f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.hap1.p_ctg.noseq.gfa",
        temp(f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.hap2.p_ctg.gfa"),
        f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.hap2.p_ctg.lowQ.bed",
        f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.hap2.p_ctg.noseq.gfa",
        temp(f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.p_utg.gfa"),
        f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.p_utg.lowQ.bed",
        f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.p_utg.noseq.gfa",
        temp(f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.r_utg.gfa"),
        f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.r_utg.lowQ.bed",
        f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.r_utg.noseq.gfa",
        temp(f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.ec.bin"),
        temp(f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.ovlp.reverse.bin"),
        temp(f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.ovlp.source.bin")
    log: f"{output_dir}/{trio_id}/logs/trio_hifiasm/{{sample}}.hifiasm.log"
    conda: "envs/hifiasm.yaml"
    params: 
        prefix = f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm",
        hap1 = father,
        hap2 = mother,
        extras = hifiasm_param
    threads: 48
    shell:
            """
            (
                hifiasm -o {params.prefix} -t {threads} {params.extras} \
                    -1 {input.pat_yak} -2 {input.mat_yak} {input.fasta} \
                && (echo -e "hap1\t{params.hap1}\nhap2\t{params.hap2}" > {output_dir}/{trio_id}/hifiasm/{wildcards.sample}.asm.key.txt) \
            ) > {log} 2>&1
            """
    

rule gfa2fa:
    input: f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.{{infix}}.gfa"
    output: f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta"
    log: f"{output_dir}/{trio_id}/logs/trio_hifiasm/{{sample}}.{{infix}}.gfa2fa.log"
    conda: "envs/gfatools.yaml"
    message: "Executing {rule}: Extracting fasta from assembly {input}."
    shell: "(gfatools gfa2fa {input} > {output}) 2> {log}"


rule bgzip_fasta:
    input: f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta"
    output: f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta.gz"
    log: f"{output_dir}/{trio_id}/logs/trio_hifiasm/{{sample}}.{{infix}}.bgzip_fasta.log"
    threads: 4
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Compressing {input}."
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


rule asm_stats:
    input: f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta.gz"
    output: f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta.stats.txt"
    log: f"{output_dir}/{trio_id}/logs/trio_hifiasm/{{sample}}.{{infix}}.asm_stats.log"
    conda: "envs/k8.yaml"
    message: "Executing {rule}: Calculating stats for {input}."
    shell: f"(k8 {workflow_dir}/calN50/calN50.js -f {config['ref']['index']} {{input}} > {{output}}) > {{log}} 2>&1"


rule yak_trioeval:
    input: 
        pat_yak = f"{output_dir}/{trio_id}/yak/{father}.yak",
        mat_yak = f"{output_dir}/{trio_id}/yak/{mother}.yak",
        fasta = f"{output_dir}/{trio_id}/hifiasm/{child}.asm.dip.{{infix}}.fasta.gz",
    output: f"{output_dir}/{trio_id}/hifiasm/{child}.asm.dip.{{infix}}.trioeval.txt"
    log: f"{output_dir}/{trio_id}/logs/trio_hifiasm/{child}.{{infix}}.yak_trioeval.log"
    conda: "envs/yak.yaml"
    threads: 16
    shell: "(yak trioeval -t {threads} {input.pat_yak} {input.mat_yak} {input.fasta} > {output}) > {log} 2>&1"


rule align_hifiasm:
    input:
        ref = config['ref']['fasta'],
        assembly = f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.{{hap}}.p_ctg.fasta.gz"
    output: f"{output_dir}/{trio_id}/hifiasm/{{sample}}.{{hap}}.{ref}.bam"
    log: f"{output_dir}/{trio_id}/logs/trio_hifiasm/{{sample}}.{{hap}}.align_hifiasm.log"
    params:
        minimap2_args = "-L --secondary=no --eqx -ax asm5",
        minimap2_threads = 10,
        readgroup = f"@RG\\tID:{{sample}}_{{hap}}\\tSM:{{sample}}",
        samtools_threads = 3,
        samtools_mem = "-m8G"
    threads: 16  # minimap2 + samtools(+1) + 2x awk + seqtk + cat
    conda: "envs/align_hifiasm.yaml"
    shell:
        """
        (minimap2 -t {params.minimap2_threads} {params.minimap2_args} {input.ref} \
                -R '{params.readgroup}' {input.assembly} - \
                | samtools sort -@ {params.samtools_threads} {params.samtools_mem} > {output}) > {log} 2>&1
        """


rule samtools_index_bam:
    input: f"{output_dir}/{trio_id}/hifiasm/{{prefix}}.bam"
    output: f"{output_dir}/{trio_id}/hifiasm/{{prefix}}.bam.bai"
    log: f"{output_dir}/{trio_id}/logs/trio_hifiasm/index/{{prefix}}.samtools_index_bam.log"
    threads: 4
    conda: "envs/samtools.yaml"
    message: "Executing {rule}: Indexing {input}."
    shell: "(samtools index -@ 3 {input}) > {log} 2>&1"


# rule align_hifiasm_chunk:
#     input:
#         target = config['ref']['fasta'],
#         query = [f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.dip.{infix}.fasta.gz"
#                  for infix in ["hap1.p_ctg", "hap2.p_ctg"]]
#     output: f"{output_dir}/{trio_id}/hifiasm/{{sample}}.asm.{ref}.bam"
#     log: f"{output_dir}/{trio_id}/logs/align_hifiasm/{{sample}}.asm.{ref}.log"
#     params:
#         max_chunk = 200000,
#         minimap2_args = "-L --secondary=no --eqx -ax asm5",
#         minimap2_threads = 10,
#         readgroup = f"@RG\\tID:{{sample}}_hifiasm\\tSM:{{sample}}",
#         samtools_threads = 3
#     threads: 16  # minimap2 + samtools(+1) + 2x awk + seqtk + cat
#     conda: "envs/align_hifiasm.yaml"
#     message: "Executing {rule}: Aligning {input.query} to {input.target}."
#     shell:
#         """
#         (cat {input.query} \
#             | seqtk seq -l {params.max_chunk} - \
#             | awk '{{ if ($1 ~ />/) {{ n=$1; i=0; }} else {{ i++; print n "." i; print $0; }} }}' \
#             | minimap2 -t {params.minimap2_threads} {params.minimap2_args} \
#                 -R '{params.readgroup}' {input.target} - \
#                 | awk '{{ if ($1 !~ /^@/) \
#                                 {{ Rct=split($1,R,"."); N=R[1]; for(i=2;i<Rct;i++) {{ N=N"."R[i]; }} print $0 "\tTG:Z:" N; }} \
#                               else {{ print; }} }}' \
#                 | samtools sort -@ {params.samtools_threads} > {output}) > {log} 2>&1
#         """
