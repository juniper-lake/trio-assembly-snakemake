rule samtools_fasta:
    input: f"reads/{proband}/{{movie}}.bam"
    output: f"results/{proband}/fasta/{proband}/{{movie}}.fasta"
    log: f"results/{proband}/logs/{proband}/{{movie}}.fasta.log"
    threads: 4
    conda: "envs/samtools.yaml"
    message: "Executing {rule}: Converting {input} to {output}."
    shell: "(samtools fasta -@ 3 {input} > {output}) > {log} 2>&1"


rule samtools_fasta_sr:
    input: f"reads/{{sample}}/{{sample}}.bam"
    output: 
        r1=f"results/{proband}/fasta/{{sample}}/{{sample}}.R1.fasta",
        r2=f"results/{proband}/fasta/{{sample}}/{{sample}}.R2.fasta"
    log: f"results/{proband}/logs/{{sample}}.fasta.log"
    threads: 4
    conda: "envs/samtools.yaml"
    shell: 
        """
        (samtools collate -u -O {input} \
        | samtools fasta -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n) > {log} 2>&1
        """


rule yak_count_sr:
    input: 
        r1=f"results/{proband}/fasta/{{sample}}/{{sample}}.R1.fasta",
        r2=f"results/{proband}/fasta/{{sample}}/{{sample}}.R2.fasta"
    output: temp(f"results/{proband}/yak/{{sample}}.yak")
    log: f"results/{proband}/logs/yak/{{sample}}.yak.log"
    benchmark: f"results/{proband}/benchmarks/yak/{{sample}}.yak.tsv"
    conda: "envs/yak.yaml"
    params: "-b37"
    threads: 32
    message: "Executing {rule}: Counting k-mers in {input}."
    shell: "(yak count -t {threads} {params} -o {output} <(zcat {input.r1}) <(zcat {input.r2})) > {log} 2>&1"


rule hifiasm_assemble:
    input: 
        fasta = lambda wildcards: expand(f"results/{proband}/fasta/{proband}/{{movie}}.fasta", movie=ubam_dict[proband]),
        parent1_yak = f"results/{proband}/yak/{parents[0]}.yak",
        parent2_yak = f"results/{proband}/yak/{parents[1]}.yak"
    output:
        temp(f"results/{proband}/hifiasm/{{sample}}.asm.dip.hap1.p_ctg.gfa"),
        f"results/{proband}/hifiasm/{{sample}}.asm.dip.hap1.p_ctg.lowQ.bed",
        f"results/{proband}/hifiasm/{{sample}}.asm.dip.hap1.p_ctg.noseq.gfa",
        temp(f"results/{proband}/hifiasm/{{sample}}.asm.dip.hap2.p_ctg.gfa"),
        f"results/{proband}/hifiasm/{{sample}}.asm.dip.hap2.p_ctg.lowQ.bed",
        f"results/{proband}/hifiasm/{{sample}}.asm.dip.hap2.p_ctg.noseq.gfa",
        temp(f"results/{proband}/hifiasm/{{sample}}.asm.dip.p_utg.gfa"),
        temp(f"results/{proband}/hifiasm/{{sample}}.asm.dip.p_utg.lowQ.bed"),
        temp(f"results/{proband}/hifiasm/{{sample}}.asm.dip.p_utg.noseq.gfa"),
        temp(f"results/{proband}/hifiasm/{{sample}}.asm.dip.r_utg.gfa"),
        temp(f"results/{proband}/hifiasm/{{sample}}.asm.dip.r_utg.lowQ.bed"),
        temp(f"results/{proband}/hifiasm/{{sample}}.asm.dip.r_utg.noseq.gfa"),
        temp(f"results/{proband}/hifiasm/{{sample}}.asm.ec.bin"),
        temp(f"results/{proband}/hifiasm/{{sample}}.asm.ovlp.reverse.bin"),
        temp(f"results/{proband}/hifiasm/{{sample}}.asm.ovlp.source.bin")
    log: f"results/{proband}/logs/hifiasm/{{sample}}.hifiasm.log"
    benchmark: f"results/{proband}/benchmarks/hifiasm/{{sample}}.hifiasm.tsv"
    conda: "envs/hifiasm.yaml"
    params: 
        prefix = f"results/{proband}/hifiasm/{{sample}}.asm",
        parent1 = parents[0],
        parent2 = parents[1]
    threads: 48
    message: "Executing {rule}: Assembling sample {wildcards.sample} from {input.fasta} and parental k-mers."
    shell:
            """
            (
                hifiasm -o {params.prefix} -t {threads} \
                    -1 {input.parent1_yak} -2 {input.parent2_yak} {input.fasta} \
                && (echo -e "hap1\t{params.parent1}\nhap2\t{params.parent2}" > results/{proband}/hifiasm/{wildcards.sample}.asm.key.txt) \
            ) > {log} 2>&1
            """
    

rule gfa2fa:
    input: f"results/{proband}/hifiasm/{{sample}}.asm.dip.{{infix}}.gfa"
    output: f"results/{proband}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta"
    log: f"results/{proband}/logs/gfa2fa/{{sample}}.asm.dip.{{infix}}.log"
    benchmark: f"results/{proband}/benchmarks/gfa2fa/{{sample}}.asm.dip.{{infix}}.tsv"
    conda: "envs/gfatools.yaml"
    message: "Executing {rule}: Extracting fasta from assembly {input}."
    shell: "(gfatools gfa2fa {input} > {output}) 2> {log}"


rule bgzip_fasta:
    input: f"results/{proband}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta"
    output: f"results/{proband}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta.gz"
    log: f"results/{proband}/logs/bgzip/{{sample}}.asm.dip.{{infix}}.log"
    benchmark: f"results/{proband}/benchmarks/bgzip/{{sample}}.asm.dip.{{infix}}.tsv"
    threads: 4
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Compressing {input}."
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


rule asm_stats:
    input: f"results/{proband}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta.gz"
    output: f"results/{proband}/hifiasm/{{sample}}.asm.dip.{{infix}}.fasta.stats.txt"
    log: f"results/{proband}/logs/asm_stats/{{sample}}.asm.dip.{{infix}}.fasta.log"
    benchmark: f"results/{proband}/benchmarks/asm_stats/{{sample}}.asm.dip.{{infix}}.fasta.tsv"
    conda: "envs/k8.yaml"
    message: "Executing {rule}: Calculating stats for {input}."
    shell: f"(k8 workflow/scripts/calN50/calN50.js -f {config['ref']['index']} {{input}} > {{output}}) > {{log}} 2>&1"


rule align_hifiasm:
    input:
        target = config['ref']['fasta'],
        query = [f"results/{proband}/hifiasm/{{sample}}.asm.dip.{infix}.fasta.gz"
                 for infix in ["hap1.p_ctg", "hap2.p_ctg"]]
    output: f"results/{proband}/hifiasm/{{sample}}.asm.{ref}.bam"
    log: f"results/{proband}/logs/align_hifiasm/{{sample}}.asm.{ref}.log"
    benchmark: f"results/{proband}/benchmarks/align_hifiasm/{{sample}}.asm.{ref}.tsv"
    params:
        max_chunk = 200000,
        minimap2_args = "-L --secondary=no --eqx -ax asm5",
        minimap2_threads = 10,
        readgroup = f"@RG\\tID:{{sample}}_hifiasm\\tSM:{{sample}}",
        samtools_threads = 3
    threads: 16  # minimap2 + samtools(+1) + 2x awk + seqtk + cat
    conda: "envs/align_hifiasm.yaml"
    message: "Executing {rule}: Aligning {input.query} to {input.target}."
    shell:
        """
        (cat {input.query} \
            | seqtk seq -l {params.max_chunk} - \
            | awk '{{ if ($1 ~ />/) {{ n=$1; i=0; }} else {{ i++; print n "." i; print $0; }} }}' \
            | minimap2 -t {params.minimap2_threads} {params.minimap2_args} \
                -R '{params.readgroup}' {input.target} - \
                | awk '{{ if ($1 !~ /^@/) \
                                {{ Rct=split($1,R,"."); N=R[1]; for(i=2;i<Rct;i++) {{ N=N"."R[i]; }} print $0 "\tTG:Z:" N; }} \
                              else {{ print; }} }}' \
                | samtools sort -@ {params.samtools_threads} > {output}) > {log} 2>&1
        """


rule samtools_index_bam:
    input: f"results/{proband}/hifiasm/{{prefix}}.bam"
    output: f"results/{proband}/hifiasm/{{prefix}}.bam.bai"
    log: f"results/{proband}/logs/samtools/index/{{prefix}}.log"
    benchmark: f"results/{proband}/logs/samtools/index/{{prefix}}.tsv"
    threads: 4
    conda: "envs/samtools.yaml"
    message: "Executing {rule}: Indexing {input}."
    shell: "(samtools index -@ 3 {input}) > {log} 2>&1"