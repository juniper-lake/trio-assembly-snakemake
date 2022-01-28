rule yak_count:
    input: lambda wildcards: expand(f"{output_dir}/{trio_id}/fasta/{wildcards.sample}/{{prefix}}.fasta", prefix=hifi_dict[wildcards.sample]),
    output: temp(f"{output_dir}/{trio_id}/yak/{{sample}}.yak")
    log: f"{output_dir}/{trio_id}/logs/yak/{{sample}}.yak.log"
    conda: "envs/yak.yaml"
    params: yak_param
    threads: 32
    shell: "(yak count -t {threads} {params} -o {output} {input}) > {log} 2>&1"
