rule yak_triobin:
    input:
        fasta = f"{output_dir}/{trio_id}/fasta/{child}/{{movie}}.fasta",
        pat_yak = f"{output_dir}/{trio_id}/yak/{config[trio_id]['father']['id']}.yak",
        mat_yak = f"{output_dir}/{trio_id}/yak/{config[trio_id]['mother']['id']}.yak"
    output: f"{output_dir}/{trio_id}/yak/{child}.{{movie}}.triobin.txt"
    log: f"{output_dir}/{trio_id}/logs/yak_triobin/{child}.{{movie}}.log"
    conda: "envs/yak.yaml"
    params: extra = hifiasm_param
    threads: 16
    shell: "yak triobin {params.extra} -t {threads} {input.pat_yak} {input.mat_yak} {input.fasta} > {output} 2> {log}"
