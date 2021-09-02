rule remove_vcf_header:
    input:
        vcf = f"{DATA}/{{sample}}.annovar.vcf"
    output:
        vcf_no_header = temp(f"{OUTDIR}/{{sample}}.annovar_noheader.vcf")
    threads:
        get_resource('remove_vcf_header', 'threads')
    resources:
        mem=get_resource('remove_vcf_header', 'mem'),
        walltime=get_resource('remove_vcf_header', 'walltime')
    log:
        f"{LOGDIR}/remove_vcf_header/{{sample}}.log"
    priority: 1
    shell: "sed '/^#/d' {input.vcf} > {output.vcf_no_header}"


# TODO ¿Revisar longitud de los sufijos. Cuantos pedazos va a haber para el primer cromosoma?
checkpoint split_chr:
    input:
        vcf_no_header = rules.remove_vcf_header.output.vcf_no_header
    output:
        split_folder = directory(f"{OUTDIR}/split/{{sample}}")
    threads:
        get_resource('split_chr', 'threads')
    resources:
        mem=get_resource('split_chr', 'mem'),
        walltime=get_resource('split_chr', 'walltime')
    params:
        split_size = get_params('split_chr', 'num_lines')
    log:
        f"{LOGDIR}/split_chr/{{sample}}.log"
    priority: 2
    shell: 'mkdir {output.split_folder}; split -l {params.split_size} "{input.vcf_no_header}" "{output.split_folder}/{wildcards.sample}.part-"'


rule paste_header:
    input:
        split_file = f"{OUTDIR}/split/{{sample}}/{{sample}}.part-{{part}}",
    output:
        headed_file = temp(f"{OUTDIR}/headed/{{sample}}/{{sample}}.part-{{part}}")
    threads:
        get_resource('default', 'threads')
    resources:
        mem=get_resource('default', 'mem'),
        walltime=get_resource('default', 'walltime')
    params:
        header = get_params('paste_header', 'header_vcf')
    log:
        f"{LOGDIR}/paste_header/{{sample}}/{{sample}}.part-{{part}}.log"
    priority: 3
    shell: 'cat {params.header} {input.split_file} > {output.headed_file}'


#rule annotate:
#    input:
#        headed_file = rules.paste_header.output.headed_file
#    output:
#        annotate_part = temp(f"{OUTDIR}/annotate/{{sample}}/{{sample}}.part-{{part}}")
#    threads:
#        get_resource('annotate', 'threads')
#    resources:
#        mem=get_resource('annotate', 'mem'),
#        walltime=get_resource('annotate', 'walltime')
#    params:
#    conda:
#        "../envs/vep_annotation.yaml"
#    log:
#        f"{LOGDIR}/annotate/{{sample}}/{{sample}}.part-{{part}}.log"
#    priority: 4
#    shell: "vep --cache -i {input.headed_file} -e -o {output.annotate_part} --force_overwrite \
#    --assembly GRCh37 --af_gnomad --no_stats --tab --offline --verbose --per_gene --canonical"


rule annotate_variants:
    input:
        calls=rules.paste_header.output.headed_file,  # .vcf, .vcf.gz or .bcf
        cache=f"{VEP_CACHE}", # can be omitted if fasta and gff are specified
        # optionally add reference genome fasta
        # fasta="genome.fasta",
        # fai="genome.fasta.fai", # fasta index
        # gff="annotation.gff",
        # csi="annotation.gff.csi", # tabix index
    output:
        calls=temp(f"{OUTDIR}/annotate/{{sample}}/{{sample}}.part-{{part}}"),  # .vcf, .vcf.gz or .bcf
        stats=f"{OUTDIR}/stats/{{sample}}_variants.html"
    threads:
        get_resource('annotate', 'threads')
    resources:
        mem=get_resource('annotate', 'mem'),
        walltime=get_resource('annotate', 'walltime')
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        extra=get_params('annotate','extra')  # optional: extra arguments
    log:
        "logs/vep/annotate.log"
    priority: 4
    wrapper:
        "0.77.0/bio/vep/annotate"


rule remove_txt_header:
    input:
        annotate_part = rules.annotate.output.annotate_part
    output:
        txt_no_header = f"{OUTDIR}/annotate_no_header/{{sample}}/{{sample}}.part-{{part}}"
    threads:
        get_resource('default', 'threads')
    resources:
        mem=get_resource('default', 'mem'),
        walltime=get_resource('default', 'walltime')
    log:
        f"{LOGDIR}/remove_txt_header/{{sample}}/{{sample}}.part-{{part}}.log"
    priority: 10
    shell: "sed '/^#/d' > {output.txt_no_header}"


def aggregate_input(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the scatter step
    '''
    checkpoint_output = checkpoints.split_chr.get(**wildcards).output[0]
    return expand(f"{OUTDIR}/annotate_no_header/{wildcards.sample}/{wildcards.sample}.part-{{i}}",
           i=glob_wildcards(os.path.join(checkpoint_output, f"{wildcards.sample}.part-{{i}}")).i)


rule join_files:
    input:
        aggregate_input
    output:
        combined = f"{OUTDIR}/final/{{sample}}.txt",
    threads:
        get_resource('default', 'threads')
    resources:
        mem=get_resource('default', 'mem'),
        walltime=get_resource('default', 'walltime')
    params:
        header = get_params('join_files', 'header_txt')
    log:
        f"{LOGDIR}/join_files/{{sample}}.log"
    priority: 10
    shell: "cat {params.header} {input} > {output.combined}"

