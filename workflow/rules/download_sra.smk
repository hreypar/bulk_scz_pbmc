#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Download SRA file 
"""
rule prefetch_sra:
    conda: 
        "../envs/sra_data.yaml"
    output:
        temp("results/sra/{sampid}/{runid}.sra")
    log:
        "results/sra/{sampid}/prefetch_{runid}.log"
    shell:
        '''
	prefetch\
            --max-size u\
	    {wildcards.runid}\
	    --output-file {output[0]}\
	    &> {log}
	'''

# the full directory is removed since there are no logs
#fastq-dump adds SRA ID to each read in the file, to avoid it, use –-origfmt
"""
dump SRA to FASTQ
"""
rule sra_to_fastq_paired:
    conda:
        "../envs/sra_data.yaml"
    output:
        temp("results/fastq/{sampid}/{runid}_1.fastq.gz"),
        temp("results/fastq/{sampid}/{runid}_2.fastq.gz")
    input:
        "results/sra/{sampid}/{runid}.sra"
    threads: 
        snakemake.utils.available_cpu_count()
    wildcard_constraints:
        sra_run = lambda wc: runs.loc[wc.runid]['layout'] == "PAIRED"
    log: 
        "results/fastq/{sampid}/fqdump_{runid}.log"
    shell:
        '''
        tdir=$(mktemp -d {config[local_tmp]}/{rule}.{wildcards.runid}.XXXXXX)
        
        parallel-fastq-dump\
            --threads {threads}\
            --tmpdir $tdir\
            --outdir $tdir\
            --split-3\
            --origfmt\
            --sra-id {input}\
            &> {log}

        mkdir -p $(dirname {output[0]})
        pigz -p {threads} -c $tdir/{wildcards.runid}_1.fastq > {output[0]}
        pigz -p {threads} -c $tdir/{wildcards.runid}_2.fastq > {output[1]}

        rm -rf $tdir

        #chmod 0440 {output}
        #chmod 0550 $(dirname {output[0]})
	'''

"""
dump SRA to FASTQ but single 
"""
rule sra_to_fastq_single:
    conda:
        "../envs/sra_data.yaml"
    output:
        temp("results/fastq/{sampid}/{runid}.fastq.gz")
    input:
        "results/sra/{sampid}/{runid}.sra"
    threads:
        snakemake.utils.available_cpu_count()
    wildcard_constraints:
        sra_run = lambda wc: runs.loc[wc.runid]['layout'] == "SINGLE"
    log:
        "results/fastq/{sampid}/fqdump_{runid}.log"
    shell:
        '''
        tdir=$(mktemp -d {config[local_tmp]}/{rule}.{wildcards.runid}.XXXXXX)
        
        parallel-fastq-dump\
            --threads {threads}\
            --tmpdir $tdir\
            --outdir $tdir\
            --origfmt\
            --skip-technical\
            --sra-id {input}\
            &> {log}

        mkdir -p $(dirname {output[0]})
        pigz -p {threads} -c $tdir/{wildcards.runid}.fastq > {output[0]}

        rm -rf $tdir
        '''

