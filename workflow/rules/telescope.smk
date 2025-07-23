#! /usr/bin/env python
# -*- coding: utf-8 -*-

rule telescope:
    output:
        "results/telescope/{sampid}/report.tsv",
        "results/telescope/{sampid}/updated.bam",
        "results/telescope/{sampid}/other.bam"
    input:
        mbam = rules.align_multi_star.output[0],
        annot = config['telescope']['annotation']
    log:
        "results/telescope/{sampid}/telescope.log"        
    conda:
        "../envs/telescope.v1.yaml"
    shell:
        '''
        tdir=$(mktemp -d {config[local_tmp]}/{rule}.{wildcards.sampid}.XXXXXX)

        telescope assign\
            --exp_tag inform\
            --theta_prior 200000\
            --max_iter 1000\
            --updated_sam\
            --outdir $tdir\
            {input.mbam}\
            {input.annot}\
            2>&1 | tee {log[0]}

        mkdir -p $(dirname {output[0]})
        mv $tdir/inform-telescope_report.tsv {output[0]}
        mv $tdir/inform-updated.bam {output[1]}
        mv $tdir/inform-other.bam {output[2]}

        chmod 0440 {output}
        rm -rf $tdir
        '''

