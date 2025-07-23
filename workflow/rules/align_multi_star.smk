#! /usr/bin/env python
# -*- coding: utf-8 -*-

#from snakemake import WorkflowError
from itertools import chain

# takes a sample id and determines if its layout
# returns read1 and read2 for paired and just fastq for single
def sample_fastq(sampid):
    # get the list of SRA accessions for this sample
    sample_runs = samples.loc[sampid, 'runid']
    
    # check the layout for the runs
    paired_runs = [r for r in sample_runs if runs.loc[r, 'layout'] == 'PAIRED']
    single_runs = [r for r in sample_runs if runs.loc[r, 'layout'] == 'SINGLE']

    if paired_runs and single_runs:
        raise WorkflowError(f'Mixed layouts for sample "{sampid}"')

    if paired_runs:
        return {
            'read1': expand("results/fastq/{s}/{run}_1.fastq.gz", s=sampid, run=paired_runs),
            'read2': expand("results/fastq/{s}/{run}_2.fastq.gz", s=sampid, run=paired_runs)
       }
    elif single_runs:
        return {
            'readS': expand("results/fastq/{s}/{run}.fastq.gz", s=sampid, run=single_runs)
        }
    else:
        raise WorkflowError(f'No read files found for sample "{sampid}"')
    
def readfiles_input(wildcards):
    return list(chain(*sample_fastq(wildcards.sampid).values()))

def readfiles_param(wildcards):
    d = sample_fastq(wildcards.sampid)
    if 'read1' in d and 'read2' in d:
        return f"{','.join(d['read1'])} {','.join(d['read2'])}"
    elif 'readS' in d:
        return ','.join(d['readS'])
    else:
        raise WorkflowError(f'Invalid sample_fastq {d}')
    return

def star_common_args(wildcards):
    return ' '.join(f'--{k} {v}' for k,v in config['align_multi_star'].items())

rule align_multi_star:
    output:
        "results/align_multi/{sampid}/Aligned.out.bam",
        "results/align_multi/{sampid}/ReadsPerGene.out.tab",
        "results/align_multi/{sampid}/Log.final.out",
        "results/align_multi/{sampid}/SJ.out.tab"
    input:
        readfiles_input,
        genomeDir = config['align_multi_star']['genomeDir']
    params:
        readfiles = readfiles_param,
        common_args = star_common_args
    threads: config['alignment_threads']
    conda: "../envs/star.yaml"
    shell:
        """
        STAR\
            --readFilesIn {params.readfiles}\
            --outFileNamePrefix $(dirname {output[0]})/\
            --runThreadN {threads}\
            {params.common_args}
        """

