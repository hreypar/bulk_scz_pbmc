#! /usr/bin/env python
# -*- coding: utf-8 -*-

import csv

localrules: summarize_counts
rule summarize_counts:
    input:
        telescope_tsv = expand('results/telescope/{sampid}/report.tsv', sampid=meta_table.index.to_list()),
        star_tsv = expand('results/align_multi/{sampid}/ReadsPerGene.out.tab', sampid = meta_table.index.to_list())
    output:
        counts_rds = 'results/dataset/counts.rds',
        features_rds = 'results/dataset/features.rds',
    conda: '../envs/scopetools.yaml'
    params:
        removePAR = True,
        useGeneNames = True,    
        gid2gname_rds = 'resources/gencode.v38/metadata.gid_gname.rds'
    script: '../scripts/summarize_dataset.R'

#localrules: sample_table
#rule sample_table:
#    output:
#        'results/dataset/sample_metadata.tsv'
#    run:
#        selcols = ['sm.sample_id', 'sm.run_id', 'sm.layout']
#        selcols += config['sample_table']['selected_columns']
#        samptable = (
#            runtable[selcols]
#                .rename(columns = dict((c,c[3:]) for c in selcols if c.startswith('sm.')))
#                .groupby('sample_id')
#                .agg(lambda x: ','.join(map(str, sorted(set(x)))))
#                
#        )
#        samptable.to_csv(
#            output[0],
#            sep='\t',
#            header=True,
#            index=True,
#            quoting = csv.QUOTE_NONE
#        )
