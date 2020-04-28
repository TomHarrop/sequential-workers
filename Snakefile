#!/usr/bin/env python3

gt_pipeline = ('shub://TomHarrop/honeybee-genotype-pipeline:'
               'honeybee_genotype_pipeline_v0.0.11')

rule target:
    input:
        'output/010_genotypes/calls.vcf.gz'

checkpoint genotype:
    input:
        csv = 'data/sample_key.csv',
        ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
    output:
        # disable until the pipeline works again
        # bam = 'output/010_genotypes/{ref}/merged.bam',
        # cutoffs = 'output/010_genotypes/{ref}/040_stats/ldepth.mean_cutoffs.csv',
        # fai = 'output/010_genotypes/{ref}/015_ref/ref.fasta.fai',
        # ref = 'output/010_genotypes/{ref}/015_ref/ref.fasta',
        vcf = 'output/010_genotypes/calls.vcf.gz',
    params:
        wd = 'output/010_genotypes',
        ploidy = '20'
    log:
        'output/logs/genotype.log'
    threads:
        workflow.cores
    singularity:
        gt_pipeline
    shell:
        'honeybee_genotype_pipeline '
        '--ref {input.ref} '
        '--samples_csv {input.csv} '
        '--outdir {params.wd} '
        '--ploidy {params.ploidy} '
        '--threads {threads} '
        '--restart_times 1 '
        '&> {log}'