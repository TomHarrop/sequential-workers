#!/usr/bin/env python3

from pathlib import Path

gt_pipeline = ('shub://TomHarrop/honeybee-genotype-pipeline:'
               'honeybee_genotype_pipeline_v0.0.11')
samtools = 'shub://TomHarrop/align-utils:samtools_1.9'


rule target:
    input:
        'output/020_filtering/calls.filtered.vcf.gz'


# populations is filtering out some of the SNPs that I want for bayescan /
# plots, use bcftools filter to make a VCF for that
rule filter_vcf:
    input:
        vcf = 'output/010_genotypes/calls.vcf.gz',
    output:
        'output/020_filtering/calls.filtered.vcf'
    params:
        min_maf = 0.05,
        f_missing = 0.2
    log:
        'output/logs/filter_vcf.log'
    singularity:
        samtools
    shell:
        'bcftools view '
        '-v snps '          # SNPs only
        '-m2 -M2 '          # biallelic sites only
        '--min-af {params.min_maf}:nonmajor '
        '--exclude "F_MISSING>{params.f_missing}" '
        '{input.vcf} '
        '> {output} '
        '2> {log}'


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


# generic vcf index
rule generic_index_vcf:
    input:
        Path('{folder}', '{file}.vcf')
    wildcard_constraints:
        file = '(?!populations_sorted).*'   # this is a "not" to exclude files
    output:
        gz = Path('{folder}', '{file}.vcf.gz'),
        tbi = Path('{folder}', '{file}.vcf.gz.tbi')
    singularity:
        samtools
    shell:
        'bgzip -c {input} > {output.gz} '
        '; '
        'tabix -p vcf {output.gz}'
