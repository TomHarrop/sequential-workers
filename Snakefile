#!/usr/bin/env python3

from pathlib import Path


def resolve_path(x):
    return Path(x).resolve().as_posix()


gt_pipeline = ('shub://TomHarrop/honeybee-genotype-pipeline:'
               'honeybee_genotype_pipeline_v0.0.11')
plink = 'shub://MarissaLL/singularity-containers:plink_1.9'
samtools = ('shub://TomHarrop/align-utils:samtools_1.10'
            '@d52e0b68cf74f659181d95beac31427c1ade7947')
r = 'shub://TomHarrop/r-containers:r_3.6.3'

rule target:
    input:
        'output/040_maf/maf.Rds',
        'output/020_filtering/calls.filtered.stats.txt',
        'output/010_genotypes/calls.stats.txt',
        


rule generate_maf_table:
    input:
        vcf = 'output/020_filtering/calls.filtered.vcf.gz'
    output:
        maf_mat = 'output/040_maf/maf.Rds',
        maf_dt = 'output/040_maf/maf.csv'
    log:
        'output/logs/generate_maf_table.R'
    singularity:
        r
    script:
        'src/generate_maf_table.R'


# prune LD (doesn't work with this format)
rule prune_vcf:
    input:
        vcf = 'output/020_filtering/calls.filtered.vcf.gz',
        prune = 'output/030_pruning/calls.prune.in'
    output:
        vcf = 'output/030_pruning/calls.pruned.vcf'
    log:
        'output/logs/prune_vcf.log'
    singularity:
        samtools
    shell:
        'bcftools view '
        '-i \'ID=@{input.prune}\' '
        '{input.vcf} '
        '> {output.vcf} '
        '2> {log}'


rule list_pruned_snps:
    input:
        vcf = 'output/020_filtering/calls.filtered.vcf.gz'
    output:
        'output/030_pruning/calls.prune.in'
    params:
        vcf = lambda wildcards, input: resolve_path(input.vcf),
        wd = 'output/030_pruning',
        indep = '50 10 0.1'     # 50 kb window, 10 SNPs, r2 < 0.1
    log:
        resolve_path('output/logs/list_pruned_snps.log')
    singularity:
        plink
    shell:
        'cd {params.wd} || exit 1 ; '
        'plink '
        '--vcf {params.vcf} '
        '--double-id '
        '--allow-extra-chr '
        '--set-missing-var-ids @:# '
        '--indep-pairwise {params.indep} '
        '--out calls '
        '&> {log}'


# filter sites
rule filter_vcf:
    input:
        vcf = 'output/tmp/calls.annotated.vcf',
    output:
        temp('output/020_filtering/calls.filtered.vcf')
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
        '--exclude "F_MISSING>{params.f_missing} || FMT/DP>30 || QUAL<30" '
        '{input.vcf} '
        '> {output} '
        '2> {log}'

rule annotate_loci:
    input:
        'output/010_genotypes/calls.vcf.gz'
    output:
        pipe('output/tmp/calls.annotated.vcf')
    log:
        'output/logs/annotate_loci.log'
    singularity:
        samtools
    shell:
        'bcftools annotate '
        '-Ou '
        '--set-id '
        '+\'%CHROM\\_%POS\\_%REF\\_%FIRST_ALT\' '
        '{input} '
        '>> {output} '
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
        tbi = 'output/010_genotypes/calls.vcf.gz.tbi'
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


# generic vcf stats
rule generic_vcf_stats:
    input:
        vcf = Path('{folder}', '{file}.vcf.gz'),
        tbi = Path('{folder}', '{file}.vcf.gz.tbi')
    wildcard_constraints:
        file = '(?!populations_sorted).*'   # this is a "not" to exclude files
    output:
        stats = Path('{folder}', '{file}.stats.txt'),
        plot = directory(Path('{folder}', '{file}.plots'))
    singularity:
        samtools
    shell:
        'bcftools stats '
        '--verbose '
        '-S <( bcftools query -l {input.vcf} ) '
        '{input.vcf} '
        '> {output.stats} ; '
        'plot-vcfstats -s -v -p {output.plot} {output.stats} '


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
