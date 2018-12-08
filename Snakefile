################################################################################
# SETUP
################################################################################
# Modules
from os.path import join
from snakemake.utils import R

# Configuration
OUTPUT_DIR = config['output_dir'] = config.get('output_dir', 'output') # default: output

# Directories
DATA_DIR = 'data'
SRC_DIR = 'src'

LIB_DIR = join(SRC_DIR, 'lib')
RAW_DIR = join(DATA_DIR, 'raw')

# Files
ADDITIONAL_FEATURES_FILE = join(DATA_DIR, 'additional_features.ICGC-BRCA-EU_BRCA_Davies2017.tsv')
COUNTS_BRCA_FILE = join(DATA_DIR, 'counts.ICGC-BRCA-EU_BRCA_22.SBS-96.tsv')
SAMPLE_KEY_FILE = join(DATA_DIR, 'samples.ICGC-BRCA-EU_BRCA_22.tsv')
RAW_COUNTS_FILE = join(DATA_DIR, 'counts.ICGC-BRCA-EU_BRCA_22.Letouze2017.tsv')
PLOT_1 = join(OUTPUT_DIR, 'plot_1.jpg')
PLOT_2 = join(OUTPUT_DIR, 'plot_2.jpg')
PLOT_3 = join(OUTPUT_DIR, 'plot_3.jpg')
PLOT_4 = join(OUTPUT_DIR, 'plot_4.jpg')
TABLE_1 = join(OUTPUT_DIR, 'table_1.jpg')
TABLE_2 = join(OUTPUT_DIR, 'table_2.jpg')

# Scripts
MAINCODE = join(SRC_DIR, 'main_code.R')

################################################################################
# GENRAL RULES
################################################################################
rule all:
    input:
        ADDITIONAL_FEATURES_FILE,
        COUNTS_BRCA_FILE,
        SAMPLE_KEY_FILE,
        RAW_COUNTS_FILE,
        PLOT_1
    #conda:
    #    "environment.yml"

# Download and plot data
rule plot_ROCs:
    input:
        ADDITIONAL_FEATURES_FILE,
        COUNTS_BRCA_FILE,
        SAMPLE_KEY_FILE,
        RAW_COUNTS_FILE
    output:
        PLOT_1,
        PLOT_2,
        PLOT_3,
        PLOT_4,
        TABLE_1,
        TABLE_2
    shell:
        'Rscript {MAINCODE} {input[0]} {input[1]} {input[2]} {input[3]} {output[0]} {output[1]} {output[2]} {output[3]} {output[4]} {output[5]}'



rule download_additional_features:
    params:
        url='https://obj.umiacs.umd.edu/mutation-signature-explorer/publications/Davies2017/processed/additional_features.ICGC-BRCA-EU_BRCA_Davies2017.tsv'
    output:
        ADDITIONAL_FEATURES_FILE
    shell:
        'wget -O {output} {params.url}'

rule download_counts_brca:
    params:
        url='https://obj.umiacs.umd.edu/mutation-signature-explorer/publications/Davies2017/processed/counts/counts.ICGC-BRCA-EU_BRCA_22.SBS-96.tsv'
    output:
        COUNTS_BRCA_FILE
    shell:
        'wget -O {output} {params.url}'

rule download_sample_key:
    params:
        url='https://obj.umiacs.umd.edu/mutation-signature-explorer/publications/Davies2017/processed/samples/samples.ICGC-BRCA-EU_BRCA_22.tsv'
    output:
        SAMPLE_KEY_FILE
    shell:
        'wget -O {output} {params.url}'

rule download_raw_counts:
    params:
        url='https://obj.umiacs.umd.edu/mutation-signature-explorer/publications/Nik-Zainal2016/processed/sv_counts/counts.ICGC-BRCA-EU_BRCA_22.Letouze2017.tsv'
    output:
        RAW_COUNTS_FILE
    shell:
        'wget -O {output} {params.url}'
