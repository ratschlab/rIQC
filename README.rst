===
riqc
===

..
   
   .. image:: https://img.shields.io/pypi/v/riqc.svg
           :target: https://pypi.python.org/pypi/riqc

   .. image:: https://readthedocs.org/projects/riqc/badge/?version=latest
           :target: https://riqc.readthedocs.io/en/latest/?badge=latest
           :alt: Documentation Status




Alignment-free RNA degradation tool


* Free software: MIT license


Installation
------------

The simplest is to install via:
    pip install riqc

If you downloaded the code directly, you can also run:
    python setup.py install

Usage
-----

Degradation Detection in Aligned Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example:
    riqc --bam_dir='...' --anno_fn='...' --out_dir='...' --log='...'

Degradation Detection in Non-Aligned Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example:
    riqc --fastq_dir='...' --genome='...' --anno_fn='...'
    --out_dir='...' --out_fn='...' --pickle_all='...' --pickle_filt='...' --log='...'

All Options for Degradation Detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Input options
    + --bam_dir : Directory of bam files where we get a degradation score for every single bam file and one for all (default='-')
    + --bam_fn : Specifies a single bam file as input   (default='-')
    + --fastq_dir : Directory of fastq files that will be analysed as one (so for individual degradation scores one has to give a directory with only one file or file-pair)    (default='-')
    + --cnt_dir : Directory of pre-produced tab delimited count files (mainly for experimental purposes; the preprocessed file contains all necessary information from the annotation file as well as counts for first and last exon), default='-'

    - --anno_fn : Path to annotation file (supported formats: gtf, gff, gff3)   (default='-')
    - --genome : Path to genome file (supported format: fasta)   (default='-')
    - --gene_list : File with gene-names to use     (default='-')


    + --separate_files_ON : Consider all input files individually   (default=False)
    + sparse_bam_ON : Input BAM files are in sparse hdf5 format     (default=False)

* Output Options
    + --out_dir : Directory to store output in  (default='.')
    + --out_fn : Prefix for output files  (default='out')
    + --anno_tmp_fn : Name of file for temporarily storing annotation information (if the name is '' it will automatically be set (in libs/annotation.py) to reflect whether the protein-coding-genes-filter has been used and whether the legacy option was set to True)    (default='')
    + --pickle_all : Name of pickle file for temporarily storing all kmers (if it is None, the name will automatically be set (in libs/kmer.py) to reflect the chosen length of kmers) (default=None)
    + --pickle_filt : Name of pickle file for temporarily storing filtered/cleaned kmers (if it is None, the name will automatically be set (in libs/kmer.py) to reflect the chosen length of kmers)', (default=None)

* General Options
    + --quant : What type of quantification to use (options: rpkm,raw)  (default='raw')
    + --pseudo_count_ON : Add Pseudocounts to ratio (to also consider genes where we have 0 count on at least one end)  (default=False)
    + --length : Only consider the 25% longest (uq), 50% medium (mq), or 25% shortest (lq) genes (default='uq')

    - --score_on_bases_ON : Calculate degradation score not from last and first exon but from a certain amount of bases at beginning and end that can be set via --base_number (normalized length)  (default=False)
    - --base_number : Number of bases at beginning/end for calculating score (only relevant if --score_on_bases_ON is set)  (default=100)

    + --log : Name of log file  (default='out.log')
    + --verbose_ON : Set logger to verbose  (default=False)
    + --plot_ON : Plot figures  (default=False)
    + --fn_sample_ratio : Sample Ratios in relation to yours (default=os.path.join(os.path.realpath(__file__).rsplit('/', 1)[:-1][0], 'data', 'sampleRatios/TCGA_sample_a_ratio_uq.tsv'))
    + --mask_filter : Mask all read-counts below this integer   (default='0')

    - --protein_coding_filter_OFF' : Only consider genes that are protein-coding     (default=True)
    - --length_filter_OFF : Only consider genes of certain length (specified via --length)  (default=True)
    - --save_counts_ON : Store the exon counts in .npy and .tsv files for later use (via --cnt_dir)     (default=False)
    - --legacy : Switch on some legacy behavior (currently only targeting alternate calculation of transcript length in libs/annotation.py) (default=False)

* Kmer Options (only relevant if using fastq as input)
    + --kmer_length : Length of k-mer for alignmentfree counting    (default=27)
    + --reads_kmer : Required active reads per sample or if in [0, 1] then fraction of input reads considered   (default=50000)
    + --step_k : Step-size for k-mer counting   (default=4)


Additional Options for Degradation Compensation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
+ --scale_counts_ON : Scale counts with pre-computed scaling factors for degradation compensation (gi), default=False
+ --scale_factors_dir : Directory of files containing scaling factors   (default='-')

The scaling-factor files can be generated with a command like
    python .../degradation_tool/scalingFactors.py --bam_dir='...' --anno_fn='...' --out_dir='...'
mainly using the same parameters as for the Degradation Detection with additional:
    + --bins : Number of bins for different gene lengths  (default=10)
    + --relative_binning_ON : Have relative (to number of genes) bin boundaries instead of absolute values  (default=False)
    + --average_factors_ON : Compute scaling factors by using average (instead of median) per length bin    (default=False)


