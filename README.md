# LRcaller
Genotypes variants using long reads in bam file/file of bam files.

## Getting started
### Get a static binary
Get the LRcaller binary with
```sh
wget https://github.com/DecodeGenetics/LRcaller/releases/download/v0.1/LRcaller
chmod a+x LRcaller
```

### Compile from source
You need:
* C++ compiler with C++14 support
* SeqAn v2.4
* zlib
* openmp

Make sure you have these dependencies. Then, download or clone this repository. Next, set the SEQAN_LIB variable in the Makefile to the seqan include directory and finally run

```sh
make
```

## Usage
```
    LRcaller [OPTIONS] "BAMFILE" "VCF_FILE_IN" "VCF_OUT_FILE"

DESCRIPTION
    Genotypes variants using long reads in bam file/file of bam files

REQUIRED ARGUMENTS
    BAMFILE_bam_file/file_of_bam_files STRING
    VCF_FILE_IN_-_input_vcf_file STRING
    VCF_FILE_OUT_-_genotyped_vcf_file STRING

OPTIONS
    -h, --help
          Display the help message.
    --version-check BOOL
          Turn this option off to disable version update notifications of the application. One of 1, ON, TRUE, T, YES,
          0, OFF, FALSE, F, and NO. Default: 1.
    -fa, --fa STRING
          Fasta file Default: genome.fa.
    -b2, --bam2 STRING
          Second bam file/file of bam files, use if you want to use multiple basecallings/mappings for the same reads
          Default: .
    -a, --aligner STRING
          Aligner used bwa/minimap/seqan (default seqan). (Using bwa or minimap requires programs to be in path) One
          of bwa, minimap, and seqan. Default: seqan.
    -gtm, --genotyper STRING
          Genotyper used joint/ad/va/multi (default joint). One of joint, ad, va, and multi. Default: joint.
    -faA, --faA STRING
          Fasta file alt, used when VCF file refers to an alternate fasta file
    -v, --verbose
          Verbose output
    -rb, --right_breakpoint
          Genotype right breakpoint, left breakpoint genotyped by default
    -ora, --get_ref_alt
          Output ref and alt allele, for debugging purposes
    -A, --match INTEGER
          Match score for alignment Default: 1.
    -B, --mismatch INTEGER
          Error penalty for alignment. Default: -1.
    -O, --gap_open INTEGER
          Gap open for alignment Default: -1.
    -E, --gap_extend INTEGER
          Gap extend for alignment Default: -1.
    -w, --window_size INTEGER
          Window size Default: 500.
    -ab, --align_bits DOUBLE
          Max bits for alignment Default: 10.
    -ob, --overlap_bits INTEGER
          Score for overlap (log_2) Default: 10.
    -nt, --number_of_threads INTEGER
          Number of threads Default: 1.
    -vw, --var_window INTEGER
          Var window size, look for del/ins inside this window Default: 100.
    -mapQ, --read_mapq INTEGER
          Minimum read mapQ Default: 30.
    -di, --min_del_ins INTEGER
          Minimum del/ins in cigar to consider Default: 6.
    -rt, --ref_thresh_fraction DOUBLE
          Threshold for fraction of del/ins bp to be considered ref Default: 0.1.
    -at, --alt_thresh_fraction DOUBLE
          Threshold for fraction of del/ins bp to be considered alt Default: 0.5.
    -lsf, --logScaleFactor DOUBLE
          Log scale factor for comparing alignment scores Default: 2.
```
