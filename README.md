# LRcaller
Program to estimate genotype structural variants from long read data

## Build instructions

Requirements:

  * CMake >= 3.4
  * SeqAn2 (assumed to be installed in `~/devel/seqan`, adapt paths below respectively; the `-develop` branch is required!)
  * OpenMP
  * C++20 capable compiler (might work with older ones)

Download or clone this repository to e.g. `~/devel/lrcaller` and then:

```sh
mkdir -p ~/devel/lrcaller-build/release
cd ~/devel/lrcaller-build/release
cmake -DCMAKE_CXX_COMPILER=/path/to/g++ -DCMAKE_BUILD_TYPE=Release ../../LRcaller
make
```

In deCODE network on RHEL7, use the following paths:

  * `cmake3` instead of `cmake` (`sudo yum install cmake3` to install it)
  * Add `-DCMAKE_CXX_COMPILER=/nfs/odinn/users/hannesha/bin/g++-10`

# Usage:

```
    LRcaller [OPTIONS] "BAMFILE" "VCF_FILE_IN" "VCF_OUT_FILE"

DESCRIPTION
    Genotypes variants using long reads in bam file/file of bam files

REQUIRED ARGUMENTS
    bam_file/file_of_files STRING
    vcf_in STRING
    vcf_out STRING

OPTIONS
    -h, --help
          Display the help message.
    --version-check BOOL
          Turn this option off to disable version update notifications of the application. One of 1, ON, TRUE, T, YES,
          0, OFF, FALSE, F, and NO. Default: 1.
    -fa, --fa STRING
          Fasta file Default: genome.fa.
    -b2, --bam2 STRING
          Second bam file/file of bam files Default: .
    -a, --aligner STRING
          Aligner used bwa/minimap/seqan (default seqan). One of bwa, minimap, and seqan. Default: seqan.
    -gtm, --genotyper STRING
          Genotyper used joint/ad/va/multi (default joint). One of joint, ad, va, and multi. Default: joint.
    -faA, --faA STRING
          Fasta file alt, used when VCF file refers to an alternate fasta file Default: None.
    -c, --cropread
          Use only the substring of the read within window_size from the variant (recommended)
    -v, --verbose
          Verbose output
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
