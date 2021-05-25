# LRcaller
Program to estimate genotype structural variants from long read data

## Build instructions

Requirements:

  * CMake >= 3.4
  * SeqAn2 (currently the `-develop` branch is required!)
  * OpenMP
  * C++20 capable compiler (tested with GCC≥10 and Clang≥10)

Download or clone this repository with submodules:

```sh
mkdir -p ~/devel
cd ~/devel/
git clone --recurse-submodules https://github.com/DecodeGenetics/LRcaller
```

And then build it:

```sh
mkdir -p ~/devel/lrcaller-build/release
cd ~/devel/lrcaller-build/release
cmake -DCMAKE_BUILD_TYPE=Release ../../LRcaller
make
make test # optional
```

In deCODE network on RHEL7, use the following paths:

  * `cmake3` instead of `cmake` (`sudo yum install cmake3` to install it)
  * Add `-DCMAKE_CXX_COMPILER=/nfs/odinn/users/hannesha/bin/g++-10`

## Usage

```
lrcaller [OPTIONS] "BAMFILE" "VCF_FILE_IN" "VCF_OUT_FILE"
```

For details, see `lrcaller --help`.


## Citation

LRcaller is developed at deCODE Genetics by Bjarni V. Halldorsson, Doruk Beyter, Hannes Eggertson (@hannespetur) and Hannes Hauswedell (@h-2).
Please cite the following research articel when using LRcaller in any academic work:
https://doi.org/10.1038/s41588-021-00865-4
