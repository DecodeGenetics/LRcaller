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
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-march=native" ../../LRcaller
make
make test # optional
```

The `-DCMAKE_CXX_FLAGS="-march=native"` is not required but highly recommended. Alternatively, you may specify `-DCMAKE_CXX_FLAGS="-march=x86-64-v3"` to produce a more portable binary (only works on GCC≥10 and Clang≥12).

## Usage

```
lrcaller [OPTIONS] "BAMFILE" "VCF_FILE_IN" "VCF_OUT_FILE"
```

For details, see `lrcaller --help`.


## Citation

LRcaller is developed at deCODE Genetics by Bjarni V. Halldorsson, Doruk Beyter, Hannes Eggertson (@hannespetur) and Hannes Hauswedell (@h-2).
Please cite the following research articel when using LRcaller in any academic work:
https://doi.org/10.1038/s41588-021-00865-4
