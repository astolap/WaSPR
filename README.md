
# WaSPR - light field compression

### Table of contents

 1. [Introduction](#introduction)
 2. [Installing](#installing)
 3. [Running the software](#Running)

## Introduction

This is the WaSPR (Warping and Sparse Prediction on Regions) light field compression software. The program is intended for encoding light fields such as those found in the [JPEG Pleno datasets](https://jpeg.org/plenodb/lf/pleno_lf/). 

The program is developed and maintained by [Pekka Astola](http://www.cs.tut.fi/~astolap/).

If you use this program for research, remember to cite the following publication (currently in early access),

```
@ARTICLE{8924721,
author={P. {Astola} and I. {Tabus}},
journal={IEEE Access},
title={Coding of Light Fields Using Disparity Based Sparse Prediction},
year={2019},
volume={},
number={},
keywords={Cameras;Encoding;Codecs;Lenses;Transform coding;Image coding;Light field coding;image coding;plenoptic;multi-view},
doi={10.1109/ACCESS.2019.2957934},
ISSN={2169-3536},
month={},}
```

## Installing and compiling

This software has been developed using Visual Studio on Windows 10. For Visual Studio a solution file is provided. Makefile will be added soon.

The codec relies on external utilities for various coding stages. You will need [Kakadu (JPEG 2000)](https://kakadusoftware.com/downloads/), [HM (HEVC) 16.20](https://hevc.hhi.fraunhofer.de/), and [gzip](https://www.gzip.org/).

You will also need [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), [Armadillo](http://arma.sourceforge.net/), [OpenBLAS 0.3.6](https://www.openblas.net/). On Windows OpenBLAS requires libgcc_s_seh-1.dll, libgfortran-3.dll, and libquadmath-0.dll from [mingw64-w64](http://mingw-w64.org/doku.php).

## Running the software

Download the light field data sets from [JPEG Pleno database](https://jpeg.org/plenodb/lf/pleno_lf/), and use one of the [configuration files](https://github.com/astolap/WaSPR/blob/master/configuration_files) provided. The path to Kakadu requires only the directory where the binaries of the Kakadu utilities are. For HM encoder/decoder and gzip please provide full paths to the binaries (i.e., paths should end with .exe on Windows).

For HDCA Set 2 use the [HM intra config](https://github.com/astolap/WaSPR/blob/master/configuration_files/encoder_intra_main10.cfg), and for the rest use the [HM inter config](https://github.com/astolap/WaSPR/blob/master/configuration_files/encoder_inter.cfg).

The syntax for the encoder is,
> waspr-encoder --input [INPUT DIRECTORY .PPM/.PGM --output [OUTPUT DIRECTORY .LF] --config [JSON CONFIG FILE] --kakadu [KAKADU BINARY DIRECTORY] --TAppEncoder [PATH TO HM ENCODER BINARY] --TAppDecoder [PATH TO HM DECODER BINARY] --HEVCcfg [PATH TO HM .CFG] --gzip-path  [PATH TO GZIP UTILITY BINARY].

The syntax for the decoder is,
> waspr-decoder --input [INPUT .LF] --output [OUTPUT DIRECTORY .PPM/.PGM] --kakadu [KAKADU BINARY DIRECTORY] --TAppDecoder [PATH TO HM DECODER] --gzip-path  [PATH TO GZIP UTILITY].
