# VICES

VICES (Verify Intensity Contamination from Estimated Sources) is a program that
jointly estimates contamination and its sources in genotyping arrays.

Users should follow the following steps to compile VICES 

<<< SEE https://genome.sph.umich.edu/wiki/VICES FOR DOCUMENTATION >>>

## Prerequisites

Automatic installation of VICES requires 
[cget](http://cget.readthedocs.io/en/latest/src/intro.html#installing-cget) 
and cmake v3.2. These prerequisites can be installed as follows:

Ubuntu 16.04
```
sudo apt-get install cmake python-pip python-dev
pip install cget
```
Ubuntu 14.04
```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:george-edison55/cmake-3.x
sudo apt-get update
sudo apt-get install cmake python-pip python-dev
pip install cget
```
MacOS
```
brew install cmake
sudo easy-install pip
pip install --user cget --ignore-installed six
```

## Installation
The easiest way to install VICES is to use cget as follows:
```bash
cget install --prefix <install_prefix> gjmzajac/VICES
```

Alternatively, you can setup a dev environment cmake directly.
```bash
cd vices
cget install -f ./requirements.txt                     # Install dependencies locally.
mkdir build && cd build                                # Create out of source build directory.
cmake -DCMAKE_TOOLCHAIN_FILE=./cget/cget/cget.cmake .. # Configure project with dependency paths.
make                                                   # Build.
```

## Usage
A typical VICES command line for contamination estimation is as follows
```bash
./vices -r reports_list.txt \
        -o contam_estimates.txt
```

Here reports_list.txt is a text file with paths to Illumina report files 
containing array intensity information and contam_estimates.txt is the output.
