# Cloud-J

<p>
  <a href="https://github.com/geoschem/hemco/cloud-j"><img src="https://img.shields.io/github/v/release/geoschem/cloud-j?include_prereleases&label=Latest%20Pre-Release"></a>
  <a href="https://github.com/geoschem/cloud-j/releases"><img src="https://img.shields.io/github/v/release/geoschem/cloud-j?label=Latest%20Stable%20Release"></a>
  <a href="https://github.com/geoschem/cloud-j/releases/"><img src="https://img.shields.io/github/release-date/geoschem/cloud-j"></a>
  <br />
  <a href="https://doi.org/10.5281/zenodo.13862693"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.13862693.svg" alt="DOI"></a>
  <a href="https://github.com/geoschem/Cloud-J/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg"></a>
  <a href="https://github.com/geoschem/Cloud-J/actions/workflows/ubuntu.yml"><img src="https://github.com/geoschem/Cloud-J/actions/workflows/ubuntu.yml/badge.svg" alt="Ubuntu"></a>
  <a href="https://github.com/geoschem/Cloud-J/actions/workflows/mac.yml"><img src="https://github.com/geoschem/Cloud-J/actions/workflows/mac.yml/badge.svg" alt="Mac"></a>
  <a href="https://github.com/geoschem/Cloud-J/actions/workflows/windows.yml"><img src="https://github.com/geoschem/Cloud-J/actions/workflows/windows.yml/badge.svg" alt="Windows"></a>
</p>

## Description

Cloud-J is a multi-scattering eight-stream radiative transfer model for solar radiation based on Fast-J. It was originally developed by Michael J. Prather (UCI). For information about the origins and history of Cloud-J and its predecessor Fast-J please see the [history document](https://github.com/geoschem/cloud-j/blob/main/docs/History_of_Fast-J_photolysis_code.md) in the docs subdirectory of this repository.

## How to download the model

Download Cloud-J by cloning from GitHub. The default branch when downloading is <tt>main</tt>. You may pass the name of the directory you wish to use.

```
git clone https://github.com/geoschem/cloud-j.git code.cloudj
```

## Setting your environment

Cloud-J can be built using CMake and a fortran compiler. Cloud-J has been successfully tested with the following compilers:
* Intel 2019.1.3.304, 2021.5.0, 2021.10.0
* GNU 10.2.0, 12.2.0

We recommend creating an environment script that loads all libraries you will use. You can source this file every time you build and run Cloud-J to ensure you always use the same environment.

You may need to set stack memory limit to <tt>unlimited</tt> to avoid potential segmentation fault if the system default memory limit is too low. You can set this in your environment as follows.

```
ulimit -s unlimited
```

## How to build Cloud-J standalone

You may build the model from anywhere with access to the source code. We recommend building within a dedicated <tt>build</tt> directory to keep the installation files in one place. Below is an example of building from within the Cloud-J code repository.

```
source /path/to/env
cd /path/to/code
mkdir build
cd build
cmake ..
make -j
make install
```

Another option is to build Cloud-J within a run directory. You can symbolically link to your source code and environment file within the run directory for easy access.

```
mkdir cloudj_rundir
cd cloudj_rundir
ln -s /path/to/env cloudj.env
source cloudj.env
ln -s /path/to/code CodeDir
mkdir build
cd build
cmake ../CodeDir
make -j
make install
```

If your build is successful the executable <tt>cloudj_standalone</tt> will be placed in the <tt>build/bin</tt> directory and also copied to <tt>build</tt>.

## How to run Cloud-J standalone

You may run Cloud-J anywhere you have the the executable <tt>cloudj_standalone</tt>. We recommend creating a run directory and copying the executable there.

```
cd cloudj_rundir
cp path/to/cloudj_standalone .
```

Cloud-J needs input data to run and this data is expected in a local subdirectory called <tt>tables</tt>. The Cloud-J repository comes with example tables stored in source code repository subdirectory <tt>tables</tt>. To do a test run with Cloud-J standalone you may either copy these <tt>tables</tt> from the source code repository to your run directory, or symbolically link to them. Here is an example of symbolically linking to them.

```
cd cloudj_rundir
ln -s /path/to/cloud-j/code/tables
```

Execute <tt>cloudj_standalone</tt> to run the model, sending both standard output and error messages to screen and a log file (cloudj.log). The below example runs Cloud-J in the terminal.

```
cd cloudj_rundir
./cloudj_standalone | tee cloudj.log 2>&1
```


## How to run the tests
Cloud-J ensures correctness by testing all changes against known output. This means that anytime that input data
or an algorithm changes the [reference output](test/expected_output/reference_output.txt) will need to be updated.

To run the tests, compile Cloud-J and then run `make test` or `ctest` from within the build directory.

```
make test
```

```
ctest
```

## Debugging

If you wish to build Cloud-J with compiler debug flags on simply run the following command in your build folder prior to the <tt>make</tt> command.

```
cmake . -DCMAKE_BUILD_TYPE=Debug
```

## Support

Cloud-J is supported by Michael Prather (UCI) and the GEOS-Chem Support Team. Please create a GitHub issue for help. For more information see our [Support Guidelines](https://github.com/geoschem/cloud-j/blob/main/SUPPORT.md).

## Contributing

For information about contributing to Cloud-J, please see our [Contributing Guidelines](https://github.com/geoschem/cloud-j/blob/main/CONTRIBUTING.md).
