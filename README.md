# Cloud-J

Cloud-J is a multi-scattering eight-stream radiative transfer model for solar radiation based on Fast-J. It was originally developed by Michael J. Prather (UCI). For information about the origins and history of Cloud-J and its predecessor Fast-J please see the [history document](https://github.com/geoschem/cloud-j/blob/main/docs/History_of_Fast-J_photolysis_code.md) in the docs subdirectory of this repository.

## How to build and run Cloud-J standalone

Cloud-J can be built using CMake and a fortran compiler. After downloading the model navigate to the main directory and execute the following commands. If your build is successful an executable will be placed in the build/bin directory.

```
cd /path/to/cloud-j
mkdir build
cd build
cmake ..
make -j
make install
```

You may run Cloud-J from wherever you place the executable cloudj_standalone. Cloud-J needs input data to run and this data is expected in a local subdirectory called tables. You may either copy the tables directory from the source code repository or symbolically link to it.

Because running will produce output we recommend that you create a run directory for your run. Here is an example of how to create a run directory and run the model:

```
mkdir run
cd run
cp /path/to/cloud-j/build/bin/cloudj_standalone .
ln -s /path/to/cloud-j/tables
./cloudj_standalone > cloudj.log
```

## Contributing

For information about contributing to Cloud-J, please see our [Contributing Guidelines](https://github.com/geoschem/cloud-j/blob/main/CONTRIBUTING.md).

## Support

Cloud-J is supported by Michael Prather (UCI) and the GEOS-Chem Support Team. For more information see our [Support Guidelines](https://github.com/geoschem/cloud-j/blob/main/SUPPORT.md). 