# RTKLIB TUTORIAL 01

## Installation

* for Windows:

    `https://www.rtklib.com/rtklib.htm`

    DO NOT download the beta version.


* for Linux:

    `https://github.com/tomojitakasu/RTKLIB`

    GUI is not officially supported in Linux. `https://github.com/rtklibexplorer/RTKLIB` can be a good supplement.

    Compile RTKLIB into a static library (See README.md Sec 8). It should be found as `/RTKLIB/src/build/librtklib.a`.


## User Manual

See `/RTKLIB/doc/manual_2.4.2.pdf` or download from `https://www.rtklib.com/rtklib_document.htm`


## Samples

* Download from `https://www.rtklib.com/rtklib_tutorial.htm` (The corresponding tutorial is for Windows)

* Use test data in `/RTKLIB/test/data/rinex`


## Usage of Example Codes

In `/learning_rtklib/build`, compile the example codes
```
cmake ..
make
```

To use them, run
```
cd /learning_rtklib/build
./spp_example.cpp
```