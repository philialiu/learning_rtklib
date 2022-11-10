# RTKLIB TUTORIAL 01

## Installation

* for Windows:

    <https://www.rtklib.com/rtklib.htm>

    The beta version is not recommended.


* for Linux:

    <https://github.com/tomojitakasu/RTKLIB>

    GUI is not officially supported in Linux. <https://github.com/rtklibexplorer/RTKLIB> can be a good supplement.

    Create `CMakeLists.txt` under `RTKLIB/src`

    ```
    cmake_minimum_required(VERSION 3.1)
    project(rtklib)

    add_definitions(-DTRACE)#enable trace in rtklib

    set(RTKLIB_SRC convkml.c convrnx.c datum.c download.c ephemeris.c geoid.c ionex.c lambda.c options.c pntpos.c postpos.c ppp_ar.c ppp.c preceph.c qzslex.c rcvraw.c rinex.c rtcm2.c rtcm3.c rtcm3e.c rtcm.c rtkcmn.c rtkpos.c rtksvr.c sbas.c solution.c stream.c streamsvr.c tle.c)
    add_library(rtklib STATIC ${RTKLIB_SRC})

    install(TARGETS rtklib
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

    install(FILES rtklib.h DESTINATION include)
    ```
    
    Compile RTKLIB into a static library

    ```
    mkdir build
    cd build
    cmake ..
    make
    sudo make install
    ```
    
    It should be found as `RTKLIB/src/build/librtklib.a`.


## User Manual

See `RTKLIB/doc/manual_2.4.2.pdf` or download from <https://www.rtklib.com/rtklib_document.htm>


## Samples

* Download from <https://www.rtklib.com/rtklib_tutorial.htm> (The corresponding tutorial is for Windows)

* Use test data in `RTKLIB/test/data/rinex`


## Usage of Example Codes

### Compile with CMakeLists.txt

In `learning_rtklib/build`, compile the example codes
```
mkdir build
cd build
cmake ..
make
```

To use them, run
```
./spp_example.cpp
```

### Debug Mode

To debug in VSCode, add the following line into `"args"` in the `.vscode/tasks.json`
```
"-l:librtklib.a",
```

The file path in the codes should be modified:

`"../data/rinex/daej229a00.20n"` -> `"data/rinex/daej229a00.20n"`.


## Usage of GUI apps in Windows

### Raw data conversion: RTKCONV

* If the input file is formatted in RTCM, specify the version of RTCM. 

  Windows: change the file extension.
  
  Linux: use `-r` to specify the format type.

* In __Options...__, select RINEX version and satellite systems of the output file.

### Post-processing: RTKPOST

* In __Options...__:

  * In __Setting 1__, select positioning mode (RTK: Kinematic), satellite systems.
  * In __Positions__, fill in the position of base station according to obs file. Select appropriate coordinates.
