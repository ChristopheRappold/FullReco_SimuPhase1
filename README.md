# FullReco_SimuPhase1

The full reconstruction code for Geant4 simulations of the solenoid-type experiment at FRS.

The code consists in the framework developed by the HypHI collaboration for the two first
HypHI experiments (phase0 & phase0.5). It is responsible for the track and event
reconstruction of simulated and experimental data. It uses the external library, Genfit,
for the Kalman Filter. It is written in C++14.

The framework uses several extrenal packages or libraries, either internal programming or
physics analysis purposes.
The current list of the dependencies and their purposes is :

| Libraries:    | Purposes:                                                                         |
|---------------|-----------------------------------------------------------------------------------|
| ROOT v6       | General analysis framework in particle & nuclear physics                          | 
| boost         | Library for Boost.PropertyTree                                                    |
| Genfit v2     | Library for Kalman Filter track fitting                                           |
| Eigen3        | Library for linear algebra fast SIMD computation                                  |
| spdlog        | Library for fast and multithread logging                                          |
| msgpack-c     | Library for efficient binary serialization to use in message passing              |
| zeroMQ        | Library for high-performance asynchronous messaging in distributed computation    |
| cppzmq        | Library for C++ handling of ZeroMQ framework                                      |
| sqlite3       | Library for sqlite database to use in book keeping configuration, calibration etc |
| sqlite_orm    | Library for modern C++ handling of sqlite database                                |
| KFParticle    | Library for Kalman Filter vertex fitting from CBM, STAR &  ALICE                  |
| Vc            | Library for SIMD vector data structure and calculation                            |
| TrickTrack    | Library for track finding via Cellular Automaton from CMS                         |


## Installation

OS X & Linux:

1.  first:

```sh
git clone git@gitlab.com:HypHI-GSI/FullReco_SimuPhase1.git
```

2.  then you need to activate the submodule for Genfit:

```sh
cd FullReco_SimuPhase1
git submodule init
git submodule update
```

3.  Because the sub repository of Genfit will be only on a detached head, I think it is
better to assign a branch to the current and local repository:

```sh
cd src/Genfit2_v2
git checkout -b hyphidev remotes/origin/HypHIdev
```

Now all repositories are set.

4.  The local directory should look like this:

geo/ \
src/ \
src/.deps \
lib/ \
input/ \
field/ \
config/

The following directories are mandatory : geo/, src/, src/.deps, lib/ and config/. 
If lib or .deps or config is missing just:
```sh
mkdir lib
mkdir src/.deps
mkdir config 
```
| Dir:    | Usage:                                            |
|---------|---------------------------------------------------|
| geo/    | Gathers geometry rootfiles                        | 
| src/    | Source directory                                  |
| lib/    | Installation directory for the compiled libraries |
| input/  | Gathers parameter files for calibrations          |
| field/  | Gathers field map files                           |
| config/ | Gathers the configuration files                   |

## Requirements

External: ROOT v6 + boost + Genfit library + Eigen3 + spdlog + msgpack-c + zeroMQ + cppzmq + sqlite_orm + KFParticle + Vc + TrickTrack \
Build: make + gcc > 9.3 or clang > 9 \
Libraries from package manager of the system : libboost-dev + libeigen3-dev + libsqlite3-dev \

Your $PATH must include ROOT bin directory. Example:
```sh
echo $PATH
/home/christophe/root-6/bin:/usr/local/bin:/usr/bin:/bin:/usr/bin/X11
```
Your $LD_LIBRARY_PATH must include `pwd`/lib and ROOT lib directory. Example:
```sh
echo $LD_LIBRARY_PATH
./lib:/home/christophe/root-6/lib
```

In src/Genfit2_v2: you need to have a build directory:
```sh
cd src/Genfit2_v2
mkdir build
```
Then configure Genfit build with:
```sh
cd build
cmake ..
```

Genfit can be built:
```sh
make 
make install
```

In src/spdlog: you need also to have a build directory:
```sh
cd src/spdlog
mkdir build
```

Now spdlog can configure with:
```sh
cd build
cmake .. -DSPDLOG_BUILD_EXAMPLE=OFF -DSPDLOG_BUILD_TESTS=OFF -DCMAKE_INSTALL_PREFIX=..
```

And spdlog can be built:
```sh
make
make install
```

Now msgpack-c (version 3.3.0) must be cloned, configured and built with:
```sh
cd src/
git clone https://github.com/msgpack/msgpack-c.git
cd msgpack-c
git checkout tags/cpp-3.3.0 -b v3.3.0
cmake -DMSGPACK_CXX17=ON -DMSGPACK_BUILD_EXAMPLES=OFF .
make
```

Now libzmq (version 4.3.3) & cppzmq must be cloned, configured and built with:
``` sh
cd src/
git clone https://github.com/zeromq/libzmq.git
cd libzmq
git checkout tags/v4.3.3 -b v4.3.3
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=.. ..
make
make install

cd src/
git clone https://github.com/zeromq/cppzmq
cd cppzmq
mkdir build
cd build
cmake -DZeroMQ_DIR=/#ZEROMQ_INSTALL#/lib/cmake/ZeroMQ -DCMAKE_INSTALL_PREFIX=.. ..
make
make install
```

Now TrickTrack must be cloned, configured and built with:
``` sh
cd src/
git clone https://github.com/HSF/TrickTrack
cd TrickTrack
mkdir build
cd build
export EIGEN_INCLUDE_DIR=/usr/include/eigen3
cmake -DCMAKE_INSTALL_PREFIX=.. ..
make
make install
```

Now sqlite_orm must be only cloned with:
```sh
cd src/
git clone https://github.com/fnc12/sqlite_orm
```

Now Vc and KFParticle must be cloned, configured and built with :
```sh
cd src/
git clone https://github.com/VcDevel/Vc
cd Vc
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=../install -DBUILD_TESTING=OFF ..
make -j
make install

export Vc_DIR=/#VC_INSTALL_DIR#/lib/cmake/Vc

cd src/
git clone git@gitlab.com:HypHI-GSI/KFParticle.git
cd KFParticle
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=../install ../ -DFIXTARGET=TRUE
make -j
make install
```


Then to build everything (-j# is for parallel compilation with # being number of cpu):
```sh
cd src
make -j#
cd ..
```

You will have the libraries in lib/ and the executable ./MainSimu


## Usage example

./MainSimu -e 1000 -g ./geo/GeoWasaRealSolenoidFront.root Output.root run_H3LWithA_WasaFront3_Field1.8T.root

command line options:

--help / -h : prompt the help \
--cpu / -c : number of cpus \
--num / -n : cpu id for event splitting \
--event / -e : total number of event to analyze \
--start / -s : start from this event id \
--geo / -g : path to geometry rootfile to use \
--config / -f : path to the configuration file to use \
--log / -l : setting the level of logging: -1 = quiet mode, no stdout output / 0 = warning and higher / 1 = info and higher / 2 = debug and higher (default being 1) 


Usage / help message:

!> Example of use: \
!> ./MainSimu [-f Config.cfg] [--config Config.cfg] [-g Geofile] [--geo Geofile] [-c nb_cpu] [--cpu nb_cpu] [-n fraction] [--num fraction] [-s start_ev] [--start start_ev] [-e nb_event] [--event nb_event] [-l lvllog] [--log lvllog] [-h]  OutputFile RootInputFile 

It takes first the command line options, then the output rootfile name, then one input rootfile name.


## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Merge Request
