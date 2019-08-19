# FullReco_SimuPhase1

The full reconstruction code for Geant4 simulations of the solenoid-type experiment at FRS.

The code consists in the framework developed by the HypHI collaboration for the two first
HypHI experiments (phase0 & phase0.5). It is responsible for the track and event
reconstruction of simulated and experimental data. It uses the external library, Genfit,
for the Kalman Filter. It is written in C++14.

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
field/

The following directories are mandatory : geo/, src/, src/.deps and lib/. 
If lib or .deps is missing just:
```sh
mkdir lib
mkdir src/.deps

```
| Dir:   | Usage:                                            |
|--------|---------------------------------------------------|
| geo/   | Gathers geometry rootfiles                        | 
| src/   | Source directory                                  |
| lib/   | Installation directory for the compiled libraries |
| input/ | Gathers parameter files for calibrations          |
| field/ | Gathers field map files                           |

## Requirements

External: ROOT v6 + boost + Genfit library
Build: make + gcc > 5.1 or clang > 3.8 

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
--log / -l : setting the level of logging: -1 = quiet mode, no stdout output / 0 = warning and higher / 1 = info and higher / 2 = debug and higher (default being 1) 


Usage / help message:

!> Example of use: \
!> ./MainSimu [-c nb_cpu] [--cpu nb_cpu] [-n fraction] [--num fraction] [-s start_ev] [--start start_ev] [-e nb_event] [--event nb_event] [-l lvllog] [--log lvllog] [-h] [--help] OutputFile RootInputFile

It takes first the command line options, then the output rootfile name, then one input rootfile name.


## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Merge Request
