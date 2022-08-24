# ART - Advanced Radiative Transfer code 

LTE RT code for solar ALMA applications
(Formerly known as RTma.)


## Installation

To run & compile, you will need C++ & Fortran compiler with MPI. Additional that you will need a parallel version of the HDF5 library. On MacOS (M1, or Intel) the easiest is to install all requirements via `Home Brew`. You will need the following packages:

```shell
brew install gcc openmpi hdf5-mpi
``` 

On Linux machine @ITA you most likely have all necessary compilers and libraries installed, and you just need to load the right modules (see compilation section). 

## Compilation

To compile, enter the `src` subfolder and type `make MACHINE=xxx`, where `xxx` is the file's name from the `machine` subfolder. These files are for machine-specific setup (compilers, flags etc.). You might want to create a new one for your own machine (host) or use an existing one. For example, if you have installed `gcc` & `openmpi` on Mac Os through Home Brew you can simply use:

```shell
make MACHINE=mac_gnu
```

Successful installation should generate `RTma.x` executable:

```shell
...
mpicxx -O2 -mtune=native -std=c++11 -I/opt/homebrew/include/ -g  -c main.cc -o main.o
mpicxx  -o RTma.x cop.o  eoswrap.o witt.o model.o io.o clte.o main.o -lstdc++ -lc++ -L/opt/homebrew/lib -lhdf5 -lhdf5_hl
```

## Input files

To run ART you need an input atmosphere in `*.h5` (HDF5) format. In the example folder you can find a small model with following structure:

```shell

h5ls model.h5

dens                     Dataset {1, 25, 15, 269}
dx                       Dataset {15}
dy                       Dataset {25}
temperature              Dataset {1, 25, 15, 269}
xne                      Dataset {1, 25, 15, 269}
z                        Dataset {1, 25, 15, 269}
```

Where:

	* `dens` - is density in CGS units
	* `dx`,`dy` - mesh spacing in the X-Y direction
	* temperature - ..well obvious but make sure it is in [K]
	* xne - electron density in CGS [ 1/cm^3 ]
	* z - z-axis spacing in [cm]

There is also an example of a python script for generating input models from Bifrost snapshots in raw format (or official FITS) in `utils`.

### Input file with frequencies 

Additional to the input model atmosphere, you will also need an input file to define a list of frequencies for which you wish to calculate intensities. Example of such file (you can use it for template) can be also found in `example` folder. 

```
# well obvious ..input and out files 
input_model = model.h5
output_profiles = synthetic.h5

# do not touch 
lines_file = lines.vald.input
rt_solver = 0
verbose = 1
mu = 1.0

# EOS can be witt (default) or pisk (slower)
eos = witt

# Cut the model from above until T = temperature_cut
temperature_cut = 50000.0
		
# Compute contribution function (default = 0)
get_contribution_funtion = 0


# Options (stderr, stdout, master, file, null)
logfile = master

# Finally, a list of frequencies.. wavelengths in Angstrom units // first is ~1 mm, and second ~3.1 mm
region = 10000000., 1.0, 1, 1.0, none, none
region = 31000000., 1.0, 1, 1.0, none, none
region = 6300.000, 0.01, 16, 1.0, none, none
region = 5000.000, 1.0, 1, 1.0, none, none
```	

### Running ART

Ok, we have the model and input file.. to run use `mpirun`

```shell
cd example
mpirun -np 4 ./RTma.x -i input.cfg
```

expected output:

```shell

➜ mpirun -np 4 ./RTma.x -i input.cfg

Using input file: input.cfg
info: readInput: input_model -> model.h5
info: readInput: output_profiles -> synthetic.h5
info: readInput: lines_file -> lines.vald.input
info: readInput: rt_solver -> 0
info: readInput: verbose -> 1
info: readInput: mu -> 1.0
info: readInput: eos -> witt
info: readInput: temperature_cut -> 50000.0
info: readInput: get_contribution_funtion -> 0
log target 0x20952b978
info: readInput: region -> 10000000.,1.0,1,1.0,none,none
info: readInput: region -> 31000000.,1.0,1,1.0,none,none
info: readInput: region -> 6300.000,0.01,16,1.0,none,none
info: readInput: region -> 5000.000,1.0,1,1.0,none,none
log target 0x20952cda8
log target 0x20952cda8
log target 0x20952cda8
test log 0/4
info: io, read_lines: read [3] lines
info: initReadIO: Found dims [ nt=1, ny=25, nx=15, ndep=269 ]
info: initReadIO: Found vars [ z ][ temperature ][ dens ][ xne ]
info: processing 100%
```

## Output file

After a successful run you should get an HDF5 under the name you set in the input file ( `synthetic.h5` in the example ) with the following structure:

```

➜ h5ls synthetic.h5

Stokes_I                 Dataset {1, 25, 15, 19}
Tau1                     Dataset {1, 25, 15, 19}
Wavelength               Dataset {19}
```

Where:

	*  `Stokes_I` - Calculated intensities, the first dimension is time (not used), then `x`,`y` and `wavelength`.
	*  `Wavelength` - is the axis for 4th dimension (list of wavelengths in Angstroms) 
	*  `Tau1` - is formation hight in `cm` i.e., value on `z-axis`  where `Tau_{wavelenght} ~= 1.0`.

These should be easy to plot via Python // see example folder again. 