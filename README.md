# Numerical Simulation Laboratory (NSL)
### _Alberto Lazzeri_

These are NSL exercises during A.A.2023-2024 edition. The `c++` code is used for the simulations, whereas the data analysis is in a `jupyter-notebook`.

## Structure
Each folder matches a lecture and contains a `makefile`, a `c++` main for each sub-exercise (e.g. [es_2_1.cpp](/lezione_2/es_2_1.cpp)), various `c++` libraries, simulation results in `csv` files and a single `jupyter-notebook` with data analysis.

## Run code
To run simulation code a GCC compiler is needed. STL is sufficient to build and run the simulation code. For data analysis check the dependences above each "Consegna" in the `jupyter-notebook`.
Some extra libraries are required depending on the exercise (specificated below).

To create the executable using `makefile` use for each specific exercise:

```shell
make es_x_x.exe
```
Then to run the simulation use:

```shell
./es_x_x.exe <eventual excericise parameters>
```

To remove object files and executable use:

```shell
make clean
```

### Lecture 4,6,7

For these lectures [`armadillo`](https://arma.sourceforge.net/) libraries are required to build the simulations.

Before running the a NVE simulation, equilibration is needed. In the input file set the parameters of the phase required, then with `python3` run:
```shell
python3 equilibration.py
```
This will add a line in the input file which reports the temperature which start the simulation at to reach the real input temperature at equilibrium.
To build the executable run in the folder "/NSL_SIMULATOR/SOURCE/":
```shell
make simulator.exe
```
Then to run the code phase specification is needed:
```shell
./simulator.exe <phase>
```

Other extra will be indicated below.