# Numerical Simulation Laboratory (NSL)
### _Alberto Lazzeri_

These are NSL exercises during A.A.2023-2024 edition. The `c++` code is used for the simulations, whereas the data analysis is in a `jupyter-notebook`.

## Structure
Each folder matches a lecture and contains a `makefile`, a `c++` main for each sub-exercise (e.g. [es_2_1.cpp](/lezione_2/es_2_1.cpp)), various `c++` libraries, simulation results in `csv` files and a single `jupyter-notebook` with data analysis.

## Run code
To run simulation code a GCC compiler is needed. STL is sufficient to build and run the simulation code. For data analysis check the dependences above each "Consegna" in the `jupyter-notebook`.
Some extra libraries are required depending on the exercise:
- For Lecture 4 [`armadillo`](https://arma.sourceforge.net/) libraries are required to build the simulations.
- For Lecture 4 [`ovito`](https://www.ovito.org/) `python` package is needed for data visualization.

To create the executable using `makefile` use for each specific exercise:

```shell
make es_x_x.exe
```
Then to run the simulation use:

```shell
./es_x_x.exe
```

To remove object files and executable use:

```shell
make clean
```
Other extra will be indicated below.