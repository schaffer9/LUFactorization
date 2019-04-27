This repository implements a LU solver for dense linear systems.

# Usage
To install the project, first install *Julia*. After you did so you can add the project as follows.
Open a new shell and type
```
    $ julia
    julia> using Pkg
    julia> Pkg.add("https://github.com/sschaffer92/LUFactorization")

```

After doing so you can download the `main.jl` file of the package and
call `./main.jl --help`


```
usage: main.jl [-t THREADS] [-s SIZE] [-b] [-g BLOCK-SIZE] [-p]
               [-r RUNS] [-h] path

positional arguments:
  path                  The csv file path where the result is stored

optional arguments:
  -t, --threads THREADS
                        Number of BLAS threads (type: Int64, default:
                        5)
  -s, --size SIZE       The maximum system size to solve (min 100)
                        (default: "3000")
  -b, --blocked         Run the blocked version of the algorithm
  -g, --block-size BLOCK-SIZE
                        The block size if the factorization is blocked
                        (type: Int64, default: 100)
  -p, --no-pivoting     Do not apply pivoting
  -r, --runs RUNS       The number of runs to calculate the average
                        time and accuracy (type: Int64, default: 5)
  -h, --help            show this help message and exit
```
