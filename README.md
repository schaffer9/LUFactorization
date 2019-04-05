This repository implements a LU solver for dense linear systems.

# Usage
Make sure that you have julia installed and that the following packages are installed:
- CSV
- DataFrames

You can install the packages as follows:
```
    $ julia
    julia> using Pkg
    julia> Pkg.add("CSV")
    julia> Pkg.add("DataFrames")
```

After doing this you can run the code with ``julia main.jl``. This will generate three *csv* files which provide
runtime and accuracy data for *LU factorization*, *forward substitution* and *backward substitution*.

