#!/usr/bin/env julia


using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "path"
        help = "The csv file path where the result is stored"
        required = true
    "--threads", "-t"
        help = "Number of BLAS threads"
        arg_type = Int
        default = 5
    "--size", "-s"
        help = "The maximum system size to solve (min 100)"
        arg_type = String
        default = "3000"
    "--blocked", "-b"
        help = "Run the blocked version of the algorithm"
        action = :store_true
    "--block-size", "-g"
        help = "The block size if the factorization is blocked"
        arg_type = Int
        default = 100
    "--no-pivoting", "-p"
        help = "Do not apply pivoting"
        action = :store_false
    "--runs", "-r"
        help = "The number of runs to calculate the average time and accuracy"
        arg_type = Int
        default = 5
end
PARSED_ARGS = parse_args(ARGS, s)


using LUFactorization
using LinearAlgebra
using CSV


function main(parsed_args)

    size = eval(Meta.parse(parsed_args["size"]))
    directory_path = dirname(parsed_args["path"])
    mkpath(directory_path)
    file_path = joinpath(parsed_args["path"])

    BLAS.set_num_threads(parsed_args["threads"])

    result = run_benchmark_lu_fac(
        size;
        blocked=parsed_args["blocked"],
        block_size=parsed_args["block-size"],
        pivoting=parsed_args["no-pivoting"],
        runs=parsed_args["runs"]
    )
    CSV.write(file_path, result, writeheaders=true)

end

if  !isinteractive()
    main(PARSED_ARGS)
end
