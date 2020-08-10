#! /usr/bin/env julia

using ArgParse
using Random

include("project_util.jl")


function main_cli()::Cint
    parser = ArgParseSettings()

    @add_arg_table! parser begin
        "--seed"
                arg_type = Int
                range_tester = x -> x > 0
                action = :store_arg
                help = "Seed for random number generator."
        "number_of_ids"
                arg_type = Int
                range_tester = x -> x > 0
                action = :store_arg
                default = 1
                help = "Number of batch IDs."
    end

    parsed_args = parse_args(parser)

    rng = Random.MersenneTwister()
    if isnothing(parsed_args["seed"])
        parsed_args["seed"] = Random.rand(1:typemax(Int))
    end
    Random.seed!(rng, parsed_args["seed"])

    for i in 1:parsed_args["number_of_ids"]
        batch_num_str::AbstractString = ProjectUtil.get_batch_id(rng, 9)

        write(stdout, "$batch_num_str\n")
    end

    return 0
end
    
main_cli()
