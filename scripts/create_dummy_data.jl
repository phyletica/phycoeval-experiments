#! /usr/bin/env julia

println(Base.active_project())
using ArgParse

function write_dummy_biallelic_data_file(;
        nspecies::Int = 2,
        ngenomes::Int = 10,
        ncharacters::Int = 1000,
        prefix::AbstractString = "sp",
        out::IO = Base.stdout)
    nspecies_padding = length(string(nspecies))
    write(out, "---\n")
    write(out, "markers_are_dominant: false\n")
    write(out, "population_labels: [")
    for sp_idx in 1:nspecies
        if sp_idx > 1
            write(out, ", ")
        end
        sp_label = prefix * lpad(sp_idx, nspecies_padding, '0')
        write(out, "$sp_label")
    end
    write(out, "]\n")
    write(out, "allele_count_patterns:\n")
    write(out, "    - [")
    for sp_idx in 1:nspecies
        if sp_idx > 1
            write(out, ", ")
        end
        write(out, "[1, $ngenomes]")
    end
    write(out, "]\n")
    write(out, "pattern_weights: [$ncharacters]\n")
    return nothing
end

function main_cli()::Cint
    parser = ArgParseSettings()

    @add_arg_table! parser begin
        "--nspecies", "-s"
                arg_type = Int
                range_tester = x -> x > 0
                action = :store_arg
                default = 2
                help = "The number of populations."
        "--ngenomes", "-g"
                arg_type = Int
                range_tester = x -> x > 0
                action = :store_arg
                default = 10
                help = "The number of genomes sampled per population."
        "--ncharacters", "-c"
                arg_type = Int
                range_tester = x -> x > 0
                action = :store_arg
                default = 1000
                help = "The number of biallelic characters."
        "--prefix", "-p"
                arg_type = AbstractString
                action = :store_arg
                default = "sp"
                help = "Prefix for species labels."
    end

    #= parsed_args = parse_args(ARGS, parser) =#
    parsed_args = parse_args(parser)

    write_dummy_biallelic_data_file(
            nspecies = parsed_args["nspecies"],
            ngenomes = parsed_args["ngenomes"],
            ncharacters = parsed_args["ncharacters"],
            prefix = parsed_args["prefix"],
            out = Base.stdout)
    return 0
end
    
main_cli()
