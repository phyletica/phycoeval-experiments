#! /usr/bin/env julia

using ArgParse
using YAML

include("project_util.jl")


function get_yaml_config(path::AbstractString)::Dict{AbstractString, Any}
    open(path, "r") do istream
        return YAML.load(istream)
    end
end

function convert_for_variable_sites_only!(config::Dict{AbstractString, Any})
    if haskey(config, "data") && haskey(config["data"], "constant_sites_removed")
        config["data"]["constant_sites_removed"] = true
    end
    return nothing
end


function write_variable_sites_only_config(config_path::AbstractString)
    new_path = Base.Filesystem.joinpath(
            Base.Filesystem.dirname(config_path),
            "var-only-$(Base.Filesystem.basename(config_path))")
    config::Dict{AbstractString, Any} = get_yaml_config(config_path)
    convert_for_variable_sites_only!(config)
    YAML.write_file(new_path, config)
    return nothing
end

function main_cli()::Cint
    parser = ArgParseSettings()

    @add_arg_table! parser begin
        "sim_dir"
                arg_type = AbstractString
                action = :store_arg
                nargs = '+'
                required = true
                range_tester = x -> Base.Filesystem.isdir(x)
                help = "Path to directory with simphycoeval output files."
    end

    parsed_args = parse_args(parser)

    for sim_dir in parsed_args["sim_dir"]
        config_paths = ProjectUtil.simphycoeval_config_iter(sim_dir)
        write_variable_sites_only_config.(config_paths)
    end

    return 0
end
    
main_cli()
