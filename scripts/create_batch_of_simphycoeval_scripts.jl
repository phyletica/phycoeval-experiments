#! /usr/bin/env julia

using ArgParse
using Random

include("project_util.jl")


function get_sim_name(sim_config_path::AbstractString,
        fix_sim_model::Bool = false)::AbstractString
    cfg_name::AbstractString = Base.Filesystem.basename(
            Base.Filesystem.splitext(sim_config_path)[1])
    fixed_label::AbstractString = fix_sim_model ? "fixed" : "unfixed"

    return "$cfg_name-$fixed_label"
end

function get_sim_script_path(batch_id::AbstractString,
        sim_config_path::AbstractString,
        fix_sim_model::Bool = false)::AbstractString
    script_name = "simphycoeval-$(get_sim_name(sim_config_path, fix_sim_model))-batch-$batch_id.sh"
    return Base.Filesystem.joinpath(ProjectUtil.SIM_SCRIPT_DIR, script_name)

end

function write_sim_script(
        batch_id::AbstractString,
        number_of_reps::Int,
        number_of_topo_mcmc_gens_per_rep::Int,
        sim_config_path::AbstractString,
        prior_config_paths::AbstractVector{AbstractString},
        fix_sim_model::Bool = false)
    script_path::AbstractString = get_sim_script_path(batch_id,
            sim_config_path,
            fix_sim_model)
    sim_name::AbstractString = get_sim_name(sim_config_path, fix_sim_model)
    sim_cfg_file::AbstractString = Base.Filesystem.basename(sim_config_path)
    open(script_path, "w") do ostream
        write(ostream, ProjectUtil.get_pbs_header(script_path,
                exe_name = "simphycoeval"))
        write(ostream, "rng_seed=$batch_id\n")
        write(ostream, "number_of_reps=$number_of_reps\n")
        write(ostream, "sim_name=\"$sim_name\"\n")
        write(ostream, "config_dir=\"../../configs\"\n")
        write(ostream, "config_path=\"\${config_dir}/$sim_cfg_file\"\n")
        write(ostream, "output_dir=\"../../simulations/\${sim_name}/batch-\${rng_seed}\"\n")
        write(ostream, "config_set_up_script_path=\"\${project_dir}/scripts/set_up_configs_for_simulated_data.jl\"\n")
        write(ostream, "qsub_set_up_script_path=\"\${project_dir}/scripts/set_up_phycoeval_qsubs.jl\"\n")
        write(ostream, "\n")
        write(ostream, "mkdir -p \"\$output_dir\"\n")
        write(ostream, "\n")
        write(ostream, "\$exe_path --seed=\"\$rng_seed\" -n \"\$number_of_reps\" -t \"\$number_of_topo_mcmc_gens_per_rep\" -o \"\$output_dir\"")
        if fix_sim_model
            write(ostream, " --fix-model")
        end
        for prior_cfg_path in prior_config_paths
            prior_cfg_name = Base.Filesystem.basename(prior_cfg_path)
            write(ostream, " -p \"\${config_dir}/$prior_cfg_name\"")
        end
        write(ostream, " \"\$config_path\"")
        write(ostream, " && julia --project \"\$config_set_up_script_path\" \"\$output_dir\"")
        write(ostream, " && julia --project \"\$qsub_set_up_script_path\" \"\$output_dir\"")
        write(ostream, "\n")
    end
    return nothing
end

function main_cli()::Cint
    parser = ArgParseSettings()

    @add_arg_table! parser begin
        "--nreps", "-n"
                arg_type = Int
                range_tester = x -> x > 0
                action = :store_arg
                default = 10
                help = "Number of simulation replicates."
        "--topo-mcmc-gens-per-rep", "-t"
                arg_type = Int
                range_tester = x -> x > 0
                action = :store_arg
                default = 100
                dest_name = "topo_mcmc_gens_per_rep"
                help = ("Number of topology altering MCMC generations per "
                        * "simulation replicate.")
        "--prior", "-p"
                arg_type = AbstractString
                action = :append_arg
                required = true
                range_tester = x -> Base.Filesystem.isfile(x)
                dest_name = "prior_config_path"
                help = "Path to prior configs."
        "--fixed", "-f"
                arg_type = AbstractString
                action = :append_arg
                required = false
                range_tester = x -> Base.Filesystem.isfile(x)
                dest_name = "fixed_sim_config_path"
                help = "Path to fixed simulation configs."
        "sim_config_path"
                arg_type = AbstractString
                action = :store_arg
                nargs = '+'
                required = true
                range_tester = x -> Base.Filesystem.isfile(x)
                help = "Paths to simulation configs."
        "--seed"
                arg_type = Int
                range_tester = x -> x > 0
                action = :store_arg
                help = "Seed for random number generator."
    end

    parsed_args = parse_args(parser)

    rng = Random.MersenneTwister()
    if isnothing(parsed_args["seed"])
        parsed_args["seed"] = Random.rand(1:typemax(Int))
    end
    Random.seed!(rng, parsed_args["seed"])

    num_seed_digits::Int = 9
    batch_num_str::AbstractString = Random.randstring(rng, '0':'9',
                                                      num_seed_digits)

    try
        Base.Filesystem.mkdir(ProjectUtil.SIM_SCRIPT_DIR)
    catch err
        if ! isa(err, Base.IOError)
            throw(err)
        end
    end

    for sim_cfg_path in parsed_args["fixed_sim_config_path"]
        write_sim_script(batch_num_str,
                         parsed_args["nreps"],
                         parsed_args["topo_mcmc_gens_per_rep"],
                         sim_cfg_path,
                         parsed_args["prior_config_path"],
                         true)
    end
    for sim_cfg_path in parsed_args["sim_config_path"]
        write_sim_script(batch_num_str,
                         parsed_args["nreps"],
                         parsed_args["topo_mcmc_gens_per_rep"],
                         sim_cfg_path,
                         parsed_args["prior_config_path"],
                         false)
    end

    return 0
end
    
main_cli()
