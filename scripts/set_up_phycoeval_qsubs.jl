#! /usr/bin/env julia

using ArgParse

include("project_util.jl")


function write_qsub(config_path::AbstractString,
        run_number::Int = 1,
        relax_missing_sites::Bool = false)
    qsub_prefix = Base.Filesystem.splitext(config_path)[1]
    qsub_path = "$qsub_prefix-run-$run_number-qsub.sh"
    if Base.Filesystem.ispath(qsub_path)
        write(Base.stdout, "Skipping $qsub_path\n")
        return nothing
    end
    config_file = Base.basename(config_path)
    stdout_path = "run-$run_number-$config_file.out"
    exe_var_name = "exe_path"
    open(qsub_path, "w") do out
        write(out, ProjectUtil.get_pbs_header(qsub_path,
                exe_name = "phycoeval",
                exe_var_name = exe_var_name))
        if relax_missing_sites
            write(out, "\"\$$exe_var_name\" --seed $run_number --prefix \"run-$run_number-\" --relax-constant-sites --relax-missing-sites \"$config_file\" 1>\"$stdout_path\" 2>&1\n")
        else
            write(out, "\"\$$exe_var_name\" --seed $run_number --prefix \"run-$run_number-\" --relax-constant-sites \"$config_file\" 1>\"$stdout_path\" 2>&1\n")
        end
    end
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
        "--number-of-runs", "-n"
                arg_type = Int
                range_tester = x -> x > 0
                action = :store_arg
                default = 4
                help = "Number of qsubs to generate per config."
    end

    parsed_args = parse_args(parser)

    for sim_dir in parsed_args["sim_dir"]
        config_paths = ProjectUtil.simphycoeval_config_iter(sim_dir)
        for cfg_path in config_paths
            for i in 1:parsed_args["number-of-runs"]
                write_qsub(cfg_path, i, false)
            end
        end
    end

    return 0
end

main_cli()
