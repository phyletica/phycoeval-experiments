#! /usr/bin/env julia

using Statistics
using ArgParse
using YAML
using DataDeps

include("project_util.jl")


sim_archive_names = [
        "sim-files-true-trees.tar.gz"
        "sim-files-true-parameters.tar.gz"
        "sim-files-tree-logs.tar.gz"
        "sim-files-state-logs.tar.gz"
        "sim-files-stdout.tar.gz"
]

struct SimOutputPattern
    var_only::String
    prefix::String
    config_name::String
end

function get_sim_label(s::SimOutputPattern)::String
    return s.var_only * s.config_name
end

function check_archive_paths(batch_dir::AbstractString)
    for archive_name in sim_archive_names
        archive_path = Base.Filesystem.joinpath(batch_dir, archive_name)
        if ! Base.Filesystem.ispath(archive_path)
            throw(ArgumentError("Batch dir $batch_dir does not have archive $archive_name"))
        end
    end
    return nothing
end

function extract_sim_files(batch_dir::AbstractString)
    check_archive_paths(batch_dir)
    for archive_name in sim_archive_names
        archive_path = Base.Filesystem.joinpath(batch_dir, archive_name)
        DataDeps.unpack(archive_path, keep_originals = true)
    end
    return nothing
end

function run_sumphycoeval(
        tree_log_paths::Vector{String},
        target_tree_path::String,
        burnin::Int,
        min_split_freq::AbstractFloat
        )::Dict{AbstractString, Any}
    @assert burnin > -1
    @assert min_split_freq > 0.0
    @assert min_split_freq < 1.0
    sumphy_path = Base.Filesystem.joinpath(ProjectUtil.BIN_DIR, "sumphycoeval")
    sumphy_cmd = `sh -c "$sumphy_path --newick-target -b $burnin -t $target_tree_path --min-split-freq $min_split_freq $tree_log_paths"`
    return YAML.load(read(sumphy_cmd, String))
end

function get_sim_output_patterns(
        batch_dir::AbstractString
        )::Vector{SimOutputPattern}
    tree_log_pattern_str = (
            raw"run-1-" *
            raw"(?<var_only>var-only-)?" *
            raw"(?<config_prefix>\S*simphycoeval)" *
            raw"-sim-(?<sim_num>0+)" *
            raw"-(?<config_name>\S+)" *
            raw"-config-trees-run-(?P<dummy_run_num>\d+).nex"
           )
    tree_log_pattern = Regex(join([
            raw"^",
            tree_log_pattern_str,
            raw"$"
           ]))
    sim_output_patterns::Vector{SimOutputPattern} = []
    for file_name in Base.Filesystem.readdir(batch_dir)
        m = match(tree_log_pattern, file_name)
        if ! isnothing(m)
            var_only = ""
            prefix = ""
            if ! isnothing(m[:var_only])
                var_only = m[:var_only]
            end
            if ! isnothing(m[:config_prefix])
                prefix = m[:config_prefix]
            end
            config_name = m[:config_name]
            @assert m[:dummy_run_num] == "1"
            push!(sim_output_patterns, SimOutputPattern(var_only, prefix, config_name))
        end
    end
    return sim_output_patterns
end

function parse_sim_results(
        batch_dir::AbstractString,
        burnin::Int,
        min_split_freq::AbstractFloat,
        extract_sim_archives::Bool
       )
    cd(batch_dir)
    batch_dir = "."
    if extract_sim_archives
        write(Base.stdout, "extracting archives...\n")
        extract_sim_files(batch_dir)
    end

    # Get what output files to expect
    sim_output_patterns = get_sim_output_patterns(
            batch_dir)
    # Use true tree files to loop over sim replicates
    true_tree_file_iter = ProjectUtil.flat_file_path_iter(batch_dir,
            ProjectUtil.SIM_TRUE_TREE_FILE_PATTERN)
    number_of_runs = 0
    number_of_samples = 0
    for true_tree_path in true_tree_file_iter
        m = match(ProjectUtil.SIM_TRUE_TREE_FILE_PATTERN, true_tree_path)
        if isnothing(m)
            error("Bad true tree path '$true_tree_path'")
        end
        sim_num = m[:sim_num]
        true_params_path = Base.Filesystem.joinpath(batch_dir,
                "simphycoeval-sim-$sim_num-true-parameters.txt")
        true_params = ProjectUtil.get_data_frame([true_params_path], skip = burnin)
        for sim_pattern in sim_output_patterns
            tree_log_wildcard = (
                    "run-?-$(sim_pattern.var_only)$(sim_pattern.prefix)" *
                    "-sim-$sim_num-$(sim_pattern.config_name)" *
                    "-config-trees-run-1.nex"
                   )
            tree_summary = run_sumphycoeval(
                    [tree_log_wildcard],
                    true_tree_path,
                    burnin,
                    min_split_freq)
            if number_of_runs == 0
                number_of_runs = length(tree_summary["summary_of_tree_sources"]["sources"])
            else
                @assert number_of_runs == length(tree_summary["summary_of_tree_sources"]["sources"])
            end
            if number_of_samples == 0
                number_of_samples = tree_summary["summary_of_tree_sources"]["total_number_of_trees_sampled"]
            else
                @assert number_of_samples == tree_summary["summary_of_tree_sources"]["total_number_of_trees_sampled"]
            end
            state_log_paths = [Base.Filesystem.joinpath(batch_dir,
                    ("run-$i-$(sim_pattern.var_only)$(sim_pattern.prefix)" *
                    "-sim-$sim_num-$(sim_pattern.config_name)" *
                    "-config-state-run-1.log")
                   ) for i in 1:number_of_runs]
            state_df = ProjectUtil.get_data_frame(state_log_paths, skip = burnin)
            write(Base.stdout, "$(get_sim_label(sim_pattern))\n")
            write(Base.stdout, "  sample size: $number_of_samples\n")
            write(Base.stdout, "  n sources: $number_of_runs\n")
            write(Base.stdout, "  $(size(state_df))\n")
            # Some sanity checks to ensure that the summary provided by
            # sumphycoeval matches what is in the state log files
            @assert length(state_df[!, "ln_likelihood"]) == number_of_samples
            @assert Statistics.median(state_df[!, "number_of_heights"]) == tree_summary["number_of_heights_summary"]["median"]
            @assert isapprox(
                    Statistics.mean(state_df[!, "number_of_heights"]),
                    tree_summary["number_of_heights_summary"]["mean"])
            @assert isapprox(
                    Statistics.mean(state_df[!, "pop_size_root"]),
                    tree_summary["clades"]["root"]["pop_size_mean"])
            stdout_paths = [Base.Filesystem.joinpath(batch_dir,
                    ("run-$i-$(sim_pattern.var_only)$(sim_pattern.prefix)" *
                    "-sim-$sim_num-$(sim_pattern.config_name)" *
                    "-config.yml.out")
                   ) for i in 1:number_of_runs]
            runtimes = [ProjectUtil.parse_runtime(p) for p in stdout_paths]
            @assert length(runtimes) == number_of_runs
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
        "--burnin", "-b"
                arg_type = Int
                range_tester = x -> x > 0
                action = :store_arg
                default = 501
                help = ("Number of MCMC samples from each analysis to ignore "
                        * "as burnin.")
        "--min-split-freq", "-m"
                arg_type = AbstractFloat
                range_tester = x -> ((x > 0.0) & (x < 1.0))
                action = :store_arg
                default = 0.1
                help = ("Minimum frequency of splits for calculating the " *
                        "average standard deviation of split frequencies " *
                        "among MCMC samples.")
        "--skip-archives", "-s"
                action = :store_true
                help = "Skip extraction of sim-file archives."
    end

    parsed_args = parse_args(parser)

    for sim_dir in parsed_args["sim_dir"]
        parse_sim_results(sim_dir,
                          parsed_args["burnin"],
                          parsed_args["min-split-freq"],
                          (! parsed_args["skip-archives"])
                         )
    end
    return 0
end

main_cli()