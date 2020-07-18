#! /usr/bin/env julia

using Statistics
using ArgParse
using YAML
using DataDeps
using DataFrames
using CSV

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
    sim_name = Base.Filesystem.basename(Base.Filesystem.dirname(batch_dir))
    write(Base.stdout, "$sim_name\n")
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
    results = DataFrame()
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
            analysis_label = get_sim_label(sim_pattern)
            tree_log_wildcard = (
                    "run-?-$(sim_pattern.var_only)$(sim_pattern.prefix)" *
                    "-sim-$sim_num-$(sim_pattern.config_name)" *
                    "-config-trees-run-1.nex"
                   )
            treesum = run_sumphycoeval(
                    [tree_log_wildcard],
                    true_tree_path,
                    burnin,
                    min_split_freq)
            if number_of_runs == 0
                number_of_runs = length(treesum["summary_of_tree_sources"]["sources"])
            else
                @assert number_of_runs == length(treesum["summary_of_tree_sources"]["sources"])
            end
            if number_of_samples == 0
                number_of_samples = treesum["summary_of_tree_sources"]["total_number_of_trees_sampled"]
            else
                @assert number_of_samples == treesum["summary_of_tree_sources"]["total_number_of_trees_sampled"]
            end
            state_log_paths = [Base.Filesystem.joinpath(batch_dir,
                    ("run-$i-$(sim_pattern.var_only)$(sim_pattern.prefix)" *
                    "-sim-$sim_num-$(sim_pattern.config_name)" *
                    "-config-state-run-1.log")
                   ) for i in 1:number_of_runs]
            state_dfs = [ProjectUtil.get_data_frame(
                    [p],
                    skip = burnin
                    ) for p in state_log_paths]
            # Some sanity checks to ensure that the summary provided by
            # sumphycoeval matches what is in the state log files
            n_state_samples = sum(length(d[!, "ln_likelihood"]) for d in state_dfs)
            state_nheights = vcat((d[!, "number_of_heights"] for d in state_dfs)...)
            state_root_sizes = vcat((d[!, "pop_size_root"] for d in state_dfs)...)
            median_nheights = Statistics.median(state_nheights)
            mean_nheights = Statistics.mean(state_nheights)
            mean_root_size = Statistics.mean(state_root_sizes)
            @assert n_state_samples == number_of_samples
            @assert median_nheights == treesum["number_of_heights_summary"]["median"]
            @assert isapprox(
                    mean_nheights,
                    treesum["number_of_heights_summary"]["mean"])
            @assert isapprox(
                    mean_root_size,
                    treesum["clades"]["root"]["pop_size_mean"])
            stdout_paths = [Base.Filesystem.joinpath(batch_dir,
                    ("run-$i-$(sim_pattern.var_only)$(sim_pattern.prefix)" *
                    "-sim-$sim_num-$(sim_pattern.config_name)" *
                    "-config.yml.out")
                   ) for i in 1:number_of_runs]
            runtimes_str = [ProjectUtil.parse_runtime(p) for p in stdout_paths]
            runtimes = parse.(Int, runtimes_str)
            @assert length(runtimes) == number_of_runs
            lnl_ess = ProjectUtil.effective_sample_size(state_dfs, "ln_likelihood")
            lnl_psrf = ProjectUtil.potential_scale_reduction_factor(state_dfs, "ln_likelihood")
            lnp_ess = ProjectUtil.effective_sample_size(state_dfs, "ln_prior")
            lnp_psrf = ProjectUtil.potential_scale_reduction_factor(state_dfs, "ln_prior")
            leaf_label_map = Dict{String, Int}()
            for item in treesum["leaf_label_map"]
                leaf_label_map[item.second] = item.first
            end
            row = DataFrame(
             number_of_runs = [number_of_runs],
             number_of_samples = [number_of_samples],
             run_time_mean = [Statistics.mean(runtimes)],
             run_time_median = [Statistics.median(runtimes)],
             run_time_min = [minimum(runtimes)],
             run_time_max = [maximum(runtimes)],
             ln_likelihood_ess = [lnl_ess],
             ln_likelihood_psrf = [lnl_psrf],
             ln_prior_ess = [lnp_ess],
             ln_prior_psrf = [lnp_psrf],
             sdsf_mean = treesum["summary_of_split_freq_std_deviations"]["mean"],
             sdsf_max = treesum["summary_of_split_freq_std_deviations"]["max"],
             topo_true_prob = treesum["summary_of_target_tree"]["frequency"],
             topo_true_cred_level = treesum["summary_of_target_tree"]["credibility_level"],
             topo_true_is_map = treesum["summary_of_target_tree"]["is_a_map_topology"],
             num_heights_true = treesum["summary_of_target_tree"]["number_of_heights"],
             num_heights_true_prob = treesum["summary_of_target_tree"]["number_of_heights_frequency"],
             num_heights_true_cred_level = treesum["summary_of_target_tree"]["number_of_heights_credibility_level"],
             num_heights_map = treesum["numbers_of_heights"][1]["number_of_heights"],
             num_heights_mean = treesum["number_of_heights_summary"]["mean"],
             num_heights_median = treesum["number_of_heights_summary"]["median"],
             num_heights_eti_95_lower = treesum["number_of_heights_summary"]["eti_95"][1],
             num_heights_eti_95_upper = treesum["number_of_heights_summary"]["eti_95"][2],
             num_heights_hpdi_95_lower = treesum["number_of_heights_summary"]["hpdi_95"][1],
             num_heights_hpdi_95_upper = treesum["number_of_heights_summary"]["hpdi_95"][2],
             num_heights_ess = treesum["number_of_heights_summary"]["ess"],
             tree_length_true = treesum["summary_of_target_tree"]["tree_length"],
             tree_length_true_percentile = treesum["summary_of_target_tree"]["tree_length_percentile"],
             tree_length_mean = treesum["tree_length"]["mean"],
             tree_length_median = treesum["tree_length"]["median"],
             tree_length_eti_95_lower = treesum["tree_length"]["eti_95"][1],
             tree_length_eti_95_upper = treesum["tree_length"]["eti_95"][2],
             tree_length_hpdi_95_lower = treesum["tree_length"]["hpdi_95"][1],
             tree_length_hpdi_95_upper = treesum["tree_length"]["hpdi_95"][2],
             tree_length_ess = treesum["tree_length"]["ess"],
             tree_length_psrf = treesum["tree_length"]["psrf"],
             root_height_true = treesum["summary_of_target_tree"]["clades"]["root"]["height"],
             root_height_true_percentile = treesum["summary_of_target_tree"]["clades"]["root"]["height_percentile"],
             root_height_mean = treesum["clades"]["root"]["height_mean"],
             root_height_median = treesum["clades"]["root"]["height_median"],
             root_height_eti_95_lower = treesum["clades"]["root"]["height_eti_95"][1],
             root_height_eti_95_upper = treesum["clades"]["root"]["height_eti_95"][2],
             root_height_hpdi_95_lower = treesum["clades"]["root"]["height_hpdi_95"][1],
             root_height_hpdi_95_upper = treesum["clades"]["root"]["height_hpdi_95"][2],
             root_height_ess = treesum["clades"]["root"]["height_ess"],
             root_height_psrf = treesum["clades"]["root"]["height_psrf"],
             root_pop_size_true = treesum["summary_of_target_tree"]["clades"]["root"]["pop_size"],
             root_pop_size_true_percentile = treesum["summary_of_target_tree"]["clades"]["root"]["pop_size_percentile"],
             root_pop_size_mean = treesum["clades"]["root"]["pop_size_mean"],
             root_pop_size_median = treesum["clades"]["root"]["pop_size_median"],
             root_pop_size_eti_95_lower = treesum["clades"]["root"]["pop_size_eti_95"][1],
             root_pop_size_eti_95_upper = treesum["clades"]["root"]["pop_size_eti_95"][2],
             root_pop_size_hpdi_95_lower = treesum["clades"]["root"]["pop_size_hpdi_95"][1],
             root_pop_size_hpdi_95_upper = treesum["clades"]["root"]["pop_size_hpdi_95"][2],
             root_pop_size_ess = treesum["clades"]["root"]["pop_size_ess"],
             root_pop_size_psrf = treesum["clades"]["root"]["pop_size_psrf"],
             euclidean_distance_mean = treesum["summary_of_target_tree"]["euclidean_distance"]["mean"],
             euclidean_distance_median = treesum["summary_of_target_tree"]["euclidean_distance"]["median"],
             euclidean_distance_eti_95_lower = treesum["summary_of_target_tree"]["euclidean_distance"]["eti_95"][1],
             euclidean_distance_eti_95_upper = treesum["summary_of_target_tree"]["euclidean_distance"]["eti_95"][2],
             euclidean_distance_hpdi_95_lower = treesum["summary_of_target_tree"]["euclidean_distance"]["hpdi_95"][1],
             euclidean_distance_hpdi_95_upper = treesum["summary_of_target_tree"]["euclidean_distance"]["hpdi_95"][2],
            )
            if endswith(sim_name, "generalized-tree-provided-fixed")
                clade_12 = Set{Int}([
                        leaf_label_map["sp1"],
                        leaf_label_map["sp2"]
                       ])
                clade_123 = Set{Int}([
                        leaf_label_map["sp1"],
                        leaf_label_map["sp2"],
                        leaf_label_map["sp3"]
                       ])
                clade_456 = Set{Int}([
                        leaf_label_map["sp4"],
                        leaf_label_map["sp5"],
                        leaf_label_map["sp6"]
                       ])
                clade_45 = Set{Int}([
                        leaf_label_map["sp4"],
                        leaf_label_map["sp5"]
                       ])
                clade_46 = Set{Int}([
                        leaf_label_map["sp4"],
                        leaf_label_map["sp6"]
                       ])
                clade_56 = Set{Int}([
                        leaf_label_map["sp5"],
                        leaf_label_map["sp6"]
                       ])
                clade_789 = Set{Int}([
                        leaf_label_map["sp7"],
                        leaf_label_map["sp8"],
                        leaf_label_map["sp9"]
                       ])
                clade_78 = Set{Int}([
                        leaf_label_map["sp7"],
                        leaf_label_map["sp8"]
                       ])
                clade_79 = Set{Int}([
                        leaf_label_map["sp7"],
                        leaf_label_map["sp9"]
                       ])
                clade_89 = Set{Int}([
                        leaf_label_map["sp8"],
                        leaf_label_map["sp9"]
                       ])
                clade_123456 = Set{Int}([
                        leaf_label_map["sp1"],
                        leaf_label_map["sp2"],
                        leaf_label_map["sp3"],
                        leaf_label_map["sp4"],
                        leaf_label_map["sp5"],
                        leaf_label_map["sp6"]
                       ])
                clade_123789 = Set{Int}([
                        leaf_label_map["sp1"],
                        leaf_label_map["sp2"],
                        leaf_label_map["sp3"],
                        leaf_label_map["sp7"],
                        leaf_label_map["sp8"],
                        leaf_label_map["sp9"]
                       ])
                clade_456789 = Set{Int}([
                        leaf_label_map["sp4"],
                        leaf_label_map["sp5"],
                        leaf_label_map["sp6"],
                        leaf_label_map["sp7"],
                        leaf_label_map["sp8"],
                        leaf_label_map["sp9"]
                       ])
                height_12_789 = Set{Set{Int}}([clade_12, clade_789])
                height_123_456 = Set{Set{Int}}([clade_123, clade_456])

                clade_789_height_true = NaN
                clade_789_height_mean = NaN
                clade_789_height_median = NaN
                clade_789_height_eti_95_lower = NaN
                clade_789_height_eti_95_upper = NaN
                clade_789_height_hpdi_95_lower = NaN
                clade_789_height_hpdi_95_upper = NaN
                clade_789_height_ess = NaN
                clade_456_height_true = NaN
                clade_456_height_mean = NaN
                clade_456_height_median = NaN
                clade_456_height_eti_95_lower = NaN
                clade_456_height_eti_95_upper = NaN
                clade_456_height_hpdi_95_lower = NaN
                clade_456_height_hpdi_95_upper = NaN
                clade_456_height_ess = NaN

                clade_789_prob = 0.0
                clade_78_prob = 0.0
                clade_79_prob = 0.0
                clade_89_prob = 0.0
                clade_456_prob = 0.0
                clade_45_prob = 0.0
                clade_46_prob = 0.0
                clade_56_prob = 0.0
                clade_123456_prob = 0.0
                clade_123789_prob = 0.0
                clade_456789_prob = 0.0

                height_12_789_true = NaN
                height_12_789_prob = 0.0
                height_12_789_mean = NaN
                height_12_789_median = NaN
                height_12_789_eti_95_lower = NaN
                height_12_789_eti_95_upper = NaN
                height_12_789_hpdi_95_lower = NaN
                height_12_789_hpdi_95_upper = NaN
                height_12_789_ess = NaN
                height_123_456_true = NaN
                height_123_456_prob = 0.0
                height_123_456_mean = NaN
                height_123_456_median = NaN
                height_123_456_eti_95_lower = NaN
                height_123_456_eti_95_upper = NaN
                height_123_456_hpdi_95_lower = NaN
                height_123_456_hpdi_95_upper = NaN
                height_123_456_ess = NaN

                for c in treesum["clades"]["nontrivial_clades"]
                    leaf_set = Set{Int}(c["leaf_indices"])
                    if leaf_set == clade_789
                        clade_789_prob = c["frequency"]
                        clade_789_height_mean = c["height_mean"]
                        clade_789_height_median = c["height_median"]
                        clade_789_height_eti_95_lower = c["height_eti_95"][1]
                        clade_789_height_eti_95_upper = c["height_eti_95"][2]
                        clade_789_height_hpdi_95_lower = c["height_hpdi_95"][1]
                        clade_789_height_hpdi_95_upper = c["height_hpdi_95"][2]
                        clade_789_height_ess = c["height_ess"]
                    elseif leaf_set == clade_456
                        clade_456_prob = c["frequency"]
                        clade_456_height_mean = c["height_mean"]
                        clade_456_height_median = c["height_median"]
                        clade_456_height_eti_95_lower = c["height_eti_95"][1]
                        clade_456_height_eti_95_upper = c["height_eti_95"][2]
                        clade_456_height_hpdi_95_lower = c["height_hpdi_95"][1]
                        clade_456_height_hpdi_95_upper = c["height_hpdi_95"][2]
                        clade_456_height_ess = c["height_ess"]
                    elseif leaf_set == clade_78
                        clade_78_prob = c["frequency"]
                    elseif leaf_set == clade_79
                        clade_79_prob = c["frequency"]
                    elseif leaf_set == clade_89
                        clade_89_prob = c["frequency"]
                    elseif leaf_set == clade_45
                        clade_45_prob = c["frequency"]
                    elseif leaf_set == clade_46
                        clade_46_prob = c["frequency"]
                    elseif leaf_set == clade_56
                        clade_56_prob = c["frequency"]
                    elseif leaf_set == clade_123456
                        clade_123456_prob = c["frequency"]
                    elseif leaf_set == clade_123789
                        clade_123789_prob = c["frequency"]
                    elseif leaf_set == clade_456789
                        clade_456789_prob = c["frequency"]
                    end
                end
                for c in treesum["summary_of_target_tree"]["clades"]["nontrivial_clades"]
                    leaf_set = Set{Int}(c["leaf_indices"])
                    if leaf_set == clade_789
                        clade_789_height_true = c["height"]
                    elseif leaf_set == clade_456
                        clade_456_height_true= c["height"]
                    end
                end

                for h in treesum["heights"]
                    height_set = Set{Set{Int}}()
                    for c in h["clades"]
                        clade = Set{Int}(c["leaf_indices"])
                        push!(height_set, clade)
                    end
                    if height_set == height_12_789
                        height_12_789_prob = h["frequency"]
                        height_12_789_mean = h["mean"]
                        height_12_789_median = h["median"]
                        height_12_789_eti_95_lower = h["eti_95"][1]
                        height_12_789_eti_95_upper = h["eti_95"][2]
                        height_12_789_hpdi_95_lower = h["hpdi_95"][1]
                        height_12_789_hpdi_95_upper = h["hpdi_95"][2]
                        height_12_789_ess = h["ess"]
                    elseif height_set == height_123_456
                        height_123_456_prob = h["frequency"]
                        height_123_456_mean = h["mean"]
                        height_123_456_median = h["median"]
                        height_123_456_eti_95_lower = h["eti_95"][1]
                        height_123_456_eti_95_upper = h["eti_95"][2]
                        height_123_456_hpdi_95_lower = h["hpdi_95"][1]
                        height_123_456_hpdi_95_upper = h["hpdi_95"][2]
                        height_123_456_ess = h["ess"]
                    end
                end
                for h in treesum["summary_of_target_tree"]["heights"]
                    height_set = Set{Set{Int}}()
                    for c in h["clades"]
                        clade = Set{Int}(c["leaf_indices"])
                        push!(height_set, clade)
                    end
                    if height_set == height_12_789
                        height_12_789_true = h["height"]
                    elseif height_set == height_123_456
                        height_124_456_true = h["height"]
                    end
                end
                extra_cols = DataFrame(
                 clade_789_prob = clade_789_prob,
                 clade_78_prob = clade_78_prob,
                 clade_79_prob = clade_79_prob,
                 clade_89_prob = clade_89_prob,
                 clade_456_prob = clade_456_prob,
                 clade_45_prob = clade_45_prob,
                 clade_46_prob = clade_46_prob,
                 clade_56_prob = clade_56_prob,
                 clade_123456_prob = clade_123456_prob,
                 clade_123789_prob = clade_123789_prob,
                 clade_456789_prob = clade_456789_prob,
                 clade_789_height_true = clade_789_height_true,
                 clade_789_height_mean = clade_789_height_mean,
                 clade_789_height_median = clade_789_height_median,
                 clade_789_height_eti_95_lower = clade_789_height_eti_95_lower,
                 clade_789_height_eti_95_upper = clade_789_height_eti_95_upper,
                 clade_789_height_hpdi_95_lower = clade_789_height_hpdi_95_lower,
                 clade_789_height_hpdi_95_upper = clade_789_height_hpdi_95_upper,
                 clade_456_height_true = clade_456_height_true,
                 clade_456_height_mean = clade_456_height_mean,
                 clade_456_height_median = clade_456_height_median,
                 clade_456_height_eti_95_lower = clade_456_height_eti_95_lower,
                 clade_456_height_eti_95_upper = clade_456_height_eti_95_upper,
                 clade_456_height_hpdi_95_lower = clade_456_height_hpdi_95_lower,
                 clade_456_height_hpdi_95_upper = clade_456_height_hpdi_95_upper,
                 height_12_789_true = height_12_789_true,
                 height_12_789_mean = height_12_789_mean,
                 height_12_789_median = height_12_789_median,
                 height_12_789_eti_95_lower = height_12_789_eti_95_lower,
                 height_12_789_eti_95_upper = height_12_789_eti_95_upper,
                 height_12_789_hpdi_95_lower = height_12_789_hpdi_95_lower,
                 height_12_789_hpdi_95_upper = height_12_789_hpdi_95_upper,
                 height_12_789_ess = height_12_789_ess,
                 height_123_456_true = height_123_456_true,
                 height_123_456_mean = height_123_456_mean,
                 height_123_456_median = height_123_456_median,
                 height_123_456_eti_95_lower = height_123_456_eti_95_lower,
                 height_123_456_eti_95_upper = height_123_456_eti_95_upper,
                 height_123_456_hpdi_95_lower = height_123_456_hpdi_95_lower,
                 height_123_456_hpdi_95_upper = height_123_456_hpdi_95_upper,
                 height_123_456_ess = height_123_456_ess,
                )
                row = hcat(row, extra_cols)
            end
            results = vcat(results, row)
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
