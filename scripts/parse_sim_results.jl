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

function pull_sim_archives(batch_dir::AbstractString)
    relative_batch_dir = Base.Filesystem.relpath(batch_dir, ProjectUtil.PROJECT_DIR)
    for archive_name in sim_archive_names
        archive_path = Base.Filesystem.joinpath(relative_batch_dir, archive_name)
        git_cmd = `git lfs pull --include $archive_path`
        run(git_cmd)
    end
    return nothing
end

function extract_sim_files(batch_dir::AbstractString)
    check_archive_paths(batch_dir)
    pull_sim_archives(batch_dir)
    for archive_name in sim_archive_names
        archive_path = Base.Filesystem.joinpath(batch_dir, archive_name)
        DataDeps.unpack(archive_path, keep_originals = true)
    end
    return nothing
end

function remove_extracted_files(batch_dir::AbstractString)
    check_archive_paths(batch_dir)
    for archive_name in sim_archive_names
        archive_path = Base.Filesystem.joinpath(batch_dir, archive_name)
        run(pipeline(`tar tzf $archive_path`, `xargs rm --`))
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
    sumphy_cmd = `sh -c "$sumphy_path --newick-target -b $burnin -t $target_tree_path --min-split-freq $min_split_freq --include-merged-target-heights $tree_log_paths"`
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
    return_dir = pwd()
    cd(batch_dir)
    batch_dir = "."
    if extract_sim_archives
        write(Base.stdout, "Extracting archives...\n")
        extract_sim_files(batch_dir)
    end

    extra_results_paths = Vector{AbstractString}()

    # Get what output files to expect
    sim_output_patterns = get_sim_output_patterns(
            batch_dir)
    results_paths = Dict{SimOutputPattern, String}()
    for sim_pattern in sim_output_patterns
        res_path = Base.Filesystem.joinpath(
                batch_dir,
                "results-$(sim_pattern.var_only)$(sim_pattern.config_name).tsv")
        if Base.Filesystem.ispath(res_path)
            write(Base.stdout, "Table '$(basename(res_path))' already exists; skipping!")
            continue
        end
        if Base.Filesystem.ispath(res_path * ".gz")
            write(Base.stdout, "Table '$(basename(res_path)).gz' already exists... Skipping!\n")
            continue
        end
        results_paths[sim_pattern] = res_path
    end
    if length(results_paths) < 1
        if extract_sim_archives
            write(Base.stdout, "Cleaning up extracted files...\n")
            remove_extracted_files(batch_dir)
        end
        cd(return_dir)
        return nothing
    end

    results = Dict{SimOutputPattern, DataFrame}()
    for sim_pattern in keys(results_paths) 
        results[sim_pattern] = DataFrame()
    end
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
        for sim_pattern in keys(results_paths)
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
                    treesum["splits"]["root"]["pop_size_mean"])
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

            true_polytomy_probs_path = Base.Filesystem.joinpath(
                    batch_dir,
                    "results-true-polytomy-probs-$(sim_pattern.var_only)$(sim_pattern.config_name).tsv")
            true_shared_probs_path = Base.Filesystem.joinpath(
                    batch_dir,
                    "results-true-shared-height-probs-$(sim_pattern.var_only)$(sim_pattern.config_name).tsv")
            merged_height_probs_path = Base.Filesystem.joinpath(
                    batch_dir,
                    "results-wrong-merged-height-probs-$(sim_pattern.var_only)$(sim_pattern.config_name).tsv")

            true_split_probs::Vector{Float64} = []
            true_node_probs::Vector{Float64} = []
            for true_split in treesum["summary_of_target_tree"]["splits"]["nontrivial_splits"]
                true_split_freq::Float64 = true_split["frequency"]
                push!(true_split_probs, true_split_freq)
                true_node_freq = true_split_freq
                if haskey(true_split, "node")
                    true_node_freq = true_split["node"]["frequency"]
                    num_descendants::Int = true_split["node"]["number_of_descendants"]
                    if num_descendants > 2
                        if ! Base.Filesystem.ispath(true_polytomy_probs_path)
                            open(true_polytomy_probs_path, "w") do out
                                write(out, "polytomy_node_prob\tnum_descendants\tis_root\n")
                            end
                            push!(extra_results_paths, true_polytomy_probs_path)
                        end
                        open(true_polytomy_probs_path, "a") do out
                            write(out, "$(true_node_freq)\t$(num_descendants)\t0\n")
                        end
                    end
                end
                push!(true_node_probs, true_node_freq)
            end

            n_root_descendants::Int = treesum["summary_of_target_tree"]["splits"]["root"]["node"]["number_of_descendants"]
            if n_root_descendants > 2
                root_node_freq::Float64 = treesum["summary_of_target_tree"]["splits"]["root"]["node"]["frequency"]
                if ! Base.Filesystem.ispath(true_polytomy_probs_path)
                    open(true_polytomy_probs_path, "w") do out
                        write(out, "polytomy_node_prob\tnum_descendants\tis_root\n")
                    end
                    push!(extra_results_paths, true_polytomy_probs_path)
                end
                open(true_polytomy_probs_path, "a") do out
                    write(out, "$(root_node_freq)\t$(n_root_descendants)\t1\n")
                end
            end

            true_height_probs::Vector{Float64} = []
            true_shared_height_probs::Vector{Float64} = []
            for true_height in treesum["summary_of_target_tree"]["heights"]
                push!(true_height_probs, true_height["frequency"])
                if length(true_height["splits"]) > 1
                    push!(true_shared_height_probs, true_height["frequency"])
                    if ! Base.Filesystem.ispath(true_shared_probs_path)
                        open(true_shared_probs_path, "w") do out
                            write(out, "shared_height_prob\tnum_splits\n")
                        end
                        push!(extra_results_paths, true_shared_probs_path)
                    end
                    open(true_shared_probs_path, "a") do out
                        write(out, "$(true_height["frequency"])\t$(length(true_height["splits"]))\n")
                    end
                end
            end

            # If true tree was comb, it won't have merged target heights
            if haskey(treesum, "merged_target_heights")
                for merged_target_height in treesum["merged_target_heights"]
                    if ! Base.Filesystem.ispath(merged_height_probs_path)
                        open(merged_height_probs_path, "w") do out
                            write(out, "merged_height_prob\tnum_nodes\theight_diff\theight_midpoint\n")
                        end
                        push!(extra_results_paths, merged_height_probs_path)
                    end
                    ht_freq::Float64 = merged_target_height["frequency"]
                    n_nodes::Int = merged_target_height["merged_height_number_of_nodes"]
                    younger_ht::Float64 = merged_target_height["younger_height"]
                    older_ht::Float64 = merged_target_height["older_height"]
                    ht_diff::Float64 = older_ht - younger_ht
                    ht_midpoint::Float64 = (older_ht + younger_ht) / 2.0
                    open(merged_height_probs_path, "a") do out
                        write(out, "$(ht_freq)\t$(n_nodes)\t$(ht_diff)\t$(ht_midpoint)\n")
                    end
                end
            end

            true_split_prob_mean = NaN
            true_split_prob_median = NaN
            if length(true_split_probs) > 0
                true_split_prob_mean = Statistics.mean(true_split_probs)
                true_split_prob_median = Statistics.median(true_split_probs)
            end

            true_node_prob_mean = NaN
            true_node_prob_median = NaN
            if length(true_node_probs) > 0
                true_node_prob_mean = Statistics.mean(true_node_probs)
                true_node_prob_median = Statistics.median(true_node_probs)
            end

            true_height_prob_mean = NaN
            true_height_prob_median = NaN
            if length(true_height_probs) > 0
                true_height_prob_mean = Statistics.mean(true_height_probs)
                true_height_prob_median = Statistics.median(true_height_probs)
            end

            true_shared_height_prob_mean = NaN
            true_shared_height_prob_median = NaN
            if length(true_shared_height_probs) > 0
                true_shared_height_prob_mean = Statistics.mean(true_shared_height_probs)
                true_shared_height_prob_median = Statistics.median(true_shared_height_probs)
            end

            root_node_true_is_map = treesum["summary_of_target_tree"]["splits"]["root"]["node"]["is_a_map_node_given_split"]
            max_wrong_root_prob = treesum["splits"]["root"]["nodes"][1]["frequency"]
            if root_node_true_is_map
                if length(treesum["splits"]["root"]["nodes"]) > 1
                    max_wrong_root_prob = treesum["splits"]["root"]["nodes"][2]["frequency"]
                else
                    max_wrong_root_prob = 0
                end
            end

            row = DataFrame(
             number_of_runs = number_of_runs,
             number_of_samples = number_of_samples,
             run_time_mean = Statistics.mean(runtimes),
             run_time_median = Statistics.median(runtimes),
             run_time_min = minimum(runtimes),
             run_time_max = maximum(runtimes),
             ln_likelihood_ess = lnl_ess,
             ln_likelihood_psrf = lnl_psrf,
             ln_prior_ess = lnp_ess,
             ln_prior_psrf = lnp_psrf,
             sdsf_mean = treesum["summary_of_split_freq_std_deviations"]["mean"],
             sdsf_max = treesum["summary_of_split_freq_std_deviations"]["max"],
             topo_true_prob = treesum["summary_of_target_tree"]["frequency"],
             topo_true_cred_level = treesum["summary_of_target_tree"]["credibility_level"],
             topo_true_is_map = treesum["summary_of_target_tree"]["is_a_map_topology"],
             true_split_prob_mean = true_split_prob_mean,
             true_split_prob_median = true_split_prob_median,
             true_node_prob_mean = true_node_prob_mean,
             true_node_prob_median = true_node_prob_median,
             true_height_prob_mean = true_height_prob_mean,
             true_height_prob_median = true_height_prob_median,
             true_shared_height_prob_mean = true_shared_height_prob_mean,
             true_shared_height_prob_median = true_shared_height_prob_median,
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
             root_height_true = treesum["summary_of_target_tree"]["splits"]["root"]["height"],
             root_height_true_percentile = treesum["summary_of_target_tree"]["splits"]["root"]["height_percentile"],
             root_height_mean = treesum["splits"]["root"]["height_mean"],
             root_height_median = treesum["splits"]["root"]["height_median"],
             root_height_eti_95_lower = treesum["splits"]["root"]["height_eti_95"][1],
             root_height_eti_95_upper = treesum["splits"]["root"]["height_eti_95"][2],
             root_height_hpdi_95_lower = treesum["splits"]["root"]["height_hpdi_95"][1],
             root_height_hpdi_95_upper = treesum["splits"]["root"]["height_hpdi_95"][2],
             root_height_ess = treesum["splits"]["root"]["height_ess"],
             root_height_psrf = treesum["splits"]["root"]["height_psrf"],
             root_pop_size_true = treesum["summary_of_target_tree"]["splits"]["root"]["pop_size"],
             root_pop_size_true_percentile = treesum["summary_of_target_tree"]["splits"]["root"]["pop_size_percentile"],
             root_pop_size_mean = treesum["splits"]["root"]["pop_size_mean"],
             root_pop_size_median = treesum["splits"]["root"]["pop_size_median"],
             root_pop_size_eti_95_lower = treesum["splits"]["root"]["pop_size_eti_95"][1],
             root_pop_size_eti_95_upper = treesum["splits"]["root"]["pop_size_eti_95"][2],
             root_pop_size_hpdi_95_lower = treesum["splits"]["root"]["pop_size_hpdi_95"][1],
             root_pop_size_hpdi_95_upper = treesum["splits"]["root"]["pop_size_hpdi_95"][2],
             root_pop_size_ess = treesum["splits"]["root"]["pop_size_ess"],
             root_pop_size_psrf = treesum["splits"]["root"]["pop_size_psrf"],
             root_node_true_prob = treesum["summary_of_target_tree"]["splits"]["root"]["node"]["frequency"],
             root_node_true_is_map = root_node_true_is_map,
             max_wrong_root_prob = max_wrong_root_prob,
             euclidean_distance_mean = treesum["summary_of_target_tree"]["euclidean_distance"]["mean"],
             euclidean_distance_median = treesum["summary_of_target_tree"]["euclidean_distance"]["median"],
             euclidean_distance_eti_95_lower = treesum["summary_of_target_tree"]["euclidean_distance"]["eti_95"][1],
             euclidean_distance_eti_95_upper = treesum["summary_of_target_tree"]["euclidean_distance"]["eti_95"][2],
             euclidean_distance_hpdi_95_lower = treesum["summary_of_target_tree"]["euclidean_distance"]["hpdi_95"][1],
             euclidean_distance_hpdi_95_upper = treesum["summary_of_target_tree"]["euclidean_distance"]["hpdi_95"][2],
            )

            split_12 = Set{Int}([
                    leaf_label_map["sp1"],
                    leaf_label_map["sp2"]
                   ])
            split_13 = Set{Int}([
                    leaf_label_map["sp1"],
                    leaf_label_map["sp3"]
                   ])
            split_23 = Set{Int}([
                    leaf_label_map["sp2"],
                    leaf_label_map["sp3"]
                   ])
            split_123 = Set{Int}([
                    leaf_label_map["sp1"],
                    leaf_label_map["sp2"],
                    leaf_label_map["sp3"]
                   ])
            split_456 = Set{Int}([
                    leaf_label_map["sp4"],
                    leaf_label_map["sp5"],
                    leaf_label_map["sp6"]
                   ])
            split_45 = Set{Int}([
                    leaf_label_map["sp4"],
                    leaf_label_map["sp5"]
                   ])
            split_46 = Set{Int}([
                    leaf_label_map["sp4"],
                    leaf_label_map["sp6"]
                   ])
            split_56 = Set{Int}([
                    leaf_label_map["sp5"],
                    leaf_label_map["sp6"]
                   ])
            split_789 = Set{Int}([
                    leaf_label_map["sp7"],
                    leaf_label_map["sp8"],
                    leaf_label_map["sp9"]
                   ])
            split_78 = Set{Int}([
                    leaf_label_map["sp7"],
                    leaf_label_map["sp8"]
                   ])
            split_79 = Set{Int}([
                    leaf_label_map["sp7"],
                    leaf_label_map["sp9"]
                   ])
            split_89 = Set{Int}([
                    leaf_label_map["sp8"],
                    leaf_label_map["sp9"]
                   ])
            split_123456 = Set{Int}([
                    leaf_label_map["sp1"],
                    leaf_label_map["sp2"],
                    leaf_label_map["sp3"],
                    leaf_label_map["sp4"],
                    leaf_label_map["sp5"],
                    leaf_label_map["sp6"]
                   ])
            split_123789 = Set{Int}([
                    leaf_label_map["sp1"],
                    leaf_label_map["sp2"],
                    leaf_label_map["sp3"],
                    leaf_label_map["sp7"],
                    leaf_label_map["sp8"],
                    leaf_label_map["sp9"]
                   ])
            split_456789 = Set{Int}([
                    leaf_label_map["sp4"],
                    leaf_label_map["sp5"],
                    leaf_label_map["sp6"],
                    leaf_label_map["sp7"],
                    leaf_label_map["sp8"],
                    leaf_label_map["sp9"]
                   ])
            height_12_789 = Set{Set{Int}}([split_12, split_789])
            height_123_456 = Set{Set{Int}}([split_123, split_456])

            height_12_78 = Set{Set{Int}}([split_12, split_78])
            height_45_789 = Set{Set{Int}}([split_45, split_789])
            height_456_789 = Set{Set{Int}}([split_456, split_789])

            node_7_8_9 = Set{Set{Int}}([
                    Set{Int}([leaf_label_map["sp7"]]),
                    Set{Int}([leaf_label_map["sp8"]]),
                    Set{Int}([leaf_label_map["sp9"]])])
            node_4_5_6 = Set{Set{Int}}([
                    Set{Int}([leaf_label_map["sp4"]]),
                    Set{Int}([leaf_label_map["sp5"]]),
                    Set{Int}([leaf_label_map["sp6"]])])
            node_12_3 = Set{Set{Int}}([
                    Set{Int}([leaf_label_map["sp1"], leaf_label_map["sp2"]]),
                    Set{Int}([leaf_label_map["sp3"]])])
            node_1_2_3 = Set{Set{Int}}([
                    Set{Int}([leaf_label_map["sp1"]]),
                    Set{Int}([leaf_label_map["sp2"]]),
                    Set{Int}([leaf_label_map["sp3"]])])
            node_123_456_789 = Set{Set{Int}}([
                    Set{Int}([leaf_label_map["sp1"], leaf_label_map["sp2"], leaf_label_map["sp3"]]),
                    Set{Int}([leaf_label_map["sp4"], leaf_label_map["sp5"], leaf_label_map["sp6"]]),
                    Set{Int}([leaf_label_map["sp7"], leaf_label_map["sp8"], leaf_label_map["sp9"]])])

            split_789_height_true = NaN
            split_789_height_mean = NaN
            split_789_height_median = NaN
            split_789_height_eti_95_lower = NaN
            split_789_height_eti_95_upper = NaN
            split_789_height_hpdi_95_lower = NaN
            split_789_height_hpdi_95_upper = NaN
            split_789_height_ess = NaN
            split_456_height_true = NaN
            split_456_height_mean = NaN
            split_456_height_median = NaN
            split_456_height_eti_95_lower = NaN
            split_456_height_eti_95_upper = NaN
            split_456_height_hpdi_95_lower = NaN
            split_456_height_hpdi_95_upper = NaN
            split_456_height_ess = NaN

            split_789_prob = 0.0
            split_78_prob = 0.0
            split_79_prob = 0.0
            split_89_prob = 0.0
            split_456_prob = 0.0
            split_45_prob = 0.0
            split_46_prob = 0.0
            split_56_prob = 0.0
            split_123_prob = 0.0
            split_12_prob = 0.0
            split_13_prob = 0.0
            split_23_prob = 0.0
            split_123456_prob = 0.0
            split_123789_prob = 0.0
            split_456789_prob = 0.0

            node_7_8_9_prob = 0.0
            node_4_5_6_prob = 0.0
            node_12_3_prob = 0.0
            node_1_2_3_prob = 0.0
            node_123_456_789_prob = 0.0

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

            height_12_78_prob = 0.0
            height_45_789_prob = 0.0
            height_456_789_prob = 0.0

            for n in treesum["splits"]["root"]["nodes"]
                split_set = Set{Set{Int}}([Set{Int}(split) for split in n["descendant_splits"]])
                if split_set == node_123_456_789
                    node_123_456_789_prob = n["frequency"]
                end
            end

            for c in treesum["splits"]["nontrivial_splits"]
                leaf_set = Set{Int}(c["leaf_indices"])
                if leaf_set == split_789
                    split_789_prob = c["frequency"]
                    split_789_height_mean = c["height_mean"]
                    split_789_height_median = c["height_median"]
                    split_789_height_eti_95_lower = c["height_eti_95"][1]
                    split_789_height_eti_95_upper = c["height_eti_95"][2]
                    split_789_height_hpdi_95_lower = c["height_hpdi_95"][1]
                    split_789_height_hpdi_95_upper = c["height_hpdi_95"][2]
                    split_789_height_ess = c["height_ess"]
                    for n in c["nodes"]
                        split_set = Set{Set{Int}}([Set{Int}(split) for split in n["descendant_splits"]])
                        if split_set == node_7_8_9
                            node_7_8_9_prob = n["frequency"]
                        end
                    end
                elseif leaf_set == split_456
                    split_456_prob = c["frequency"]
                    split_456_height_mean = c["height_mean"]
                    split_456_height_median = c["height_median"]
                    split_456_height_eti_95_lower = c["height_eti_95"][1]
                    split_456_height_eti_95_upper = c["height_eti_95"][2]
                    split_456_height_hpdi_95_lower = c["height_hpdi_95"][1]
                    split_456_height_hpdi_95_upper = c["height_hpdi_95"][2]
                    split_456_height_ess = c["height_ess"]
                    for n in c["nodes"]
                        split_set = Set{Set{Int}}([Set{Int}(split) for split in n["descendant_splits"]])
                        if split_set == node_4_5_6
                            node_4_5_6_prob = n["frequency"]
                        end
                    end
                elseif leaf_set == split_123
                    split_123_prob = c["frequency"]
                    for n in c["nodes"]
                        split_set = Set{Set{Int}}([Set{Int}(split) for split in n["descendant_splits"]])
                        if split_set == node_12_3
                            node_12_3_prob = n["frequency"]
                        elseif split_set == node_1_2_3
                            node_1_2_3_prob = n["frequency"]
                        end
                    end
                elseif leaf_set == split_12
                    split_12_prob = c["frequency"]
                elseif leaf_set == split_13
                    split_13_prob = c["frequency"]
                elseif leaf_set == split_23
                    split_23_prob = c["frequency"]
                elseif leaf_set == split_78
                    split_78_prob = c["frequency"]
                elseif leaf_set == split_79
                    split_79_prob = c["frequency"]
                elseif leaf_set == split_89
                    split_89_prob = c["frequency"]
                elseif leaf_set == split_45
                    split_45_prob = c["frequency"]
                elseif leaf_set == split_46
                    split_46_prob = c["frequency"]
                elseif leaf_set == split_56
                    split_56_prob = c["frequency"]
                elseif leaf_set == split_123456
                    split_123456_prob = c["frequency"]
                elseif leaf_set == split_123789
                    split_123789_prob = c["frequency"]
                elseif leaf_set == split_456789
                    split_456789_prob = c["frequency"]
                end
            end
            for c in treesum["summary_of_target_tree"]["splits"]["nontrivial_splits"]
                leaf_set = Set{Int}(c["leaf_indices"])
                if leaf_set == split_789
                    split_789_height_true = c["height"]
                elseif leaf_set == split_456
                    split_456_height_true= c["height"]
                end
            end

            for h in treesum["heights"]
                height_set = Set{Set{Int}}()
                for c in h["splits"]
                    split = Set{Int}(c["leaf_indices"])
                    push!(height_set, split)
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
                elseif height_set == height_12_78
                    height_12_78_prob = h["frequency"]
                elseif height_set == height_45_789
                    height_45_789_prob = h["frequency"]
                elseif height_set == height_456_789
                    height_456_789_prob = h["frequency"]
                end
            end
            for h in treesum["summary_of_target_tree"]["heights"]
                height_set = Set{Set{Int}}()
                for c in h["splits"]
                    split = Set{Int}(c["leaf_indices"])
                    push!(height_set, split)
                end
                if height_set == height_12_789
                    height_12_789_true = h["height"]
                elseif height_set == height_123_456
                    height_123_456_true = h["height"]
                end
            end
            extra_cols = DataFrame(
             split_789_prob = split_789_prob,
             split_78_prob = split_78_prob,
             split_79_prob = split_79_prob,
             split_89_prob = split_89_prob,
             max_789_subsplit_prob = maximum((split_78_prob, split_79_prob, split_89_prob)),
             split_456_prob = split_456_prob,
             split_45_prob = split_45_prob,
             split_46_prob = split_46_prob,
             split_56_prob = split_56_prob,
             max_456_subsplit_prob = maximum((split_45_prob, split_46_prob, split_56_prob)),
             split_123_prob = split_123_prob,
             split_12_prob = split_12_prob,
             split_13_prob = split_13_prob,
             split_23_prob = split_23_prob,
             split_123456_prob = split_123456_prob,
             split_123789_prob = split_123789_prob,
             split_456789_prob = split_456789_prob,
             split_789_height_true = split_789_height_true,
             split_789_height_mean = split_789_height_mean,
             split_789_height_median = split_789_height_median,
             split_789_height_eti_95_lower = split_789_height_eti_95_lower,
             split_789_height_eti_95_upper = split_789_height_eti_95_upper,
             split_789_height_hpdi_95_lower = split_789_height_hpdi_95_lower,
             split_789_height_hpdi_95_upper = split_789_height_hpdi_95_upper,
             split_456_height_true = split_456_height_true,
             split_456_height_mean = split_456_height_mean,
             split_456_height_median = split_456_height_median,
             split_456_height_eti_95_lower = split_456_height_eti_95_lower,
             split_456_height_eti_95_upper = split_456_height_eti_95_upper,
             split_456_height_hpdi_95_lower = split_456_height_hpdi_95_lower,
             split_456_height_hpdi_95_upper = split_456_height_hpdi_95_upper,
             node_7_8_9_prob = node_7_8_9_prob,
             node_4_5_6_prob = node_4_5_6_prob,
             node_12_3_prob = node_12_3_prob,
             node_1_2_3_prob = node_1_2_3_prob,
             node_123_456_789_prob = node_123_456_789_prob,
             height_12_789_true = height_12_789_true,
             height_12_789_prob = height_12_789_prob,
             height_12_789_mean = height_12_789_mean,
             height_12_789_median = height_12_789_median,
             height_12_789_eti_95_lower = height_12_789_eti_95_lower,
             height_12_789_eti_95_upper = height_12_789_eti_95_upper,
             height_12_789_hpdi_95_lower = height_12_789_hpdi_95_lower,
             height_12_789_hpdi_95_upper = height_12_789_hpdi_95_upper,
             height_12_789_ess = height_12_789_ess,
             height_123_456_true = height_123_456_true,
             height_123_456_prob = height_123_456_prob,
             height_123_456_mean = height_123_456_mean,
             height_123_456_median = height_123_456_median,
             height_123_456_eti_95_lower = height_123_456_eti_95_lower,
             height_123_456_eti_95_upper = height_123_456_eti_95_upper,
             height_123_456_hpdi_95_lower = height_123_456_hpdi_95_lower,
             height_123_456_hpdi_95_upper = height_123_456_hpdi_95_upper,
             height_123_456_ess = height_123_456_ess,
             height_12_78_prob = height_12_78_prob,
             height_45_789_prob = height_45_789_prob,
             height_456_789_prob = height_456_789_prob,
            )
            row = hcat(row, extra_cols)
            results[sim_pattern] = vcat(results[sim_pattern], row)
        end
    end
    for (sim_pattern, res_path) in results_paths
        CSV.write(res_path,
                  results[sim_pattern],
                  delim = '\t',
                  append = false)
        run(`gzip $res_path`)
    end
    for path in extra_results_paths
        run(`gzip $path`)
    end
    if extract_sim_archives
        write(Base.stdout, "Cleaning up extracted files...\n")
        remove_extracted_files(batch_dir)
    end
    cd(return_dir)
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
