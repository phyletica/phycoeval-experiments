#! /usr/bin/env julia

using Statistics
using DataDeps
using DataFrames
using CSV
using Plots
using StatsPlots
# using GR

include("project_util.jl")


struct SimResults
    df::DataFrame
    function SimResults()
        results = DataFrame()
        for is_fixed in [true, false]
            for sim_condition in ["bifurcating", "generalized"]
                sim_dir_name = ""
                if is_fixed  
                    sim_dir_name = "species-9-genomes-2-$sim_condition-tree-provided-fixed"
                else
                    sim_dir_name = "species-9-genomes-2-$sim_condition-tree-random-unfixed"
                end
                sim_dir = joinpath(ProjectUtil.SIM_DIR, sim_dir_name)
                for batch_dir in ProjectUtil.batch_dir_iter(sim_dir)
                    for analysis_condition in ["bifurcating", "generalized"]
                        for is_var_only in [true, false]
                            result_file_name = "results-species-9-genomes-2-$analysis_condition-tree-random.tsv.gz"
                            if is_var_only
                                result_file_name = "results-var-only-species-9-genomes-2-$analysis_condition-tree-random.tsv.gz"
                            end
                            result_path = joinpath(batch_dir, result_file_name)
                            df = ProjectUtil.get_data_frame_from_gz([result_path])
                            nrows = size(df)[1]
                            df[!, "sim_fixed"] = repeat([is_fixed], nrows)
                            df[!, "sim_generalized"] = repeat([sim_condition == "generalized"], nrows)
                            df[!, "analysis_generalized"] = repeat([analysis_condition == "generalized"], nrows)
                            df[!, "only_var_sites"] = repeat([is_var_only], nrows)
                            results = vcat(results, df)
                        end
                    end
                end
            end
        end
        new(results)
    end
end

function get_rows(r::SimResults,
        fixed::Bool,
        sim_generalized::Bool,
        analysis_generalized::Bool,
        var_only::Bool
        )::DataFrame
    return r.df[(r.df[!, :sim_fixed] .== fixed) .& (r.df[!, :sim_generalized] .== sim_generalized) .& (r.df[!, :analysis_generalized] .== analysis_generalized) .& (r.df[!, :only_var_sites] .== var_only), :]
end

function get_floats(r::SimResults,
        fixed::Bool,
        sim_generalized::Bool,
        analysis_generalized::Bool,
        var_only::Bool,
        col_header::Symbol
        )::Vector{Float64}
    return r.df[(r.df[!, :sim_fixed] .== fixed) .& (r.df[!, :sim_generalized] .== sim_generalized) .& (r.df[!, :analysis_generalized] .== analysis_generalized) .& (r.df[!, :only_var_sites] .== var_only), :][!, col_header]
end

function get_ints(r::SimResults,
        fixed::Bool,
        sim_generalized::Bool,
        analysis_generalized::Bool,
        var_only::Bool,
        col_header::Symbol
        )::Vector{Int64}
    return r.df[(r.df[!, :sim_fixed] .== fixed) .& (r.df[!, :sim_generalized] .== sim_generalized) .& (r.df[!, :analysis_generalized] .== analysis_generalized) .& (r.df[!, :only_var_sites] .== var_only), :][!, col_header]
end

function parse_results()::DataFrame
end

function get_bools(r::SimResults,
        fixed::Bool,
        sim_generalized::Bool,
        analysis_generalized::Bool,
        var_only::Bool,
        col_header::Symbol
        )::Vector{Bool}
    return r.df[(r.df[!, :sim_fixed] .== fixed) .& (r.df[!, :sim_generalized] .== sim_generalized) .& (r.df[!, :analysis_generalized] .== analysis_generalized) .& (r.df[!, :only_var_sites] .== var_only), :][!, col_header]
end

function main_cli()::Cint
    if ! Base.Filesystem.ispath(ProjectUtil.RESULTS_DIR)
        Base.Filesystem.mkdir(ProjectUtil.RESULTS_DIR)
    end

    # Set backend to GR
    # gr()
    write(stdout, "Plotting backend: $(backend())\n")

    brooks_gelman_1998_recommended_psrf = 1.2

    results = SimResults()

    root_node_probs = get_floats(results, true, true, true, false, :root_node_true_prob)
    vo_root_node_probs = get_floats(results, true, true, true, true, :root_node_true_prob)
    vln_root_node_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(root_node_probs, vo_root_node_probs),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Root node",
            #= size = (600, 400), =#
    )

    node_789_probs = get_floats(results, true, true, true, false, :node_7_8_9_prob)
    vo_node_789_probs = get_floats(results, true, true, true, true, :node_7_8_9_prob)
    vln_node_789_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(node_789_probs, vo_node_789_probs),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Node 789",
            #= size = (600, 400), =#
    )

    node_456_probs = get_floats(results, true, true, true, false, :node_4_5_6_prob)
    vo_node_456_probs = get_floats(results, true, true, true, true, :node_4_5_6_prob)
    vln_node_456_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(node_456_probs, vo_node_456_probs),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Node 456",
            #= size = (600, 400), =#
    )

    node_123_probs = get_floats(results, true, true, true, false, :node_12_3_prob)
    vo_node_123_probs = get_floats(results, true, true, true, true, :node_12_3_prob)
    write(stdout, "$node_123_probs\n")
    write(stdout, "$vo_node_123_probs\n")
    vln_node_123_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(node_123_probs, vo_node_123_probs),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Node 123",
            #= size = (600, 400), =#
    )

    height_12_789_probs = get_floats(results, true, true, true, false, :height_12_789_prob)
    vo_height_12_789_probs = get_floats(results, true, true, true, true, :height_12_789_prob)
    vln_height_12_789_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(height_12_789_probs, vo_height_12_789_probs),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Height 12-789",
            #= size = (600, 400), =#
    )

    height_123_456_probs = get_floats(results, true, true, true, false, :height_123_456_prob)
    vo_height_123_456_probs = get_floats(results, true, true, true, true, :height_123_456_prob)
    vln_height_123_456_probs = StatsPlots.violin(["All sites" "SNPs"],
            hcat(height_123_456_probs, vo_height_123_456_probs),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Height 123-456",
            #= size = (600, 400), =#
    )

    split_12_probs = get_floats(results, true, true, true, false, :split_12_prob)
    vo_split_12_probs = get_floats(results, true, true, true, true, :split_12_prob)
    vln_split_12_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(split_12_probs, vo_split_12_probs),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Split 12",
            #= size = (600, 400), =#
    )

    true_topo_probs = get_floats(results, true, true, true, false, :topo_true_prob)
    vo_true_topo_probs = get_floats(results, true, true, true, true, :topo_true_prob)
    vln_true_topo_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(true_topo_probs, vo_true_topo_probs),
            ylims = (0.0, 1.0),
            legend = false,
            title = "True topology",
            #= size = (600, 400), =#
    )

    vln_grid = Plots.plot(
            vln_root_node_probs,
            vln_node_789_probs,
            vln_node_456_probs,
            vln_node_123_probs,
            vln_height_12_789_probs,
            vln_height_123_456_probs,
            vln_split_12_probs,
            vln_true_topo_probs,
            layout = (2, 4),
            legend = false,
            size = (800, 400),
    )

    Plots.savefig(vln_grid, joinpath(ProjectUtil.RESULTS_DIR, "node-height-probs.pdf"))

    return 0
end

main_cli()
