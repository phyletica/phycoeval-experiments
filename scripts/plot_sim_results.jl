#! /usr/bin/env julia

using Random
using Distributions
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
        main_tsv_path = joinpath(ProjectUtil.RESULTS_DIR, "results.tsv")
        if Base.Filesystem.isfile(main_tsv_path * ".gz")
            results = ProjectUtil.get_data_frame_from_gz([main_tsv_path * ".gz"])
            return new(results)
        end
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
        CSV.write(main_tsv_path,
                  results,
                  delim = '\t',
                  append = false)
        run(`gzip $main_tsv_path`)
        return new(results)
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
        col_header::Union{Symbol, AbstractString}
        )::Vector{Float64}
    return r.df[(r.df[!, :sim_fixed] .== fixed) .& (r.df[!, :sim_generalized] .== sim_generalized) .& (r.df[!, :analysis_generalized] .== analysis_generalized) .& (r.df[!, :only_var_sites] .== var_only), :][!, col_header]
end

function get_ints(r::SimResults,
        fixed::Bool,
        sim_generalized::Bool,
        analysis_generalized::Bool,
        var_only::Bool,
        col_header::Union{Symbol, AbstractString}
        )::Vector{Int64}
    return r.df[(r.df[!, :sim_fixed] .== fixed) .& (r.df[!, :sim_generalized] .== sim_generalized) .& (r.df[!, :analysis_generalized] .== analysis_generalized) .& (r.df[!, :only_var_sites] .== var_only), :][!, col_header]
end

function get_bools(r::SimResults,
        fixed::Bool,
        sim_generalized::Bool,
        analysis_generalized::Bool,
        var_only::Bool,
        col_header::Union{Symbol, AbstractString}
        )::Vector{Bool}
    return r.df[(r.df[!, :sim_fixed] .== fixed) .& (r.df[!, :sim_generalized] .== sim_generalized) .& (r.df[!, :analysis_generalized] .== analysis_generalized) .& (r.df[!, :only_var_sites] .== var_only), :][!, col_header]
end

function get_true_v_est_plot(r::SimResults,
        fixed::Bool,
        sim_generalized::Bool,
        analysis_generalized::Bool,
        var_only::Bool,
        col_header_prefix::String,
        use_median::Bool = false,
        axis_limit_buffer::AbstractFloat = 0.05
        )::Plots.Plot
    x = get_floats(r,
            fixed,
            sim_generalized,
            analysis_generalized,
            var_only,
            col_header_prefix * "_true")
    if use_median
        y = get_floats(r,
                fixed,
                sim_generalized,
                analysis_generalized,
                var_only,
                col_header_prefix * "_median")
    else
        y = get_floats(r,
                fixed,
                sim_generalized,
                analysis_generalized,
                var_only,
                col_header_prefix * "_mean")
    end
    y_lower = get_floats(r,
            fixed,
            sim_generalized,
            analysis_generalized,
            var_only,
            col_header_prefix * "_eti_95_lower")
    y_upper = get_floats(r,
            fixed,
            sim_generalized,
            analysis_generalized,
            var_only,
            col_header_prefix * "_eti_95_upper")

    extremes = (minimum((minimum(x), minimum(y))),
                maximum((maximum(x), maximum(y))))
    xy_range = extremes[2] - extremes[1]
    xy_buffer = xy_range * axis_limit_buffer
    axis_limits = [extremes[1] - xy_buffer,
                   extremes[2] + xy_buffer]
    identity_line = [extremes[1] - (xy_buffer * 2.0)
                     extremes[2] + (xy_buffer * 2.0)]

    plt = Plots.plot(identity_line, identity_line,
            seriestype = :line,
            legend = false,
            xlims = axis_limits,
            ylims = axis_limits)
    Plots.plot!(plt,
            x,
            y,
            yerror = (y - y_lower, y_upper - y),
            seriestype = :scatter,
            link = :all,
            legend = false,
            xlims = axis_limits,
            ylims = axis_limits)
end

function get_shared_x_limits(
        plots::Plots.Plot...
        )::Tuple{Real, Real}
    x_min = Inf
    x_max = -Inf
    for plt in plots
        for subplt in plt.subplots
            x_min = minimum((x_min, subplt[:xaxis][:lims]...))
            x_max = maximum((x_max, subplt[:xaxis][:lims]...))
        end
    end
    return (x_min, x_max)
end

function get_shared_y_limits(
        plots::Plots.Plot...
        )::Tuple{Real, Real}
    y_min = Inf
    y_max = -Inf
    for plt in plots
        for subplt in plt.subplots
            y_min = minimum((y_min, subplt[:yaxis][:lims]...))
            y_max = maximum((y_max, subplt[:yaxis][:lims]...))
        end
    end
    return (y_min, y_max)
end

function get_shared_xy_limits(
        plots::Plots.Plot...
        )::Tuple{Real, Real}
    mn = Inf
    mx = -Inf
    for plt in plots
        for subplt in plt.subplots
            mn = minimum((mn, subplt[:xaxis][:lims]..., subplt[:yaxis][:lims]...))
            mx = maximum((mx, subplt[:xaxis][:lims]..., subplt[:yaxis][:lims]...))
        end
    end
    return (mn, mx)
end

function share_x_limits!(
        plots::Plots.Plot...
        )::Tuple{Real, Real}
    x_lims = get_shared_x_limits(plots...)
    for plt in plots
        xlims!(plt, x_lims)
    end
    return x_lims
end

function share_y_limits!(
        plots::Plots.Plot...
        )::Tuple{Real, Real}
    y_lims = get_shared_y_limits(plots...)
    for plt in plots
        ylims!(plt, y_lims)
    end
    return y_lims
end

function share_xy_limits!(
        plots::Plots.Plot...
        )::Tuple{Real, Real}
    xy_lims = get_shared_xy_limits(plots...)
    for plt in plots
        xlims!(plt, xy_lims)
        ylims!(plt, xy_lims)
    end
    return xy_lims
end

function share_xy_axes!(
        plot_grid::Plots.Plot)
    layout_dims = size(plot_grid.layout.grid)
    nrows = layout_dims[1]
    ncols = layout_dims[2]
    for col_i in 1:ncols
        for row_i in 1:nrows
            if col_i != 1
                plt_index = ((row_i - 1) * ncols) + col_i
                plot!(plot_grid[plt_index], yformatter=_->"", ylabel = "")
            end
            if row_i != nrows
                plt_index = ((row_i - 1) * ncols) + col_i
                plot!(plot_grid[plt_index], xformatter=_->"", xlabel = "")
            end
        end
    end
    return nothing
end

function prep_for_violin(
        rng::MersenneTwister,
        v::Vector{Float64})::Vector{Float64}
    s = Set(v)
    if ! ((length(s) == 1) & (1.0 in s))
        return v
    end
    uniform_dist = Uniform(0.0, 1e-8)
    return v - rand(uniform_dist, length(v))
end

function main_cli()::Cint
    if ! Base.Filesystem.ispath(ProjectUtil.RESULTS_DIR)
        Base.Filesystem.mkdir(ProjectUtil.RESULTS_DIR)
    end

    # Set backend to GR
    # gr()
    write(stdout, "Plotting backend: $(backend())\n")

    brooks_gelman_1998_recommended_psrf = 1.2

    rng = Random.MersenneTwister()
    Random.seed!(rng, 39475839)

    results = SimResults()
    nrows = size(results.df)[1]

    root_node_probs = get_floats(results, true, true, true, false, :root_node_true_prob)
    vo_root_node_probs = get_floats(results, true, true, true, true, :root_node_true_prob)
    vln_root_node_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(prep_for_violin(rng, root_node_probs),
                 prep_for_violin(rng, vo_root_node_probs)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Root node",
            #= size = (600, 400), =#
    )

    node_789_probs = get_floats(results, true, true, true, false, :node_7_8_9_prob)
    vo_node_789_probs = get_floats(results, true, true, true, true, :node_7_8_9_prob)
    vln_node_789_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(prep_for_violin(rng, node_789_probs),
                 prep_for_violin(rng, vo_node_789_probs)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Node 789",
            #= size = (600, 400), =#
    )

    node_456_probs = get_floats(results, true, true, true, false, :node_4_5_6_prob)
    vo_node_456_probs = get_floats(results, true, true, true, true, :node_4_5_6_prob)
    vln_node_456_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(prep_for_violin(rng, node_456_probs),
                 prep_for_violin(rng, vo_node_456_probs)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Node 456",
            #= size = (600, 400), =#
    )

    node_12_3_probs = get_floats(results, true, true, true, false, :node_12_3_prob)
    vo_node_12_3_probs = get_floats(results, true, true, true, true, :node_12_3_prob)
    write(stdout, "$node_12_3_probs\n")
    write(stdout, "$vo_node_12_3_probs\n")
    vln_node_12_3_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(prep_for_violin(rng, node_12_3_probs),
                 prep_for_violin(rng, vo_node_12_3_probs)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Node 12-3",
            #= size = (600, 400), =#
    )

    height_12_789_probs = get_floats(results, true, true, true, false, :height_12_789_prob)
    vo_height_12_789_probs = get_floats(results, true, true, true, true, :height_12_789_prob)
    vln_height_12_789_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(prep_for_violin(rng, height_12_789_probs),
                 prep_for_violin(rng, vo_height_12_789_probs)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Height 12-789",
            #= size = (600, 400), =#
    )

    height_123_456_probs = get_floats(results, true, true, true, false, :height_123_456_prob)
    vo_height_123_456_probs = get_floats(results, true, true, true, true, :height_123_456_prob)
    vln_height_123_456_probs = StatsPlots.violin(["All sites" "SNPs"],
            hcat(prep_for_violin(rng, height_123_456_probs),
                 prep_for_violin(rng, vo_height_123_456_probs)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Height 123-456",
            #= size = (600, 400), =#
    )

    split_12_probs = get_floats(results, true, true, true, false, :split_12_prob)
    vo_split_12_probs = get_floats(results, true, true, true, true, :split_12_prob)
    vln_split_12_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(prep_for_violin(rng, split_12_probs),
                 prep_for_violin(rng, vo_split_12_probs)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Split 12",
            #= size = (600, 400), =#
    )

    true_topo_probs = get_floats(results, true, true, true, false, :topo_true_prob)
    vo_true_topo_probs = get_floats(results, true, true, true, true, :topo_true_prob)
    vln_true_topo_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(prep_for_violin(rng, true_topo_probs),
                 prep_for_violin(rng, vo_true_topo_probs)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "True topology",
            #= size = (600, 400), =#
    )

    vln_grid = Plots.plot(
            vln_root_node_probs,
            vln_node_789_probs,
            vln_node_456_probs,
            vln_node_12_3_probs,
            vln_height_12_789_probs,
            vln_height_123_456_probs,
            vln_split_12_probs,
            vln_true_topo_probs,
            layout = (2, 4),
            legend = false,
            size = (800, 400),
            link = :all, # :none, :x, :y, :both, :all
    )
    Plots.ylabel!(vln_grid, "Posterior probability")
    share_xy_axes!(vln_grid)

    Plots.savefig(vln_grid, joinpath(ProjectUtil.RESULTS_DIR, "fixed-gen-true-tree-probs.pdf"))


    gen_max_456_subsplit_prob = get_floats(results, true, true, true, false, :max_456_subsplit_prob)
    vo_gen_max_456_subsplit_prob = get_floats(results, true, true, true, true, :max_456_subsplit_prob)
    bif_max_456_subsplit_prob = get_floats(results, true, true, false, false, :max_456_subsplit_prob)
    vo_bif_max_456_subsplit_prob = get_floats(results, true, true, false, true, :max_456_subsplit_prob)
    vln_max_456_subsplit_probs = StatsPlots.violin(
            ["gen-all" "gen-SNPs" "bif-all" "bif-SNPs"],
            hcat(prep_for_violin(rng, gen_max_456_subsplit_prob),
                 prep_for_violin(rng, vo_gen_max_456_subsplit_prob),
                 prep_for_violin(rng, bif_max_456_subsplit_prob),
                 prep_for_violin(rng, vo_bif_max_456_subsplit_prob)),
            ylims = (0.0, 1.0),
            ylabel = "Max posterior probability",
            legend = false,
    )
    Plots.savefig(vln_max_456_subsplit_probs,
            joinpath(ProjectUtil.RESULTS_DIR, "fixed-gen-max-456-subsplit-probs.pdf"))


    gen_max_789_subsplit_prob = get_floats(results, true, true, true, false, :max_789_subsplit_prob)
    vo_gen_max_789_subsplit_prob = get_floats(results, true, true, true, true, :max_789_subsplit_prob)
    bif_max_789_subsplit_prob = get_floats(results, true, true, false, false, :max_789_subsplit_prob)
    vo_bif_max_789_subsplit_prob = get_floats(results, true, true, false, true, :max_789_subsplit_prob)
    vln_max_789_subsplit_probs = StatsPlots.violin(
            ["gen-all" "gen-SNPs" "bif-all" "bif-SNPs"],
            hcat(prep_for_violin(rng, gen_max_789_subsplit_prob),
                 prep_for_violin(rng, vo_gen_max_789_subsplit_prob),
                 prep_for_violin(rng, bif_max_789_subsplit_prob),
                 prep_for_violin(rng, vo_bif_max_789_subsplit_prob)),
            ylims = (0.0, 1.0),
            ylabel = "Max posterior probability",
            legend = false,
    )
    Plots.savefig(vln_max_789_subsplit_probs,
            joinpath(ProjectUtil.RESULTS_DIR, "fixed-gen-max-789-subsplit-probs.pdf"))


    bif_gen_root_node_probs = get_floats(results, true, false, true, false, :root_node_true_prob)
    vo_bif_gen_root_node_probs = get_floats(results, true, false, true, true, :root_node_true_prob)
    bif_bif_root_node_probs = get_floats(results, true, false, false, false, :root_node_true_prob)
    vo_bif_bif_root_node_probs = get_floats(results, true, false, false, true, :root_node_true_prob)

    bif_gen_true_split_prob_mean = get_floats(results, true, false, true, false, :true_split_prob_mean)
    vo_bif_gen_true_split_prob_mean = get_floats(results, true, false, true, true, :true_split_prob_mean)
    bif_bif_true_split_prob_mean = get_floats(results, true, false, false, false, :true_split_prob_mean)
    vo_bif_bif_true_split_prob_mean = get_floats(results, true, false, false, true, :true_split_prob_mean)

    bif_gen_ht_12_78_prob = get_floats(results, true, false, true, false, :height_12_78_prob)
    vo_bif_gen_ht_12_78_prob = get_floats(results, true, false, true, true, :height_12_78_prob)

    bif_gen_ht_45_789_prob = get_floats(results, true, false, true, false, :height_45_789_prob)
    vo_bif_gen_ht_45_789_prob = get_floats(results, true, false, true, true, :height_45_789_prob)

    bif_gen_ht_456_789_prob = get_floats(results, true, false, true, false, :height_456_789_prob)
    vo_bif_gen_ht_456_789_prob = get_floats(results, true, false, true, true, :height_456_789_prob)

    bif_gen_nd_1_2_3_prob = get_floats(results, true, false, true, false, :node_1_2_3_prob)
    vo_bif_gen_nd_1_2_3_prob = get_floats(results, true, false, true, true, :node_1_2_3_prob)

    bif_gen_nd_4_5_6_prob = get_floats(results, true, false, true, false, :node_4_5_6_prob)
    vo_bif_gen_nd_4_5_6_prob = get_floats(results, true, false, true, true, :node_4_5_6_prob)

    bif_gen_nd_123_456_789_prob = get_floats(results, true, false, true, false, :node_123_456_789_prob)
    vo_bif_gen_nd_123_456_789_prob = get_floats(results, true, false, true, true, :node_123_456_789_prob)


    vln_bif_root_node_probs = StatsPlots.violin(
            ["gen-all" "gen-SNPs" "bif-all" "bif-SNPs"],
            hcat(prep_for_violin(rng, bif_gen_root_node_probs),
                 prep_for_violin(rng, vo_bif_gen_root_node_probs),
                 prep_for_violin(rng, bif_bif_root_node_probs),
                 prep_for_violin(rng, vo_bif_bif_root_node_probs)),
            ylims = (0.0, 1.0),
            ylabel = "Posterior probability",
            legend = false,
    )
    Plots.savefig(vln_bif_root_node_probs,
            joinpath(ProjectUtil.RESULTS_DIR, "fixed-bif-root-node-probs.pdf"))

    vln_bif_true_split_prob_means = StatsPlots.violin(
            ["gen-all" "gen-SNPs" "bif-all" "bif-SNPs"],
            hcat(prep_for_violin(rng, bif_gen_true_split_prob_mean),
                 prep_for_violin(rng, vo_bif_gen_true_split_prob_mean),
                 prep_for_violin(rng, bif_bif_true_split_prob_mean),
                 prep_for_violin(rng, vo_bif_bif_true_split_prob_mean)),
            ylims = (0.0, 1.0),
            legend = false,
    )
    Plots.savefig(vln_bif_true_split_prob_means,
            joinpath(ProjectUtil.RESULTS_DIR, "fixed-bif-true-split-prob-means.pdf"))

    vln_bif_gen_ht_12_78_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(prep_for_violin(rng, bif_gen_ht_12_78_prob),
                 prep_for_violin(rng, vo_bif_gen_ht_12_78_prob)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Height 12-78",
    )
    vln_bif_gen_ht_45_789_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(prep_for_violin(rng, bif_gen_ht_45_789_prob),
                 prep_for_violin(rng, vo_bif_gen_ht_45_789_prob)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Height 45-789",
    )
    vln_bif_gen_ht_456_789_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(prep_for_violin(rng, bif_gen_ht_456_789_prob),
                 prep_for_violin(rng, vo_bif_gen_ht_456_789_prob)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Height 456-789",
    )
    vln_bif_gen_nd_1_2_3_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(prep_for_violin(rng, bif_gen_nd_1_2_3_prob),
                 prep_for_violin(rng, vo_bif_gen_nd_1_2_3_prob)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Node 1-2-3",
    )
    vln_bif_gen_nd_4_5_6_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(prep_for_violin(rng, bif_gen_nd_4_5_6_prob),
                 prep_for_violin(rng, vo_bif_gen_nd_4_5_6_prob)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Node 4-5-6",
    )
    vln_bif_gen_nd_123_456_789_probs = StatsPlots.violin(
            ["All sites" "SNPs"],
            hcat(prep_for_violin(rng, bif_gen_nd_123_456_789_prob),
                 prep_for_violin(rng, vo_bif_gen_nd_123_456_789_prob)),
            ylims = (0.0, 1.0),
            legend = false,
            title = "Node 123-456-789",
    )
    vln_bif_gen_wrong_probs_grid = Plots.plot(
            vln_bif_gen_ht_12_78_probs,
            vln_bif_gen_ht_45_789_probs,
            vln_bif_gen_ht_456_789_probs,
            vln_bif_gen_nd_1_2_3_probs,
            vln_bif_gen_nd_4_5_6_probs,
            vln_bif_gen_nd_123_456_789_probs,
            layout = (2, 3),
            legend = false,
            size = (800, 350),
            link = :all,
    )
    Plots.ylabel!(vln_bif_gen_wrong_probs_grid, "Posterior probability")
    share_xy_axes!(vln_bif_gen_wrong_probs_grid)

    Plots.savefig(vln_bif_gen_wrong_probs_grid, joinpath(ProjectUtil.RESULTS_DIR, "fixed-bif-gen-wrong-probs.pdf"))


    gen_gen_tree_len = get_true_v_est_plot(
            results,
            false,
            true,
            true,
            false,
            "tree_length")
    vo_gen_gen_tree_len = get_true_v_est_plot(
            results,
            false,
            true,
            true,
            true,
            "tree_length")
    gen_bif_tree_len = get_true_v_est_plot(
            results,
            false,
            true,
            false,
            false,
            "tree_length")
    vo_gen_bif_tree_len = get_true_v_est_plot(
            results,
            false,
            true,
            false,
            true,
            "tree_length")

    bif_gen_tree_len = get_true_v_est_plot(
            results,
            false,
            false,
            true,
            false,
            "tree_length")
    vo_bif_gen_tree_len = get_true_v_est_plot(
            results,
            false,
            false,
            true,
            true,
            "tree_length")
    bif_bif_tree_len = get_true_v_est_plot(
            results,
            false,
            false,
            false,
            false,
            "tree_length")
    vo_bif_bif_tree_len = get_true_v_est_plot(
            results,
            false,
            false,
            false,
            true,
            "tree_length")

    #= Plots.savefig(gen_gen_tree_len, joinpath(ProjectUtil.RESULTS_DIR, "unfixed-gen-gen-tree-length-allsites.pdf")) =#
    #= Plots.savefig(bif_bif_tree_len, joinpath(ProjectUtil.RESULTS_DIR, "unfixed-bif-bif-tree-length-allsites.pdf")) =#

    #= write(stdout, "\n\n$(plt_tree_len.subplots[1][:xaxis][:lims])\n\n") =#
    #= write(stdout, "\n\n$(plt_tree_len.subplots[1][:yaxis][:lims])\n\n") =#

    #= Plots.savefig(plt_tree_len, joinpath(ProjectUtil.RESULTS_DIR, "tree-length-unfixed-gen-gen-all.pdf")) =#

    #= xy_lims = share_xy_axes!(plt_tree_len, plt_tree_len2) =#

    tree_len_grid = Plots.plot(
            gen_gen_tree_len,
            vo_gen_gen_tree_len,
            gen_bif_tree_len,
            vo_gen_bif_tree_len,
            bif_gen_tree_len,
            vo_bif_gen_tree_len,
            bif_bif_tree_len,
            vo_bif_bif_tree_len,
            layout = (4, 2),
            legend = false,
            size = (700, 1200),
            link = :all,
            xlabel = "True tree length",
            ylabel = "Posterior mean tree length",
            )
    xy_lims = share_xy_limits!(tree_len_grid)
    share_xy_axes!(tree_len_grid)
    Plots.savefig(tree_len_grid, joinpath(ProjectUtil.RESULTS_DIR, "unfixed-tree-length.pdf"))


    # Plot ASDSF
    fixed_gen_gen_asdsf = get_floats(results, true, true, true, false, :sdsf_mean)
    vo_fixed_gen_gen_asdsf = get_floats(results, true, true, true, true, :sdsf_mean)
    fixed_gen_bif_asdsf = get_floats(results, true, true, false, false, :sdsf_mean)
    vo_fixed_gen_bif_asdsf = get_floats(results, true, true, false, true, :sdsf_mean)
    fixed_bif_gen_asdsf = get_floats(results, true, false, true, false, :sdsf_mean)
    vo_fixed_bif_gen_asdsf = get_floats(results, true, false, true, true, :sdsf_mean)
    fixed_bif_bif_asdsf = get_floats(results, true, false, false, false, :sdsf_mean)
    vo_fixed_bif_bif_asdsf = get_floats(results, true, false, false, true, :sdsf_mean)

    write(stdout, "$fixed_bif_gen_asdsf\n")
    write(stdout, "$fixed_bif_bif_asdsf\n")

    vln_fixed_gen_asdsf = StatsPlots.violin(
            ["gen-all" "gen-SNPs" "bif-all" "bif-SNPs"],
            hcat(prep_for_violin(rng, fixed_gen_gen_asdsf),
                 prep_for_violin(rng, vo_fixed_gen_gen_asdsf),
                 prep_for_violin(rng, fixed_gen_bif_asdsf),
                 prep_for_violin(rng, vo_fixed_gen_bif_asdsf)),
            legend = false,
            title = "generalized",
    )
    vln_fixed_bif_asdsf = StatsPlots.violin(
            ["gen-all" "gen-SNPs" "bif-all" "bif-SNPs"],
            hcat(prep_for_violin(rng, fixed_bif_gen_asdsf),
                 prep_for_violin(rng, vo_fixed_bif_gen_asdsf),
                 prep_for_violin(rng, fixed_bif_bif_asdsf),
                 prep_for_violin(rng, vo_fixed_bif_bif_asdsf)),
            legend = false,
            title = "bifurcating",
    )
    vln_fixed_asdsf_grid = Plots.plot(
            vln_fixed_gen_asdsf,
            vln_fixed_bif_asdsf,
            layout = (1, 2),
            legend = false,
            size = (800, 300),
            link = :all,
    )
    Plots.ylabel!(vln_fixed_asdsf_grid, "ASDSF")
    share_xy_axes!(vln_fixed_asdsf_grid)
    Plots.savefig(vln_fixed_asdsf_grid,
            joinpath(ProjectUtil.RESULTS_DIR, "fixed-asdsf.pdf"))

    unfixed_gen_gen_asdsf = get_floats(results, false, true, true, false, :sdsf_mean)
    vo_unfixed_gen_gen_asdsf = get_floats(results, false, true, true, true, :sdsf_mean)
    unfixed_gen_bif_asdsf = get_floats(results, false, true, false, false, :sdsf_mean)
    vo_unfixed_gen_bif_asdsf = get_floats(results, false, true, false, true, :sdsf_mean)
    unfixed_bif_gen_asdsf = get_floats(results, false, false, true, false, :sdsf_mean)
    vo_unfixed_bif_gen_asdsf = get_floats(results, false, false, true, true, :sdsf_mean)
    unfixed_bif_bif_asdsf = get_floats(results, false, false, false, false, :sdsf_mean)
    vo_unfixed_bif_bif_asdsf = get_floats(results, false, false, false, true, :sdsf_mean)

    vln_unfixed_gen_asdsf = StatsPlots.violin(
            ["gen-all" "gen-SNPs" "bif-all" "bif-SNPs"],
            hcat(prep_for_violin(rng, unfixed_gen_gen_asdsf),
                 prep_for_violin(rng, vo_unfixed_gen_gen_asdsf),
                 prep_for_violin(rng, unfixed_gen_bif_asdsf),
                 prep_for_violin(rng, vo_unfixed_gen_bif_asdsf)),
            legend = false,
            title = "generalized",
    )
    vln_unfixed_bif_asdsf = StatsPlots.violin(
            ["gen-all" "gen-SNPs" "bif-all" "bif-SNPs"],
            hcat(prep_for_violin(rng, unfixed_bif_gen_asdsf),
                 prep_for_violin(rng, vo_unfixed_bif_gen_asdsf),
                 prep_for_violin(rng, unfixed_bif_bif_asdsf),
                 prep_for_violin(rng, vo_unfixed_bif_bif_asdsf)),
            legend = false,
            title = "bifurcating",
    )
    vln_unfixed_asdsf_grid = Plots.plot(
            vln_unfixed_gen_asdsf,
            vln_unfixed_bif_asdsf,
            layout = (1, 2),
            legend = false,
            size = (800, 300),
            link = :all,
    )
    Plots.ylabel!(vln_unfixed_asdsf_grid, "ASDSF")
    share_xy_axes!(vln_unfixed_asdsf_grid)
    Plots.savefig(vln_unfixed_asdsf_grid,
            joinpath(ProjectUtil.RESULTS_DIR, "unfixed-asdsf.pdf"))

    return 0
end

main_cli()
