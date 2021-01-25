#! /usr/bin/env julia

using Random
using Distributions
using Statistics
using DataDeps
using DataFrames
using CSV
using Plots
using StatsPlots
using HypothesisTests
using LaTeXStrings
# using GR

include("project_util.jl")

dark_blue_col = RGB(2/255, 75/255, 120/255)
l_dark_blue_col = RGB(3/255, 138/255, 221/255)
comp_blue_col= RGB(2/255, 219/255, 240/255)
dark_orange_col = RGB(186/255, 88/255, 0/255)
l_dark_orange_col = RGB(255/255, 135/255, 31/255)
comp_orange_col = RGB(255/255, 181/255, 91/255)
highlight_col = "red"

gen_col = dark_blue_col
vo_gen_col = comp_blue_col
#= vo_gen_col = dark_blue_col =#
bif_col = dark_orange_col 
vo_bif_col = comp_orange_col 
#= vo_bif_col = dark_orange_col =# 

marker_alpha = 0.8
fill_alpha = 0.5

gen_marker_alpha = marker_alpha
#= vo_gen_marker_alpha = marker_alpha - 0.3 =#
vo_gen_marker_alpha = marker_alpha
gen_fill_alpha = fill_alpha
#= vo_gen_fill_alpha = fill_alpha - 0.3 =#
vo_gen_fill_alpha = fill_alpha

bif_marker_alpha = marker_alpha
vo_bif_marker_alpha = marker_alpha
bif_fill_alpha = fill_alpha
vo_bif_fill_alpha = fill_alpha

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
                for locus_size in [1, 100]
                    locus_conditions = [""]
                    if locus_size > 1
                        locus_conditions = ["locus-$locus_size-linked-", "locus-$locus_size-unlinked-"]
                    end
                    for locus_condition in locus_conditions
                        sim_dir_name = ""
                        if is_fixed
                            sim_dir_name = "$(locus_condition)species-9-genomes-2-$sim_condition-tree-provided-fixed"
                        else
                            sim_dir_name = "$(locus_condition)species-9-genomes-2-$sim_condition-tree-random-unfixed"
                        end
                        sim_dir = joinpath(ProjectUtil.SIM_DIR, sim_dir_name)
                        for batch_dir in ProjectUtil.batch_dir_iter(sim_dir)
                            for analysis_condition in ["bifurcating", "generalized"]
                                var_only_bools = [true, false]
                                if endswith(locus_condition, "-linked-")
                                    var_only_bools = [false]
                                elseif endswith(locus_condition, "-unlinked-")
                                    var_only_bools = [true]
                                end
                                for is_var_only in var_only_bools
                                    result_file_name = "results-species-9-genomes-2-$analysis_condition-tree-random.tsv.gz"
                                    if is_var_only & (locus_condition != "locus-100-unlinked-")
                                        result_file_name = "results-var-only-species-9-genomes-2-$analysis_condition-tree-random.tsv.gz"
                                    end
                                    result_path = joinpath(batch_dir, result_file_name)
                                    df = ProjectUtil.get_data_frame_from_gz([result_path])
                                    nrows = size(df)[1]
                                    df[!, "sim_fixed"] = repeat([is_fixed], nrows)
                                    df[!, "sim_generalized"] = repeat([sim_condition == "generalized"], nrows)
                                    df[!, "analysis_generalized"] = repeat([analysis_condition == "generalized"], nrows)
                                    df[!, "only_var_sites"] = repeat([is_var_only], nrows)
                                    df[!, "locus_size"] = repeat([locus_size], nrows)
                                    results = vcat(results, df)
                                end
                            end
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
    function SimResults(data_frame::DataFrame)
        return new(data_frame)
    end
end

function get_rows(r::SimResults,
        fixed::Bool,
        sim_generalized::Bool,
        analysis_generalized::Bool,
        var_only::Bool,
        locus_size::Int = 1
        )::DataFrame
    return r.df[(r.df[!, :sim_fixed] .== fixed) .& (r.df[!, :sim_generalized] .== sim_generalized) .& (r.df[!, :analysis_generalized] .== analysis_generalized) .& (r.df[!, :only_var_sites] .== var_only) .& (r.df[!, :locus_size] .== locus_size), :]
end

function get_floats(r::SimResults,
        fixed::Bool,
        sim_generalized::Bool,
        analysis_generalized::Bool,
        var_only::Bool,
        col_header::Union{Symbol, AbstractString},
        locus_size::Int = 1
        )::Vector{Float64}
    return r.df[(r.df[!, :sim_fixed] .== fixed) .& (r.df[!, :sim_generalized] .== sim_generalized) .& (r.df[!, :analysis_generalized] .== analysis_generalized) .& (r.df[!, :only_var_sites] .== var_only) .& (r.df[!, :locus_size] .== locus_size), :][!, col_header]
end

function get_ints(r::SimResults,
        fixed::Bool,
        sim_generalized::Bool,
        analysis_generalized::Bool,
        var_only::Bool,
        col_header::Union{Symbol, AbstractString},
        locus_size::Int = 1
        )::Vector{Int64}
    return r.df[(r.df[!, :sim_fixed] .== fixed) .& (r.df[!, :sim_generalized] .== sim_generalized) .& (r.df[!, :analysis_generalized] .== analysis_generalized) .& (r.df[!, :only_var_sites] .== var_only) .& (r.df[!, :locus_size] .== locus_size), :][!, col_header]
end

function get_bools(r::SimResults,
        fixed::Bool,
        sim_generalized::Bool,
        analysis_generalized::Bool,
        var_only::Bool,
        col_header::Union{Symbol, AbstractString},
        locus_size::Int = 1
        )::Vector{Bool}
    return r.df[(r.df[!, :sim_fixed] .== fixed) .& (r.df[!, :sim_generalized] .== sim_generalized) .& (r.df[!, :analysis_generalized] .== analysis_generalized) .& (r.df[!, :only_var_sites] .== var_only) .& (r.df[!, :locus_size] .== locus_size), :][!, col_header]
end

function get_groups_by_y(
        y,
        y_lower,
        y_upper,
        labels,
        colors,
        alphas;
        y_buffer = 0.02,
        error_bar_alpha = 0.2,
        show_group_means = true,
        group_mean_line_alpha = 0.8,
        group_mean_line_width = 2.0,
        show_labels_on_x = false,
        x_label_size = 8)::Plots.Plot
    @assert ndims(y) == 2
    @assert size(y) == size(y_lower) 
    @assert size(y) == size(y_upper) 
    @assert size(y)[2] == size(colors)[2]
    @assert size(y)[2] == size(labels)[2]
    nsamples = size(y)[1]
    ngroups = size(y)[2]
    x = Array{Int64}(undef, size(y)...)
    x_val = 0
    for i in eachindex(x)
        x_val += 1
        x[i] = x_val
    end
    y_extremes = (minimum(y_lower), maximum(y_upper))
    y_buf = (y_extremes[2] - y_extremes[1]) * y_buffer
    y_limits = [y_extremes[1] - y_buf, y_extremes[2] + y_buf]
    x_limits = [0, x_val + 1]
    plt = Plots.plot(
            x[:,1],
            y[:,1],
            yerror = (y[:,1] - y_lower[:,1], y_upper[:,1] - y[:,1]),
            seriestype = :scatter,
            legend = false,
            xlims = x_limits,
            ylims = y_limits,
            label = labels[:,1],
            markercolor = colors[:,1],
            markeralpha = alphas[:,1],
            markerstrokecolor = colors[:,1],
            markerstrokealpha = error_bar_alpha,
            linecolor = colors[:,1],
            linealpha = error_bar_alpha)
    if show_group_means
        y_mean = Statistics.mean(y[:,1])
        #= write(stdout, "$y_mean\n") =#
        Plots.plot!(plt,
                [x[1, 1],  x[nsamples, 1]],
                [y_mean,  y_mean],
                seriestype = :line,
                legend = false,
                linecolor = colors[:,1],
                linealpha = group_mean_line_alpha,
                linewidth = group_mean_line_width)
    end
    if show_labels_on_x
        x_position = Statistics.mean(x[:,1])
        y_position = y_limits[1]
        annotate!(plt, x_position, y_position,
                  text(labels[1,1],
                       :center,
                       :top,
                       x_label_size),
                  annotation_clip = false)
    end
    for i in 2:ngroups
        Plots.plot!(plt,
                x[:,i],
                y[:,i],
                yerror = (y[:,i] - y_lower[:,i], y_upper[:,i] - y[:,i]),
                seriestype = :scatter,
                legend = false,
                label = labels[:,i],
                markercolor = colors[:,i],
                markeralpha = alphas[:,i],
                markerstrokecolor = colors[:,i],
                markerstrokealpha = error_bar_alpha,
                linecolor = colors[:,i],
                linealpha = error_bar_alpha)
        if show_group_means
            y_mean = Statistics.mean(y[:,i])
            #= write(stdout, "$y_mean\n") =#
            Plots.plot!(plt,
                    [x[1, i],  x[nsamples, i]],
                    [y_mean,  y_mean],
                    seriestype = :line,
                    legend = false,
                    linecolor = colors[:,i],
                    linealpha = group_mean_line_alpha,
                    linewidth = group_mean_line_width)
        end
        if show_labels_on_x
            x_position = Statistics.mean(x[:,i])
            y_position = y_limits[1]
            annotate!(plt, x_position, y_position,
                      text(labels[1,i],
                           :center,
                           :top,
                           x_label_size),
                      annotation_clip = false)
        end
    end
    Plots.xaxis!(plt, false)
    return plt
end

function get_violin_plot(
        values;
        xlabels,
        fill_colors,
        marker_colors,
        fill_alphas,
        marker_alphas,
        include_dots = true)::Plots.Plot
    vln = StatsPlots.violin(
            xlabels,
            values,
            legend = false,
            fillcolor = fill_colors,
            linecolor = fill_colors,
            fillalpha = fill_alphas,
    )
    if include_dots
        StatsPlots.dotplot!(vln,
                xlabels,
                left_values,
                legend = false,
                markercolor = marker_colors,
                markeralpha = marker_alphas
        )
    end
    return vln
end

function get_split_violin_plot(
        left_values,
        right_values;
        xlabels,
        left_fill_colors,
        left_marker_colors,
        left_fill_alphas,
        left_marker_alphas,
        left_labels,
        right_fill_colors,
        right_marker_colors,
        right_fill_alphas,
        right_marker_alphas,
        right_labels,
        legend = false,
        dot_legend = false)::Plots.Plot
    vln = StatsPlots.violin(
            xlabels,
            left_values,
            legend = legend,
            side = :left,
            fillcolor = left_fill_colors,
            fillalpha = left_fill_alphas,
            label = left_labels
    )
    StatsPlots.dotplot!(vln,
            xlabels,
            left_values,
            legend = dot_legend,
            side = :left,
            markercolor = left_marker_colors,
            markeralpha = left_marker_alphas,
            label = left_labels
    )
    StatsPlots.violin!(vln,
            xlabels,
            right_values,
            legend = legend,
            side = :right,
            fillcolor = right_fill_colors,
            fillalpha = right_fill_alphas,
            label = right_labels
    )
    StatsPlots.dotplot!(vln,
            xlabels,
            right_values,
            legend = dot_legend,
            side = :right,
            markercolor = right_marker_colors,
            markeralpha = right_marker_alphas,
            label = right_labels
    )
    return vln
end


function get_true_v_est_plot(r::SimResults,
        fixed::Bool,
        sim_generalized::Bool,
        analysis_generalized::Bool,
        var_only::Bool,
        col_header_prefix::String,
        locus_size::Int = 1,
        use_median::Bool = false,
        axis_limit_buffer::AbstractFloat = 0.05,
        psrf_threshold::AbstractFloat = -1.0,
        ess_threshold::AbstractFloat = -1.0,
        highlight_color = "red"
        )::Plots.Plot
    x = get_floats(r,
            fixed,
            sim_generalized,
            analysis_generalized,
            var_only,
            col_header_prefix * "_true",
            locus_size)
    if use_median
        y = get_floats(r,
                fixed,
                sim_generalized,
                analysis_generalized,
                var_only,
                col_header_prefix * "_median",
                locus_size)
    else
        y = get_floats(r,
                fixed,
                sim_generalized,
                analysis_generalized,
                var_only,
                col_header_prefix * "_mean",
                locus_size)
    end
    y_lower = get_floats(r,
            fixed,
            sim_generalized,
            analysis_generalized,
            var_only,
            col_header_prefix * "_eti_95_lower",
            locus_size)
    y_upper = get_floats(r,
            fixed,
            sim_generalized,
            analysis_generalized,
            var_only,
            col_header_prefix * "_eti_95_upper",
            locus_size)

    extremes = (minimum((minimum(x), minimum(y))),
                maximum((maximum(x), maximum(y))))
    xy_range = extremes[2] - extremes[1]
    xy_buffer = xy_range * axis_limit_buffer
    axis_limits = [extremes[1] - xy_buffer,
                   extremes[2] + xy_buffer]
    identity_line = [extremes[1] - (xy_buffer * 20.0)
                     extremes[2] + (xy_buffer * 20.0)]

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
    if (psrf_threshold > 0) | (ess_threshold > 0)
        highlight_x = Vector{Float64}()
        highlight_y = Vector{Float64}()
        highlight_y_lower = Vector{Float64}()
        highlight_y_upper = Vector{Float64}()
        psrf = get_floats(r,
                fixed,
                sim_generalized,
                analysis_generalized,
                var_only,
                col_header_prefix * "_psrf",
                locus_size)
        ess = get_floats(r,
                fixed,
                sim_generalized,
                analysis_generalized,
                var_only,
                col_header_prefix * "_ess",
                locus_size)
        for i in 1:length(x)
            if (((psrf_threshold > 0) & (psrf[i] > psrf_threshold)) |
                    ((ess_threshold > 0) & (ess[i] < ess_threshold)))
                append!(highlight_x, x[i])
                append!(highlight_y, y[i])
                append!(highlight_y_lower, y_lower[i])
                append!(highlight_y_upper, y_upper[i])
            end
        end
        if length(highlight_x) > 0
            Plots.plot!(plt,
                    highlight_x,
                    highlight_y,
                    yerror = (highlight_y - highlight_y_lower, highlight_y_upper - highlight_y),
                    seriestype = :scatter,
                    markercolor = highlight_color)
        end
    end
    return plt
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

# A hack to fix error with StatsPlots.violin that causes it to choke when the
# values are all equal.
function prep_for_violin(
        rng::MersenneTwister,
        v::Vector{Float64},
        max_jiggle::Float64 = 1e-3)::Vector{Float64}
    s = Set(v)
    if (length(s) == 1)
        uniform_dist = Uniform(0.0, max_jiggle)
        if (1.0 in s)
            return 1.0 .- rand(uniform_dist, 10000)
        end
        if (0.0 in s)
            return 0.0 .+ rand(uniform_dist, 10000)
        end
    end
    return v
end

function relative_xy(plt::Plots.Plot, x, y)
        x_lims = Plots.xlims(plt)
        y_lims = Plots.ylims(plt)
        coords = x_lims[1] + (x * (x_lims[2] - x_lims[1])),
                 y_lims[1] + (y * (y_lims[2] - y_lims[1]))
        return coords
end


function main_cli()::Cint
    if ! Base.Filesystem.ispath(ProjectUtil.RESULTS_DIR)
        Base.Filesystem.mkdir(ProjectUtil.RESULTS_DIR)
    end

    # Set backend to GR
    # gr()
    write(stdout, "Plotting backend: $(backend())\n")

    brooks_gelman_1998_recommended_psrf = 1.2
    ess_threshold = 200.0

    violin_probability_jiggle = 1.0e-3
    violin_asdsf_jiggle = 1e-8

    prob_bins = collect(0.0:0.1:1.0)
    prob_bins[length(prob_bins)] += 0.0001

    rng = Random.MersenneTwister()
    Random.seed!(rng, 39475839)

    results = SimResults()
    nrows = size(results.df)[1]

    locus_sizes = [1, 100]

    for locus_size in locus_sizes
        locus_prefix = ""
        if locus_size > 1
            locus_prefix = "locus-$locus_size-"
        end

        sum_stats = DataFrame()
        bools = [true, false]
        for sim_fixed in bools
            for sim_generalized in bools
                for analysis_generalized in bools
                    for only_var_sites in bools
                        tree_len_perc = get_floats(results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :tree_length_true_percentile,
                                locus_size)
                        tree_len_cover = Statistics.mean(
                                0.025 .<= tree_len_perc .<= 0.975)
                        root_height_perc = get_floats(results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :root_height_true_percentile,
                                locus_size)
                        root_height_cover = Statistics.mean(
                                0.025 .<= root_height_perc .<= 0.975)
                        root_pop_size_perc = get_floats(results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :root_pop_size_true_percentile,
                                locus_size)
                        root_pop_size_cover = Statistics.mean(
                                0.025 .<= root_pop_size_perc .<= 0.975)

                        tree_cred_level = get_floats(results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :topo_true_cred_level,
                                locus_size)
                        tree_cover = Statistics.mean(tree_cred_level .>= 0.95)

                        freq_tree_correct = Statistics.mean(get_bools(results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :topo_true_is_map,
                                locus_size))

                        num_heights_cred_level = get_floats(results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :num_heights_true_cred_level,
                                locus_size)
                        num_heights_cover = Statistics.mean(num_heights_cred_level .>= 0.95)

                        freq_num_heights_correct = Statistics.mean(
                                get_ints(results,
                                        sim_fixed,
                                        sim_generalized,
                                        analysis_generalized,
                                        only_var_sites,
                                        :num_heights_true,
                                        locus_size) .==
                                get_ints(results,
                                        sim_fixed,
                                        sim_generalized,
                                        analysis_generalized,
                                        only_var_sites,
                                        :num_heights_map,
                                        locus_size))

                        row = DataFrame(
                                sim_fixed = sim_fixed,
                                sim_generalized = sim_generalized,
                                analysis_generalized = analysis_generalized,
                                only_var_sites = only_var_sites,
                                locus_size = locus_size,
                                coverage_tree_length = tree_len_cover,
                                coverage_root_height = root_height_cover,
                                coverage_root_pop_size = root_pop_size_cover,
                                coverage_tree = tree_cover,
                                coverage_num_heights = num_heights_cover,
                                freq_correct_tree = freq_tree_correct,
                                freq_correct_num_heights = freq_num_heights_correct,
                                )
                        sum_stats = vcat(sum_stats, row)
                    end
                end
            end
        end

        sum_results = SimResults(sum_stats)
        for sim_fixed in bools
            for sim_generalized in bools
                for analysis_generalized in bools
                    for only_var_sites in bools
                        write(stdout, "\n\n")
                        write(stdout, "$(sim_fixed ? "Fixed" : "Free")")
                        write(stdout, " $(sim_generalized ? "gen" : "bif")")
                        write(stdout, " $(analysis_generalized ? "gen" : "bif")")
                        write(stdout, " $(only_var_sites ? "var-only" : "all-sites")")
                        write(stdout, " $(locus_size)\n")
                        write(stdout, "tree-length-coverage: $(get_floats(sum_results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :coverage_tree_length,
                                locus_size))\n")
                        write(stdout, "root-height-coverage: $(get_floats(sum_results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :coverage_root_height,
                                locus_size))\n")
                        write(stdout, "root-pop-size-coverage: $(get_floats(sum_results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :coverage_root_pop_size,
                                locus_size))\n")
                        write(stdout, "topo-coverage: $(get_floats(sum_results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :coverage_tree,
                                locus_size))\n")
                        write(stdout, "nheights-coverage: $(get_floats(sum_results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :coverage_num_heights,
                                locus_size))\n")
                        write(stdout, "freq-tree-correct: $(get_floats(sum_results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :freq_correct_tree,
                                locus_size))\n")
                        write(stdout, "freq-nheights-correct: $(get_floats(sum_results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :freq_correct_num_heights,
                                locus_size))\n")
                    end
                end
            end
        end


        root_node_probs = get_floats(results, true, true, true, false, :root_node_true_prob, locus_size)
        vo_root_node_probs = get_floats(results, true, true, true, true, :root_node_true_prob, locus_size)
        node_456_probs = get_floats(results, true, true, true, false, :node_4_5_6_prob, locus_size)
        vo_node_456_probs = get_floats(results, true, true, true, true, :node_4_5_6_prob, locus_size)
        node_789_probs = get_floats(results, true, true, true, false, :node_7_8_9_prob, locus_size)
        vo_node_789_probs = get_floats(results, true, true, true, true, :node_7_8_9_prob, locus_size)
        node_12_3_probs = get_floats(results, true, true, true, false, :node_12_3_prob, locus_size)
        vo_node_12_3_probs = get_floats(results, true, true, true, true, :node_12_3_prob, locus_size)
        height_12_789_probs = get_floats(results, true, true, true, false, :height_12_789_prob, locus_size)
        vo_height_12_789_probs = get_floats(results, true, true, true, true, :height_12_789_prob, locus_size)
        height_123_456_probs = get_floats(results, true, true, true, false, :height_123_456_prob, locus_size)
        vo_height_123_456_probs = get_floats(results, true, true, true, true, :height_123_456_prob, locus_size)
        split_12_probs = get_floats(results, true, true, true, false, :split_12_prob, locus_size)
        vo_split_12_probs = get_floats(results, true, true, true, true, :split_12_prob, locus_size)
        true_topo_probs = get_floats(results, true, true, true, false, :topo_true_prob, locus_size)
        vo_true_topo_probs = get_floats(results, true, true, true, true, :topo_true_prob, locus_size)

        v = get_split_violin_plot(
                root_node_probs,
                vo_root_node_probs,
                xlabels = [ "" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (200, 200))
        Plots.ylims!(v, (-0.02, 1.02))
        #= Plots.ylabel!(v, "Posterior probability") =#
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-root-vln.pdf"))

        v = get_violin_plot(
                prep_for_violin(rng, root_node_probs, violin_probability_jiggle),
                xlabels = [ "" ],
                fill_colors = gen_col,
                marker_colors = gen_col,
                fill_alphas = gen_fill_alpha,
                marker_alphas = gen_marker_alpha,
                include_dots = false)
        Plots.plot!(v, size = (120, 120), xaxis = false, grid = :y)
        #= Plots.ylims!(v, (-0.02, 1.02)) =#
        Plots.ylims!(v, (0.0, 1.0))
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-root-all-sites-vln.pdf"))

        h = Plots.histogram(root_node_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, grid = :x, xticks = (0:0.5:1))
        Plots.savefig(h,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-root-all-sites-hist.pdf"))

        v = get_split_violin_plot(
                node_789_probs,
                vo_node_789_probs,
                xlabels = [ "" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (200, 200))
        Plots.ylims!(v, (-0.02, 1.02))
        #= Plots.ylabel!(v, "Posterior probability") =#
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-789-vln.pdf"))

        v = get_violin_plot(
                prep_for_violin(rng, node_789_probs, violin_probability_jiggle),
                xlabels = [ "" ],
                fill_colors = gen_col,
                marker_colors = gen_col,
                fill_alphas = gen_fill_alpha,
                marker_alphas = gen_marker_alpha,
                include_dots = false)
        Plots.plot!(v, size = (120, 120), xaxis = false, grid = :y)
        #= Plots.ylims!(v, (-0.02, 1.02)) =#
        Plots.ylims!(v, (0.0, 1.0))
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-789-all-sites-vln.pdf"))

        h = Plots.histogram(node_789_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, grid = :x, xticks = (0:0.5:1))
        Plots.savefig(h,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-789-all-sites-hist.pdf"))

        v = get_split_violin_plot(
                node_456_probs,
                vo_node_456_probs,
                xlabels = [ "" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (200, 200))
        Plots.ylims!(v, (-0.02, 1.02))
        #= Plots.ylabel!(v, "Posterior probability") =#
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-456-vln.pdf"))

        v = get_violin_plot(
                prep_for_violin(rng, node_456_probs, violin_probability_jiggle),
                xlabels = [ "" ],
                fill_colors = gen_col,
                marker_colors = gen_col,
                fill_alphas = gen_fill_alpha,
                marker_alphas = gen_marker_alpha,
                include_dots = false)
        Plots.plot!(v, size = (120, 120), xaxis = false, grid = :y)
        #= Plots.ylims!(v, (-0.02, 1.02)) =#
        Plots.ylims!(v, (0.0, 1.0))
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-456-all-sites-vln.pdf"))

        h = Plots.histogram(node_456_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, grid = :x, xticks = (0:0.5:1))
        Plots.savefig(h,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-456-all-sites-hist.pdf"))

        v = get_split_violin_plot(
                node_12_3_probs,
                vo_node_12_3_probs,
                xlabels = [ "" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (200, 200))
        Plots.ylims!(v, (-0.02, 1.02))
        #= Plots.ylabel!(v, "Posterior probability") =#
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-12-3-vln.pdf"))

        v = get_violin_plot(
                prep_for_violin(rng, node_12_3_probs, 0.02),
                xlabels = [ "" ],
                fill_colors = gen_col,
                marker_colors = gen_col,
                fill_alphas = gen_fill_alpha,
                marker_alphas = gen_marker_alpha,
                include_dots = false)
        Plots.plot!(v, size = (120, 120), xaxis = false, grid = :y)
        #= Plots.ylims!(v, (-0.02, 1.02)) =#
        Plots.ylims!(v, (0.0, 1.0))
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-12-3-all-sites-vln.pdf"))

        h = Plots.histogram(node_12_3_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, grid = :x, xticks = (0:0.5:1))
        Plots.savefig(h,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-12-3-all-sites-hist.pdf"))

        v = get_split_violin_plot(
                height_12_789_probs,
                vo_height_12_789_probs,
                xlabels = [ "" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (200, 200))
        Plots.ylims!(v, (-0.02, 1.02))
        #= Plots.ylabel!(v, "Posterior probability") =#
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-height-probs-12-789-vln.pdf"))

        v = get_violin_plot(
                prep_for_violin(rng, height_12_789_probs, violin_probability_jiggle),
                xlabels = [ "" ],
                fill_colors = gen_col,
                marker_colors = gen_col,
                fill_alphas = gen_fill_alpha,
                marker_alphas = gen_marker_alpha,
                include_dots = false)
        Plots.plot!(v, size = (120, 120), xaxis = false, grid = :y)
        #= Plots.ylims!(v, (-0.02, 1.02)) =#
        Plots.ylims!(v, (0.0, 1.0))
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-height-probs-12-789-all-sites-vln.pdf"))

        h = Plots.histogram(height_12_789_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, grid = :x, xticks = (0:0.5:1))
        Plots.savefig(h,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-height-probs-12-789-all-sites-hist.pdf"))

        v = get_split_violin_plot(
                height_123_456_probs,
                vo_height_123_456_probs,
                xlabels = [ "" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (200, 200))
        Plots.ylims!(v, (-0.02, 1.02))
        #= Plots.ylabel!(v, "Posterior probability") =#
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-height-probs-123-456-vln.pdf"))

        v = get_violin_plot(
                prep_for_violin(rng, height_123_456_probs, violin_probability_jiggle),
                xlabels = [ "" ],
                fill_colors = gen_col,
                marker_colors = gen_col,
                fill_alphas = gen_fill_alpha,
                marker_alphas = gen_marker_alpha,
                include_dots = false)
        Plots.plot!(v, size = (120, 120), xaxis = false, grid = :y)
        #= Plots.ylims!(v, (-0.02, 1.02)) =#
        Plots.ylims!(v, (0.0, 1.0))
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-height-probs-123-456-all-sites-vln.pdf"))

        h = Plots.histogram(height_123_456_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, grid = :x, xticks = (0:0.5:1))
        Plots.savefig(h,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-height-probs-123-456-all-sites-hist.pdf"))

        v = get_split_violin_plot(
                split_12_probs,
                vo_split_12_probs,
                xlabels = [ "" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (200, 200))
        Plots.ylims!(v, (-0.02, 1.02))
        #= Plots.ylabel!(v, "Posterior probability") =#
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-split-probs-12-vln.pdf"))

        v = get_violin_plot(
                prep_for_violin(rng, split_12_probs, violin_probability_jiggle),
                xlabels = [ "" ],
                fill_colors = gen_col,
                marker_colors = gen_col,
                fill_alphas = gen_fill_alpha,
                marker_alphas = gen_marker_alpha,
                include_dots = false)
        Plots.plot!(v, size = (120, 120), xaxis = false, grid = :y)
        #= Plots.ylims!(v, (-0.02, 1.02)) =#
        Plots.ylims!(v, (0.0, 1.0))
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-split-probs-12-all-sites-vln.pdf"))

        h = Plots.histogram(split_12_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, grid = :x, xticks = (0:0.5:1))
        Plots.savefig(h,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-split-probs-12-all-sites-hist.pdf"))

        v = get_split_violin_plot(
                true_topo_probs,
                vo_true_topo_probs,
                xlabels = [ "" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (200, 200))
        Plots.ylims!(v, (-0.02, 1.02))
        #= Plots.ylabel!(v, "Posterior probability") =#
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-topo-probs-vln.pdf"))

        h = Plots.histogram(true_topo_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, grid = :x, xticks = (0:0.5:1))
        Plots.savefig(h,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-topo-probs-all-sites-hist.pdf"))

        v = get_split_violin_plot(
                hcat(root_node_probs,
                     node_789_probs,
                     node_456_probs,
                     node_12_3_probs,
                     height_12_789_probs,
                     height_123_456_probs,
                     split_12_probs,
                     true_topo_probs),
                hcat(vo_root_node_probs,
                     vo_node_789_probs,
                     vo_node_456_probs,
                     vo_node_12_3_probs,
                     vo_height_12_789_probs,
                     vo_height_123_456_probs,
                     vo_split_12_probs,
                     vo_true_topo_probs),
                xlabels = [ "Root node" "Node 789" "Node 456" "Node 12-3" "Height 12-789" "Height 123-456" "Split 12" "True topology" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (1200, 300))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-true-tree-probs.pdf"))

        # vln_root_node_probs = StatsPlots.violin(
        #         ["All sites" "SNPs"],
        #         hcat(prep_for_violin(rng, root_node_probs, violin_probability_jiggle),
        #              prep_for_violin(rng, vo_root_node_probs, violin_probability_jiggle)),
        #         ylims = (0.0, 1.0),
        #         legend = false,
        #         title = "Root node",
        #         #= size = (600, 400), =#
        # )
        # vln_node_456_probs = StatsPlots.violin(
        #         ["All sites" "SNPs"],
        #         hcat(prep_for_violin(rng, node_456_probs, violin_probability_jiggle),
        #              prep_for_violin(rng, vo_node_456_probs, violin_probability_jiggle)),
        #         ylims = (0.0, 1.0),
        #         legend = false,
        #         title = "Node 456",
        #         #= size = (600, 400), =#
        # )
        # vln_node_789_probs = StatsPlots.violin(
        #         ["All sites" "SNPs"],
        #         hcat(prep_for_violin(rng, node_789_probs, violin_probability_jiggle),
        #              prep_for_violin(rng, vo_node_789_probs, violin_probability_jiggle)),
        #         ylims = (0.0, 1.0),
        #         legend = false,
        #         title = "Node 789",
        #         #= size = (600, 400), =#
        # )
        # vln_node_12_3_probs = StatsPlots.violin(
        #         ["All sites" "SNPs"],
        #         hcat(prep_for_violin(rng, node_12_3_probs, violin_probability_jiggle),
        #              prep_for_violin(rng, vo_node_12_3_probs, violin_probability_jiggle)),
        #         ylims = (0.0, 1.0),
        #         legend = false,
        #         title = "Node 12-3",
        #         #= size = (600, 400), =#
        # )
        # vln_height_12_789_probs = StatsPlots.violin(
        #         ["All sites" "SNPs"],
        #         hcat(prep_for_violin(rng, height_12_789_probs, violin_probability_jiggle),
        #              prep_for_violin(rng, vo_height_12_789_probs, violin_probability_jiggle)),
        #         ylims = (0.0, 1.0),
        #         legend = false,
        #         title = "Height 12-789",
        #         #= size = (600, 400), =#
        # )
        # vln_height_123_456_probs = StatsPlots.violin(["All sites" "SNPs"],
        #         hcat(prep_for_violin(rng, height_123_456_probs, violin_probability_jiggle),
        #              prep_for_violin(rng, vo_height_123_456_probs, violin_probability_jiggle)),
        #         ylims = (0.0, 1.0),
        #         legend = false,
        #         title = "Height 123-456",
        #         #= size = (600, 400), =#
        # )
        # vln_split_12_probs = StatsPlots.violin(
        #         ["All sites" "SNPs"],
        #         hcat(prep_for_violin(rng, split_12_probs, violin_probability_jiggle),
        #              prep_for_violin(rng, vo_split_12_probs, violin_probability_jiggle)),
        #         ylims = (0.0, 1.0),
        #         legend = false,
        #         title = "Split 12",
        #         #= size = (600, 400), =#
        # )
        # vln_true_topo_probs = StatsPlots.violin(
        #         ["All sites" "SNPs"],
        #         hcat(prep_for_violin(rng, true_topo_probs, violin_probability_jiggle),
        #              prep_for_violin(rng, vo_true_topo_probs, violin_probability_jiggle)),
        #         ylims = (0.0, 1.0),
        #         legend = false,
        #         title = "True topology",
        #         #= size = (600, 400), =#
        # )
        # vln_grid = Plots.plot(
        #         vln_root_node_probs,
        #         vln_node_789_probs,
        #         vln_node_456_probs,
        #         vln_node_12_3_probs,
        #         vln_height_12_789_probs,
        #         vln_height_123_456_probs,
        #         vln_split_12_probs,
        #         vln_true_topo_probs,
        #         layout = (2, 4),
        #         legend = false,
        #         size = (800, 400),
        #         link = :all, # :none, :x, :y, :both, :all
        # )
        # Plots.ylabel!(vln_grid, "Posterior probability")
        # share_xy_axes!(vln_grid)

        # write(stdout, "Writing to fixed-gen-true-tree-probs.pdf\n")
        # Plots.savefig(vln_grid, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-true-tree-probs.pdf"))


        gen_max_456_subsplit_prob = get_floats(results, true, true, true, false, :max_456_subsplit_prob, locus_size)
        vo_gen_max_456_subsplit_prob = get_floats(results, true, true, true, true, :max_456_subsplit_prob, locus_size)
        bif_max_456_subsplit_prob = get_floats(results, true, true, false, false, :max_456_subsplit_prob, locus_size)
        vo_bif_max_456_subsplit_prob = get_floats(results, true, true, false, true, :max_456_subsplit_prob, locus_size)

        v = get_split_violin_plot(
                hcat(gen_max_456_subsplit_prob, vo_gen_max_456_subsplit_prob),
                hcat(bif_max_456_subsplit_prob, vo_bif_max_456_subsplit_prob),
                xlabels = ["All sites" "Variable sites"],
                left_fill_colors = [gen_col vo_gen_col],
                left_marker_colors = [gen_col vo_gen_col],
                left_fill_alphas = [gen_fill_alpha vo_gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha vo_gen_marker_alpha],
                left_labels = "Generalized",
                right_fill_colors = [bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha],
                right_labels = "Bifurcating",
                legend = false,
                dot_legend = false)
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-max-456-subsplit-probs.pdf"))


        gen_max_789_subsplit_prob = get_floats(results, true, true, true, false, :max_789_subsplit_prob, locus_size)
        vo_gen_max_789_subsplit_prob = get_floats(results, true, true, true, true, :max_789_subsplit_prob, locus_size)
        bif_max_789_subsplit_prob = get_floats(results, true, true, false, false, :max_789_subsplit_prob, locus_size)
        vo_bif_max_789_subsplit_prob = get_floats(results, true, true, false, true, :max_789_subsplit_prob, locus_size)

        v = get_split_violin_plot(
                hcat(gen_max_789_subsplit_prob, vo_gen_max_789_subsplit_prob),
                hcat(bif_max_789_subsplit_prob, vo_bif_max_789_subsplit_prob),
                xlabels = ["All sites" "Variable sites"],
                left_fill_colors = [gen_col vo_gen_col],
                left_marker_colors = [gen_col vo_gen_col],
                left_fill_alphas = [gen_fill_alpha vo_gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha vo_gen_marker_alpha],
                left_labels = "Generalized",
                right_fill_colors = [bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha],
                right_labels = "Bifurcating",
                legend = false,
                dot_legend = false)
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-max-789-subsplit-probs.pdf"))


        bif_gen_root_node_probs = get_floats(results, true, false, true, false, :root_node_true_prob, locus_size)
        vo_bif_gen_root_node_probs = get_floats(results, true, false, true, true, :root_node_true_prob, locus_size)
        bif_bif_root_node_probs = get_floats(results, true, false, false, false, :root_node_true_prob, locus_size)
        vo_bif_bif_root_node_probs = get_floats(results, true, false, false, true, :root_node_true_prob, locus_size)

        bif_gen_true_topo_probs = get_floats(results, true, false, true, false, :topo_true_prob, locus_size)
        vo_bif_gen_true_topo_probs = get_floats(results, true, false, true, true, :topo_true_prob, locus_size)
        bif_bif_true_topo_probs = get_floats(results, true, false, false, false, :topo_true_prob, locus_size)
        vo_bif_bif_true_topo_probs = get_floats(results, true, false, false, true, :topo_true_prob, locus_size)

        bif_gen_true_split_prob_mean = get_floats(results, true, false, true, false, :true_split_prob_mean, locus_size)
        vo_bif_gen_true_split_prob_mean = get_floats(results, true, false, true, true, :true_split_prob_mean, locus_size)
        bif_bif_true_split_prob_mean = get_floats(results, true, false, false, false, :true_split_prob_mean, locus_size)
        vo_bif_bif_true_split_prob_mean = get_floats(results, true, false, false, true, :true_split_prob_mean, locus_size)

        bif_gen_ht_12_78_prob = get_floats(results, true, false, true, false, :height_12_78_prob, locus_size)
        vo_bif_gen_ht_12_78_prob = get_floats(results, true, false, true, true, :height_12_78_prob, locus_size)

        bif_gen_ht_45_789_prob = get_floats(results, true, false, true, false, :height_45_789_prob, locus_size)
        vo_bif_gen_ht_45_789_prob = get_floats(results, true, false, true, true, :height_45_789_prob, locus_size)

        bif_gen_ht_456_789_prob = get_floats(results, true, false, true, false, :height_456_789_prob, locus_size)
        vo_bif_gen_ht_456_789_prob = get_floats(results, true, false, true, true, :height_456_789_prob, locus_size)

        bif_gen_nd_1_2_3_prob = get_floats(results, true, false, true, false, :node_1_2_3_prob, locus_size)
        vo_bif_gen_nd_1_2_3_prob = get_floats(results, true, false, true, true, :node_1_2_3_prob, locus_size)

        bif_gen_nd_4_5_6_prob = get_floats(results, true, false, true, false, :node_4_5_6_prob, locus_size)
        vo_bif_gen_nd_4_5_6_prob = get_floats(results, true, false, true, true, :node_4_5_6_prob, locus_size)

        bif_gen_nd_123_456_789_prob = get_floats(results, true, false, true, false, :node_123_456_789_prob, locus_size)
        vo_bif_gen_nd_123_456_789_prob = get_floats(results, true, false, true, true, :node_123_456_789_prob, locus_size)

        v = get_split_violin_plot(
                hcat(bif_gen_root_node_probs,
                     vo_bif_gen_root_node_probs),
                hcat(bif_bif_root_node_probs,
                     vo_bif_bif_root_node_probs),
                xlabels = [ "All sites" "Variable sites" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "Generalized" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Bifurcating" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (1200, 300))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        write(stdout, "Writing to fixed-bif-root-node-probs.pdf\n")
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-root-node-probs.pdf"))

        v = get_split_violin_plot(
                hcat(bif_gen_true_topo_probs,
                     vo_bif_gen_true_topo_probs),
                hcat(bif_bif_true_topo_probs,
                     vo_bif_bif_true_topo_probs),
                xlabels = [ "All sites" "Variable sites" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "Generalized" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Bifurcating" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (1200, 300))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        write(stdout, "Writing to fixed-bif-true-topo-probs.pdf\n")
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-true-topo-probs.pdf"))

        v = get_split_violin_plot(
                hcat(bif_gen_true_split_prob_mean,
                     vo_bif_gen_true_split_prob_mean),
                hcat(bif_bif_true_split_prob_mean,
                     vo_bif_bif_true_split_prob_mean),
                xlabels = [ "All sites" "Variable sites" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "Generalized" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Bifurcating" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (1200, 300))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        write(stdout, "Writing to fixed-bif-true-split-prob-means.pdf\n")
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-true-split-prob-means.pdf"))

        v = get_split_violin_plot(
                hcat(bif_gen_ht_12_78_prob,
                     bif_gen_ht_45_789_prob,
                     bif_gen_ht_456_789_prob,
                     bif_gen_nd_1_2_3_prob,
                     bif_gen_nd_4_5_6_prob,
                     bif_gen_nd_123_456_789_prob),
                hcat(vo_bif_gen_ht_12_78_prob,
                     vo_bif_gen_ht_45_789_prob,
                     vo_bif_gen_ht_456_789_prob,
                     vo_bif_gen_nd_1_2_3_prob,
                     vo_bif_gen_nd_4_5_6_prob,
                     vo_bif_gen_nd_123_456_789_prob),
                xlabels = [ "Height 12-78" "Height 45-789" "Height 456-789" "Poly 123" "Poly 456" "Root poly" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (1200, 300))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        write(stdout, "Writing to fixed-bif-gen-wrong-probs.pdf\n")
        Plots.savefig(v, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-gen-wrong-probs.pdf"))


        fixed_bif_bif_dist_mean = get_floats(results,
                true, false, false, false,
                :euclidean_distance_mean, locus_size)
        fixed_bif_bif_dist_lower = get_floats(results,
                true, false, false, false,
                :euclidean_distance_eti_95_lower, locus_size)
        fixed_bif_bif_dist_upper = get_floats(results,
                true, false, false, false,
                :euclidean_distance_eti_95_upper, locus_size)
        vo_fixed_bif_bif_dist_mean = get_floats(results,
                true, false, false, true,
                :euclidean_distance_mean, locus_size)
        vo_fixed_bif_bif_dist_lower = get_floats(results,
                true, false, false, true,
                :euclidean_distance_eti_95_lower, locus_size)
        vo_fixed_bif_bif_dist_upper = get_floats(results,
                true, false, false, true,
                :euclidean_distance_eti_95_upper, locus_size)

        fixed_bif_gen_dist_mean = get_floats(results,
                true, false, true, false,
                :euclidean_distance_mean, locus_size)
        fixed_bif_gen_dist_lower = get_floats(results,
                true, false, true, false,
                :euclidean_distance_eti_95_lower, locus_size)
        fixed_bif_gen_dist_upper = get_floats(results,
                true, false, true, false,
                :euclidean_distance_eti_95_upper, locus_size)
        vo_fixed_bif_gen_dist_mean = get_floats(results,
                true, false, true, true,
                :euclidean_distance_mean, locus_size)
        vo_fixed_bif_gen_dist_lower = get_floats(results,
                true, false, true, true,
                :euclidean_distance_eti_95_lower, locus_size)
        vo_fixed_bif_gen_dist_upper = get_floats(results,
                true, false, true, true,
                :euclidean_distance_eti_95_upper, locus_size)

        fixed_gen_bif_dist_mean = get_floats(results,
                true, true, false, false,
                :euclidean_distance_mean, locus_size)
        fixed_gen_bif_dist_lower = get_floats(results,
                true, true, false, false,
                :euclidean_distance_eti_95_lower, locus_size)
        fixed_gen_bif_dist_upper = get_floats(results,
                true, true, false, false,
                :euclidean_distance_eti_95_upper, locus_size)
        vo_fixed_gen_bif_dist_mean = get_floats(results,
                true, true, false, true,
                :euclidean_distance_mean, locus_size)
        vo_fixed_gen_bif_dist_lower = get_floats(results,
                true, true, false, true,
                :euclidean_distance_eti_95_lower, locus_size)
        vo_fixed_gen_bif_dist_upper = get_floats(results,
                true, true, false, true,
                :euclidean_distance_eti_95_upper, locus_size)

        fixed_gen_gen_dist_mean = get_floats(results,
                true, true, true, false,
                :euclidean_distance_mean, locus_size)
        fixed_gen_gen_dist_lower = get_floats(results,
                true, true, true, false,
                :euclidean_distance_eti_95_lower, locus_size)
        fixed_gen_gen_dist_upper = get_floats(results,
                true, true, true, false,
                :euclidean_distance_eti_95_upper, locus_size)
        vo_fixed_gen_gen_dist_mean = get_floats(results,
                true, true, true, true,
                :euclidean_distance_mean, locus_size)
        vo_fixed_gen_gen_dist_lower = get_floats(results,
                true, true, true, true,
                :euclidean_distance_eti_95_lower, locus_size)
        vo_fixed_gen_gen_dist_upper = get_floats(results,
                true, true, true, true,
                :euclidean_distance_eti_95_upper, locus_size)


        unfixed_bif_bif_dist_mean = get_floats(results,
                false, false, false, false,
                :euclidean_distance_mean, locus_size)
        unfixed_bif_bif_dist_lower = get_floats(results,
                false, false, false, false,
                :euclidean_distance_eti_95_lower, locus_size)
        unfixed_bif_bif_dist_upper = get_floats(results,
                false, false, false, false,
                :euclidean_distance_eti_95_upper, locus_size)
        vo_unfixed_bif_bif_dist_mean = get_floats(results,
                false, false, false, true,
                :euclidean_distance_mean, locus_size)
        vo_unfixed_bif_bif_dist_lower = get_floats(results,
                false, false, false, true,
                :euclidean_distance_eti_95_lower, locus_size)
        vo_unfixed_bif_bif_dist_upper = get_floats(results,
                false, false, false, true,
                :euclidean_distance_eti_95_upper, locus_size)

        unfixed_bif_gen_dist_mean = get_floats(results,
                false, false, true, false,
                :euclidean_distance_mean, locus_size)
        unfixed_bif_gen_dist_lower = get_floats(results,
                false, false, true, false,
                :euclidean_distance_eti_95_lower, locus_size)
        unfixed_bif_gen_dist_upper = get_floats(results,
                false, false, true, false,
                :euclidean_distance_eti_95_upper, locus_size)
        vo_unfixed_bif_gen_dist_mean = get_floats(results,
                false, false, true, true,
                :euclidean_distance_mean, locus_size)
        vo_unfixed_bif_gen_dist_lower = get_floats(results,
                false, false, true, true,
                :euclidean_distance_eti_95_lower, locus_size)
        vo_unfixed_bif_gen_dist_upper = get_floats(results,
                false, false, true, true,
                :euclidean_distance_eti_95_upper, locus_size)

        unfixed_gen_bif_dist_mean = get_floats(results,
                false, true, false, false,
                :euclidean_distance_mean, locus_size)
        unfixed_gen_bif_dist_lower = get_floats(results,
                false, true, false, false,
                :euclidean_distance_eti_95_lower, locus_size)
        unfixed_gen_bif_dist_upper = get_floats(results,
                false, true, false, false,
                :euclidean_distance_eti_95_upper, locus_size)
        vo_unfixed_gen_bif_dist_mean = get_floats(results,
                false, true, false, true,
                :euclidean_distance_mean, locus_size)
        vo_unfixed_gen_bif_dist_lower = get_floats(results,
                false, true, false, true,
                :euclidean_distance_eti_95_lower, locus_size)
        vo_unfixed_gen_bif_dist_upper = get_floats(results,
                false, true, false, true,
                :euclidean_distance_eti_95_upper, locus_size)

        unfixed_gen_gen_dist_mean = get_floats(results,
                false, true, true, false,
                :euclidean_distance_mean, locus_size)
        unfixed_gen_gen_dist_lower = get_floats(results,
                false, true, true, false,
                :euclidean_distance_eti_95_lower, locus_size)
        unfixed_gen_gen_dist_upper = get_floats(results,
                false, true, true, false,
                :euclidean_distance_eti_95_upper, locus_size)
        vo_unfixed_gen_gen_dist_mean = get_floats(results,
                false, true, true, true,
                :euclidean_distance_mean, locus_size)
        vo_unfixed_gen_gen_dist_lower = get_floats(results,
                false, true, true, true,
                :euclidean_distance_eti_95_lower, locus_size)
        vo_unfixed_gen_gen_dist_upper = get_floats(results,
                false, true, true, true,
                :euclidean_distance_eti_95_upper, locus_size)

        write(stdout, "\n")
        wsr_test = HypothesisTests.SignedRankTest(
                fixed_gen_gen_dist_mean,
                fixed_gen_bif_dist_mean)
        pval = HypothesisTests.pvalue(wsr_test, tail = :both)
        write(stdout, "WSR test, fixed gen, all sites:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(fixed_gen_gen_dist_mean - fixed_gen_bif_dist_mean))\n")
        write(stdout, "\n")
        wsr_test = HypothesisTests.SignedRankTest(
                vo_fixed_gen_gen_dist_mean,
                vo_fixed_gen_bif_dist_mean)
        pval = HypothesisTests.pvalue(wsr_test, tail = :both)
        write(stdout, "WSR test, fixed gen, var only:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(vo_fixed_gen_gen_dist_mean - vo_fixed_gen_bif_dist_mean))\n")
        write(stdout, "\n")
        wsr_test = HypothesisTests.SignedRankTest(
                fixed_bif_gen_dist_mean,
                fixed_bif_bif_dist_mean)
        pval = HypothesisTests.pvalue(wsr_test, tail = :both)
        write(stdout, "WSR test, fixed bif, all sites:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(fixed_bif_gen_dist_mean - fixed_bif_bif_dist_mean))\n")
        write(stdout, "\n")
        wsr_test = HypothesisTests.SignedRankTest(
                vo_fixed_bif_gen_dist_mean,
                vo_fixed_bif_bif_dist_mean)
        pval = HypothesisTests.pvalue(wsr_test, tail = :both)
        write(stdout, "WSR test, fixed bif, var only:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(vo_fixed_bif_gen_dist_mean - vo_fixed_bif_bif_dist_mean))\n")
        write(stdout, "\n")

        wsr_test = HypothesisTests.SignedRankTest(
                unfixed_gen_gen_dist_mean,
                unfixed_gen_bif_dist_mean)
        pval = HypothesisTests.pvalue(wsr_test, tail = :both)
        write(stdout, "WSR test, unfixed gen, all sites:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(unfixed_gen_gen_dist_mean - unfixed_gen_bif_dist_mean))\n")
        write(stdout, "\n")
        wsr_test = HypothesisTests.SignedRankTest(
                vo_unfixed_gen_gen_dist_mean,
                vo_unfixed_gen_bif_dist_mean)
        pval = HypothesisTests.pvalue(wsr_test, tail = :both)
        write(stdout, "WSR test, unfixed gen, var only:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(vo_unfixed_gen_gen_dist_mean - vo_unfixed_gen_bif_dist_mean))\n")
        write(stdout, "\n")
        wsr_test = HypothesisTests.SignedRankTest(
                unfixed_bif_gen_dist_mean,
                unfixed_bif_bif_dist_mean)
        pval = HypothesisTests.pvalue(wsr_test, tail = :both)
        write(stdout, "WSR test, unfixed bif, all sites:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(unfixed_bif_gen_dist_mean - unfixed_bif_bif_dist_mean))\n")
        write(stdout, "\n")
        wsr_test = HypothesisTests.SignedRankTest(
                vo_unfixed_bif_gen_dist_mean,
                vo_unfixed_bif_bif_dist_mean)
        pval = HypothesisTests.pvalue(wsr_test, tail = :both)
        write(stdout, "WSR test, unfixed bif, var only:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(vo_unfixed_bif_gen_dist_mean - vo_unfixed_bif_bif_dist_mean))\n")
        write(stdout, "\n")

        p = get_groups_by_y(
                hcat(fixed_gen_gen_dist_mean,
                     fixed_gen_bif_dist_mean,
                     vo_fixed_gen_gen_dist_mean,
                     vo_fixed_gen_bif_dist_mean),
                hcat(fixed_gen_gen_dist_lower,
                     fixed_gen_bif_dist_lower,
                     vo_fixed_gen_gen_dist_lower,
                     vo_fixed_gen_bif_dist_lower),
                hcat(fixed_gen_gen_dist_upper,
                     fixed_gen_bif_dist_upper,
                     vo_fixed_gen_gen_dist_upper,
                     vo_fixed_gen_bif_dist_upper),
                [ "Generalized model\n(true model)\nAll sites" "Independent bifurcating\nmodel\nAll sites" "Generalized model\n(true model)\nOnly SNPs" "Independent bifurcating\nmodel\nOnly SNPs" ],
                [ gen_col bif_col vo_gen_col vo_bif_col ],
                [ gen_marker_alpha bif_marker_alpha vo_gen_marker_alpha vo_bif_marker_alpha ],
                y_buffer = 0.02,
                show_labels_on_x = true)
        write(stdout, "Writing to fixed-gen-euclidean-distances.pdf\n")
        Plots.savefig(p, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-euclidean-distances.pdf"))

        p = get_groups_by_y(
                hcat(fixed_bif_gen_dist_mean,
                     fixed_bif_bif_dist_mean,
                     vo_fixed_bif_gen_dist_mean,
                     vo_fixed_bif_bif_dist_mean),
                hcat(fixed_bif_gen_dist_lower,
                     fixed_bif_bif_dist_lower,
                     vo_fixed_bif_gen_dist_lower,
                     vo_fixed_bif_bif_dist_lower),
                hcat(fixed_bif_gen_dist_upper,
                     fixed_bif_bif_dist_upper,
                     vo_fixed_bif_gen_dist_upper,
                     vo_fixed_bif_bif_dist_upper),
                [ "Generalized model\nAll sites" "Independent bifurcating\nmodel (true model)\nAll sites" "Generalized model\nOnly SNPs" "Independent bifurcating\nmodel (true model)\nOnly SNPs" ],
                [ gen_col bif_col vo_gen_col vo_bif_col ],
                [ gen_marker_alpha bif_marker_alpha vo_gen_marker_alpha vo_bif_marker_alpha ],
                y_buffer = 0.02,
                show_labels_on_x = true)
        Plots.ylabel!(p, "Euclidean distance from true tree")
        write(stdout, "Writing to fixed-bif-euclidean-distances.pdf\n")
        Plots.savefig(p, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-euclidean-distances.pdf"))

        p = get_groups_by_y(
                hcat(unfixed_gen_gen_dist_mean,
                     unfixed_gen_bif_dist_mean,
                     vo_unfixed_gen_gen_dist_mean,
                     vo_unfixed_gen_bif_dist_mean),
                hcat(unfixed_gen_gen_dist_lower,
                     unfixed_gen_bif_dist_lower,
                     vo_unfixed_gen_gen_dist_lower,
                     vo_unfixed_gen_bif_dist_lower),
                hcat(unfixed_gen_gen_dist_upper,
                     unfixed_gen_bif_dist_upper,
                     vo_unfixed_gen_gen_dist_upper,
                     vo_unfixed_gen_bif_dist_upper),
                [ "Generalized model\n(true model)\nAll sites" "Independent bifurcating\nmodel\nAll sites" "Generalized model\n(true model)\nOnly SNPs" "Independent bifurcating\nmodel\nOnly SNPs" ],
                [ gen_col bif_col vo_gen_col vo_bif_col ],
                [ gen_marker_alpha bif_marker_alpha vo_gen_marker_alpha vo_bif_marker_alpha ],
                y_buffer = 0.02,
                show_labels_on_x = true)
        write(stdout, "Writing to unfixed-gen-euclidean-distances.pdf\n")
        Plots.savefig(p, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-euclidean-distances.pdf"))

        p = get_groups_by_y(
                hcat(unfixed_gen_gen_dist_mean,
                     unfixed_gen_bif_dist_mean),
                hcat(unfixed_gen_gen_dist_lower,
                     unfixed_gen_bif_dist_lower),
                hcat(unfixed_gen_gen_dist_upper,
                     unfixed_gen_bif_dist_upper),
                [ "Generalized model\n(true model)\nAll sites" "Independent bifurcating\nmodel\nAll sites" ],
                [ gen_col bif_col ],
                [ gen_marker_alpha bif_marker_alpha ],
                y_buffer = 0.02,
                show_labels_on_x = true)
        write(stdout, "Writing to unfixed-gen-all-sites-euclidean-distances.pdf\n")
        Plots.savefig(p, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-all-sites-euclidean-distances.pdf"))

        p = get_groups_by_y(
                hcat(vo_unfixed_gen_gen_dist_mean,
                     vo_unfixed_gen_bif_dist_mean),
                hcat(vo_unfixed_gen_gen_dist_lower,
                     vo_unfixed_gen_bif_dist_lower),
                hcat(vo_unfixed_gen_gen_dist_upper,
                     vo_unfixed_gen_bif_dist_upper),
                [ "Generalized model\n(true model)\nOnly SNPs" "Independent bifurcating\nmodel\nOnly SNPs" ],
                [ vo_gen_col vo_bif_col ],
                [ vo_gen_marker_alpha vo_bif_marker_alpha ],
                y_buffer = 0.02,
                show_labels_on_x = true)
        write(stdout, "Writing to unfixed-gen-var-only-euclidean-distances.pdf\n")
        Plots.savefig(p, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-var-only-euclidean-distances.pdf"))

        p = get_groups_by_y(
                hcat(unfixed_bif_gen_dist_mean,
                     unfixed_bif_bif_dist_mean,
                     vo_unfixed_bif_gen_dist_mean,
                     vo_unfixed_bif_bif_dist_mean),
                hcat(unfixed_bif_gen_dist_lower,
                     unfixed_bif_bif_dist_lower,
                     vo_unfixed_bif_gen_dist_lower,
                     vo_unfixed_bif_bif_dist_lower),
                hcat(unfixed_bif_gen_dist_upper,
                     unfixed_bif_bif_dist_upper,
                     vo_unfixed_bif_gen_dist_upper,
                     vo_unfixed_bif_bif_dist_upper),
                [ "Generalized model\nAll sites" "Independent bifurcating\nmodel (true model)\nAll sites" "Generalized model\nOnly SNPs" "Independent bifurcating\nmodel (true model)\nOnly SNPs" ],
                [ gen_col bif_col vo_gen_col vo_bif_col ],
                [ gen_marker_alpha bif_marker_alpha vo_gen_marker_alpha vo_bif_marker_alpha ],
                y_buffer = 0.02,
                show_labels_on_x = true)
        Plots.ylabel!(p, "Euclidean distance from true tree")
        write(stdout, "Writing to unfixed-bif-euclidean-distances.pdf\n")
        Plots.savefig(p, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-bif-euclidean-distances.pdf"))

        p = get_groups_by_y(
                hcat(unfixed_bif_gen_dist_mean,
                     unfixed_bif_bif_dist_mean),
                hcat(unfixed_bif_gen_dist_lower,
                     unfixed_bif_bif_dist_lower),
                hcat(unfixed_bif_gen_dist_upper,
                     unfixed_bif_bif_dist_upper),
                [ "Generalized model\nAll sites" "Independent bifurcating\nmodel (true model)\nAll sites" ],
                [ gen_col bif_col ],
                [ gen_marker_alpha bif_marker_alpha ],
                y_buffer = 0.02,
                show_labels_on_x = true)
        Plots.ylabel!(p, "Euclidean distance from true tree")
        write(stdout, "Writing to unfixed-bif-all-sites-euclidean-distances.pdf\n")
        Plots.savefig(p, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-bif-all-sites-euclidean-distances.pdf"))

        p = get_groups_by_y(
                hcat(vo_unfixed_bif_gen_dist_mean,
                     vo_unfixed_bif_bif_dist_mean),
                hcat(vo_unfixed_bif_gen_dist_lower,
                     vo_unfixed_bif_bif_dist_lower),
                hcat(vo_unfixed_bif_gen_dist_upper,
                     vo_unfixed_bif_bif_dist_upper),
                [ "Generalized model\nOnly SNPs" "Independent bifurcating\nmodel (true model)\nOnly SNPs" ],
                [ vo_gen_col vo_bif_col ],
                [ vo_gen_marker_alpha vo_bif_marker_alpha ],
                y_buffer = 0.02,
                show_labels_on_x = true)
        Plots.ylabel!(p, "Euclidean distance from true tree")
        write(stdout, "Writing to unfixed-bif-var-only-euclidean-distances.pdf\n")
        Plots.savefig(p, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-bif-var-only-euclidean-distances.pdf"))


        gen_gen_tree_len = get_true_v_est_plot(
                results,
                false,
                true,
                true,
                false,
                "tree_length",
                locus_size,
                false,
                0.05,
                brooks_gelman_1998_recommended_psrf,
                ess_threshold,
                "red")
        tree_len_coverage = get_floats(
                sum_results,
                false,
                true,
                true,
                false,
                :coverage_tree_length,
                locus_size)[1]
        #= sss = "p(TL " * "\\in" * " CI) = $(tree_len_coverage)" =#
        #= sss = "p(TL $(L"\in") CI) = $(tree_len_coverage)" =#
        annotate!(gen_gen_tree_len,
                relative_xy(gen_gen_tree_len, 0.05, 0.99)...,
                #= text(L"$p(\textrm{TL} \in \textrm{CI}) = %$tree_len_coverage$", =#
                #= text(sss, =#
                #= text("\$\\textrm{\\sffamily\\it p}(\\textrm{\\sffamily TL} \\in \\textrm{\\sffamily CI}) = $(tree_len_coverage)\$", =#
                text(L"$\textrm{\sffamily\it p}(\textrm{\sffamily TL} \in \textrm{\sffamily CI}) = %$(tree_len_coverage)$",
                        :left,
                        :top,
                        8),
                annotation_clip = false)
        vo_gen_gen_tree_len = get_true_v_est_plot(
                results,
                false,
                true,
                true,
                true,
                "tree_length",
                locus_size)
        gen_bif_tree_len = get_true_v_est_plot(
                results,
                false,
                true,
                false,
                false,
                "tree_length",
                locus_size)
        vo_gen_bif_tree_len = get_true_v_est_plot(
                results,
                false,
                true,
                false,
                true,
                "tree_length",
                locus_size)

        bif_gen_tree_len = get_true_v_est_plot(
                results,
                false,
                false,
                true,
                false,
                "tree_length",
                locus_size)
        vo_bif_gen_tree_len = get_true_v_est_plot(
                results,
                false,
                false,
                true,
                true,
                "tree_length",
                locus_size)
        bif_bif_tree_len = get_true_v_est_plot(
                results,
                false,
                false,
                false,
                false,
                "tree_length",
                locus_size)
        vo_bif_bif_tree_len = get_true_v_est_plot(
                results,
                false,
                false,
                false,
                true,
                "tree_length",
                locus_size)

        #= Plots.savefig(gen_gen_tree_len, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-gen-tree-length-allsites.pdf")) =#
        #= Plots.savefig(bif_bif_tree_len, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-bif-bif-tree-length-allsites.pdf")) =#

        #= write(stdout, "\n\n$(plt_tree_len.subplots[1][:xaxis][:lims])\n\n") =#
        #= write(stdout, "\n\n$(plt_tree_len.subplots[1][:yaxis][:lims])\n\n") =#

        #= Plots.savefig(plt_tree_len, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)tree-length-unfixed-gen-gen-all.pdf")) =#

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
        write(stdout, "Writing to unfixed-tree-length.pdf\n")
        Plots.savefig(tree_len_grid, joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-tree-length.pdf"))


        # Plot ASDSF
        fixed_gen_gen_asdsf = get_floats(results, true, true, true, false, :sdsf_mean, locus_size)
        vo_fixed_gen_gen_asdsf = get_floats(results, true, true, true, true, :sdsf_mean, locus_size)
        fixed_gen_bif_asdsf = get_floats(results, true, true, false, false, :sdsf_mean, locus_size)
        vo_fixed_gen_bif_asdsf = get_floats(results, true, true, false, true, :sdsf_mean, locus_size)
        fixed_bif_gen_asdsf = get_floats(results, true, false, true, false, :sdsf_mean, locus_size)
        vo_fixed_bif_gen_asdsf = get_floats(results, true, false, true, true, :sdsf_mean, locus_size)
        fixed_bif_bif_asdsf = get_floats(results, true, false, false, false, :sdsf_mean, locus_size)
        vo_fixed_bif_bif_asdsf = get_floats(results, true, false, false, true, :sdsf_mean, locus_size)

        vln_fixed_gen_asdsf = get_split_violin_plot(
                hcat(fixed_gen_gen_asdsf, vo_fixed_gen_gen_asdsf),
                hcat(fixed_gen_bif_asdsf, vo_fixed_gen_bif_asdsf),
                xlabels = ["All sites" "Variable sites"],
                left_fill_colors = [gen_col vo_gen_col],
                left_marker_colors = [gen_col vo_gen_col],
                left_fill_alphas = [gen_fill_alpha vo_gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha vo_gen_marker_alpha],
                left_labels = "Generalized",
                right_fill_colors = [bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha],
                right_labels = "Bifurcating",
                legend = false,
                dot_legend = false)
        Plots.savefig(vln_fixed_gen_asdsf,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed_gen_asdsf.pdf"))

        vln_fixed_bif_asdsf = get_split_violin_plot(
                hcat(fixed_bif_gen_asdsf, vo_fixed_bif_gen_asdsf),
                hcat(fixed_bif_bif_asdsf, vo_fixed_bif_bif_asdsf),
                xlabels = ["All sites" "Variable sites"],
                left_fill_colors = [gen_col vo_gen_col],
                left_marker_colors = [gen_col vo_gen_col],
                left_fill_alphas = [gen_fill_alpha vo_gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha vo_gen_marker_alpha],
                left_labels = "Generalized",
                right_fill_colors = [bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha],
                right_labels = "Bifurcating",
                legend = false,
                dot_legend = false)
        Plots.savefig(vln_fixed_bif_asdsf,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed_bif_asdsf.pdf"))

        Plots.title!(vln_fixed_gen_asdsf, "Generalized")
        Plots.title!(vln_fixed_bif_asdsf, "Bifurcating")

        vln_fixed_asdsf_grid = Plots.plot(
                vln_fixed_gen_asdsf,
                vln_fixed_bif_asdsf,
                layout = (1, 2),
                #= legend = false, =#
                size = (800, 300),
                link = :all,
        )
        Plots.ylabel!(vln_fixed_asdsf_grid, "ASDSF")
        share_xy_axes!(vln_fixed_asdsf_grid)
        write(stdout, "Writing to fixed-asdsf.pdf\n")
        Plots.savefig(vln_fixed_asdsf_grid,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-asdsf.pdf"))


        unfixed_gen_gen_asdsf = get_floats(results, false, true, true, false, :sdsf_mean, locus_size)
        vo_unfixed_gen_gen_asdsf = get_floats(results, false, true, true, true, :sdsf_mean, locus_size)
        unfixed_gen_bif_asdsf = get_floats(results, false, true, false, false, :sdsf_mean, locus_size)
        vo_unfixed_gen_bif_asdsf = get_floats(results, false, true, false, true, :sdsf_mean, locus_size)
        unfixed_bif_gen_asdsf = get_floats(results, false, false, true, false, :sdsf_mean, locus_size)
        vo_unfixed_bif_gen_asdsf = get_floats(results, false, false, true, true, :sdsf_mean, locus_size)
        unfixed_bif_bif_asdsf = get_floats(results, false, false, false, false, :sdsf_mean, locus_size)
        vo_unfixed_bif_bif_asdsf = get_floats(results, false, false, false, true, :sdsf_mean, locus_size)

        vln_unfixed_gen_asdsf = get_split_violin_plot(
                hcat(unfixed_gen_gen_asdsf, vo_unfixed_gen_gen_asdsf),
                hcat(unfixed_gen_bif_asdsf, vo_unfixed_gen_bif_asdsf),
                xlabels = ["All sites" "Variable sites"],
                left_fill_colors = [gen_col vo_gen_col],
                left_marker_colors = [gen_col vo_gen_col],
                left_fill_alphas = [gen_fill_alpha vo_gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha vo_gen_marker_alpha],
                left_labels = "Generalized",
                right_fill_colors = [bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha],
                right_labels = "Bifurcating",
                legend = false,
                dot_legend = false)
        Plots.savefig(vln_unfixed_gen_asdsf,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed_gen_asdsf.pdf"))

        vln_unfixed_bif_asdsf = get_split_violin_plot(
                hcat(unfixed_bif_gen_asdsf, vo_unfixed_bif_gen_asdsf),
                hcat(unfixed_bif_bif_asdsf, vo_unfixed_bif_bif_asdsf),
                xlabels = ["All sites" "Variable sites"],
                left_fill_colors = [gen_col vo_gen_col],
                left_marker_colors = [gen_col vo_gen_col],
                left_fill_alphas = [gen_fill_alpha vo_gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha vo_gen_marker_alpha],
                left_labels = "Generalized",
                right_fill_colors = [bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha],
                right_labels = "Bifurcating",
                legend = false,
                dot_legend = false)
        Plots.savefig(vln_unfixed_bif_asdsf,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed_bif_asdsf.pdf"))

        Plots.title!(vln_unfixed_gen_asdsf, "Generalized")
        Plots.title!(vln_unfixed_bif_asdsf, "Bifurcating")

        vln_unfixed_asdsf_grid = Plots.plot(
                vln_unfixed_gen_asdsf,
                vln_unfixed_bif_asdsf,
                layout = (1, 2),
                #= legend = false, =#
                size = (800, 300),
                link = :all,
        )
        Plots.ylabel!(vln_unfixed_asdsf_grid, "ASDSF")
        share_xy_axes!(vln_unfixed_asdsf_grid)
        write(stdout, "Writing to unfixed-asdsf.pdf\n")
        Plots.savefig(vln_unfixed_asdsf_grid,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-asdsf.pdf"))


        # Plot topology-related probs

        true_topo_probs = get_floats(results, false, true, true, false, :topo_true_prob, locus_size)
        vo_true_topo_probs = get_floats(results, false, true, true, true, :topo_true_prob, locus_size)
        root_node_probs = get_floats(results, false, true, true, false, :root_node_true_prob, locus_size)
        vo_root_node_probs = get_floats(results, false, true, true, true, :root_node_true_prob, locus_size)
        true_split_prob_means = get_floats(results, false, true, true, false, :true_split_prob_mean, locus_size)
        vo_true_split_prob_means = get_floats(results, false, true, true, true, :true_split_prob_mean, locus_size)
        true_node_prob_means = get_floats(results, false, true, true, false, :true_node_prob_mean, locus_size)
        vo_true_node_prob_means = get_floats(results, false, true, true, true, :true_node_prob_mean, locus_size)
        true_height_prob_means = get_floats(results, false, true, true, false, :true_height_prob_mean, locus_size)
        vo_true_height_prob_means = get_floats(results, false, true, true, true, :true_height_prob_mean, locus_size)

        v = get_floats(results, false, true, true, false, :true_shared_height_prob_mean, locus_size)
        true_shared_height_prob_means = v[.!any.(isnan, v)]
        v = get_floats(results, false, true, true, true, :true_shared_height_prob_mean, locus_size)
        vo_true_shared_height_prob_means = v[.!any.(isnan, v)]

        v = get_split_violin_plot(
                true_shared_height_prob_means,
                vo_true_shared_height_prob_means,
                xlabels = [ "" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (200, 200))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-gen-true-shared-height-probs.pdf"))

        v = get_split_violin_plot(
                hcat(true_topo_probs,
                     root_node_probs,
                     true_split_prob_means,
                     true_node_prob_means,
                     true_height_prob_means),
                hcat(vo_true_topo_probs,
                     vo_root_node_probs,
                     vo_true_split_prob_means,
                     vo_true_node_prob_means,
                     vo_true_height_prob_means),
                xlabels = [ "True topology" "Root node" "True split mean" "True node mean" "True height mean" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (1200, 300))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-gen-true-tree-probs.pdf"))


        true_topo_probs = get_floats(results, false, false, false, false, :topo_true_prob, locus_size)
        vo_true_topo_probs = get_floats(results, false, false, false, true, :topo_true_prob, locus_size)
        root_node_probs = get_floats(results, false, false, false, false, :root_node_true_prob, locus_size)
        vo_root_node_probs = get_floats(results, false, false, false, true, :root_node_true_prob, locus_size)
        true_split_prob_means = get_floats(results, false, false, false, false, :true_split_prob_mean, locus_size)
        vo_true_split_prob_means = get_floats(results, false, false, false, true, :true_split_prob_mean, locus_size)
        true_node_prob_means = get_floats(results, false, false, false, false, :true_node_prob_mean, locus_size)
        vo_true_node_prob_means = get_floats(results, false, false, false, true, :true_node_prob_mean, locus_size)
        true_height_prob_means = get_floats(results, false, false, false, false, :true_height_prob_mean, locus_size)
        vo_true_height_prob_means = get_floats(results, false, false, false, true, :true_height_prob_mean, locus_size)

        v = get_split_violin_plot(
                hcat(true_topo_probs,
                     root_node_probs,
                     true_split_prob_means,
                     true_node_prob_means,
                     true_height_prob_means),
                hcat(vo_true_topo_probs,
                     vo_root_node_probs,
                     vo_true_split_prob_means,
                     vo_true_node_prob_means,
                     vo_true_height_prob_means),
                xlabels = [ "True topology" "Root node" "True split mean" "True node mean" "True height mean" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (1200, 300))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-bif-bif-true-tree-probs.pdf"))


        true_topo_probs = get_floats(results, false, true, false, false, :topo_true_prob, locus_size)
        vo_true_topo_probs = get_floats(results, false, true, false, true, :topo_true_prob, locus_size)
        root_node_probs = get_floats(results, false, true, false, false, :root_node_true_prob, locus_size)
        vo_root_node_probs = get_floats(results, false, true, false, true, :root_node_true_prob, locus_size)
        true_split_prob_means = get_floats(results, false, true, false, false, :true_split_prob_mean, locus_size)
        vo_true_split_prob_means = get_floats(results, false, true, false, true, :true_split_prob_mean, locus_size)
        true_node_prob_means = get_floats(results, false, true, false, false, :true_node_prob_mean, locus_size)
        vo_true_node_prob_means = get_floats(results, false, true, false, true, :true_node_prob_mean, locus_size)
        true_height_prob_means = get_floats(results, false, true, false, false, :true_height_prob_mean, locus_size)
        vo_true_height_prob_means = get_floats(results, false, true, false, true, :true_height_prob_mean, locus_size)

        v = get_floats(results, false, true, false, false, :true_shared_height_prob_mean, locus_size)
        true_shared_height_prob_means = v[.!any.(isnan, v)]
        v = get_floats(results, false, true, false, true, :true_shared_height_prob_mean, locus_size)
        vo_true_shared_height_prob_means = v[.!any.(isnan, v)]

        v = get_split_violin_plot(
                true_shared_height_prob_means,
                vo_true_shared_height_prob_means,
                xlabels = [ "" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (200, 200))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-bif-true-shared-height-probs.pdf"))

        v = get_split_violin_plot(
                hcat(true_topo_probs,
                     root_node_probs,
                     true_split_prob_means,
                     true_node_prob_means,
                     true_height_prob_means),
                hcat(vo_true_topo_probs,
                     vo_root_node_probs,
                     vo_true_split_prob_means,
                     vo_true_node_prob_means,
                     vo_true_height_prob_means),
                xlabels = [ "True topology" "Root node" "True split mean" "True node mean" "True height mean" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (1200, 300))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-bif-true-tree-probs.pdf"))


        true_topo_probs = get_floats(results, false, false, true, false, :topo_true_prob, locus_size)
        vo_true_topo_probs = get_floats(results, false, false, true, true, :topo_true_prob, locus_size)
        root_node_probs = get_floats(results, false, false, true, false, :root_node_true_prob, locus_size)
        vo_root_node_probs = get_floats(results, false, false, true, true, :root_node_true_prob, locus_size)
        true_split_prob_means = get_floats(results, false, false, true, false, :true_split_prob_mean, locus_size)
        vo_true_split_prob_means = get_floats(results, false, false, true, true, :true_split_prob_mean, locus_size)
        true_node_prob_means = get_floats(results, false, false, true, false, :true_node_prob_mean, locus_size)
        vo_true_node_prob_means = get_floats(results, false, false, true, true, :true_node_prob_mean, locus_size)
        true_height_prob_means = get_floats(results, false, false, true, false, :true_height_prob_mean, locus_size)
        vo_true_height_prob_means = get_floats(results, false, false, true, true, :true_height_prob_mean, locus_size)

        v = get_split_violin_plot(
                hcat(true_topo_probs,
                     root_node_probs,
                     true_split_prob_means,
                     true_node_prob_means,
                     true_height_prob_means),
                hcat(vo_true_topo_probs,
                     vo_root_node_probs,
                     vo_true_split_prob_means,
                     vo_true_node_prob_means,
                     vo_true_height_prob_means),
                xlabels = [ "True topology" "Root node" "True split mean" "True node mean" "True height mean" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (1200, 300))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        Plots.savefig(v,
                joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-bif-gen-true-tree-probs.pdf"))

    end

    return 0
end

main_cli()
