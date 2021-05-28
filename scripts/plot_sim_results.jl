#! /usr/bin/env julia

using Random
using Distributions
using Statistics
using DataDeps
using DataFrames
using CSV
using Plots
using Plots.PlotMeasures
using StatsPlots
using HypothesisTests
using LaTeXStrings
using Printf
# using GR

# Set backend to GR
# gr()
pgfplotsx()
default(tex_output_standalone = true)
#= font_str = "Computer Modern Sans Serif" =#
#= font_obj = Plots.font(font_str) =#
#= default(guidefont = (18, font_str), =#
#=         xtickfont = (12, font_str), =#
#=         ytickfont = (12, font_str), =#
#=         tickfont = (12, font_str), =#
#=         legendfont = (18, font_str), =#
#=         titlefont = (20, font_str), =#
#=         fontfamily = font_str, =#
#=         tickfontfamily = font_str, =#
#=         guidefontfamily = font_str) =#
# pyplot()
# plotlyjs()

#= write(stdout, "$(plotattr(:Plot))\n\n") =#
#= write(stdout, "$(plotattr(:Axis))\n\n") =#
#= write(stdout, "$(plotattr("fontfamily"))\n\n") =#

#= default(fontfamily = "sans-serif") =#

#= write(stdout, "$(plotattr("fontfamily"))\n\n") =#

include("project_util.jl")

sans_tex_lines = [
        "\\usepackage[cm]{sfmath}",
        "\\usepackage[T1]{fontenc}",
        "\\renewcommand{\\familydefault}{\\sfdefault}",
        ]
sans_tex_insert = r"^\\usepackage\{pgfplots\}$"

dark_blue_col = RGB(2/255, 75/255, 120/255)
l_dark_blue_col = RGB(3/255, 138/255, 221/255)
comp_blue_col= RGB(2/255, 219/255, 240/255)
dark_orange_col = RGB(186/255, 88/255, 0/255)
l_dark_orange_col = RGB(255/255, 135/255, 31/255)
comp_orange_col = RGB(255/255, 181/255, 91/255)
highlight_col = "red"

# terminal.sexy
#= dblue_col = RGB(95/255, 129/255, 157/255) =#
dblue_col = RGB(90/255, 122/255, 148/255)
#= lblue_col = RGB(129/255, 162/255, 190/255) =#
lblue_col = RGB(142/255, 178/255, 209/255)
#= dpurp_col = RGB(133/255, 103/255, 143/255) =#
#= dpurp_col = RGB(117/255, 90/255, 126/255) =#
dpurp_col = RGB(113/255, 74/255, 126/255)
#= lpurp_col = RGB(178/255, 148/255, 187/255) =#
#= lpurp_col = RGB(189/255, 157/255, 199/255) =#
lpurp_col = RGB(182/255, 132/255, 199/255)
dgreen_col = RGB(94/255, 141/255, 135/255)
lgreen_col = RGB(138/255, 190/255, 183/255)

# How I got these colors from matplotlib:
# from matplotlib import cm
# v = cm.get_cmap("viridis")
# viridis015 = v(0.15, bytes = True)
viridis000 = RGB(68/255, 1/255, 84/255)
viridis005 = RGB(71/255, 18/255, 101/255)
viridis010 = RGB(72/255, 35/255, 116/255)
viridis015 = RGB(69/255, 52/255, 127/255)
viridis020 = RGB(64/255, 67/255, 135/255)
viridis025 = RGB(58/255, 82/255, 139/255)
viridis030 = RGB(52/255, 94/255, 141/255)
viridis035 = RGB(46/255, 107/255, 142/255)
viridis040 = RGB(41/255, 120/255, 142/255)
viridis045 = RGB(36/255, 132/255, 141/255)
viridis050 = RGB(32/255, 144/255, 140/255)
viridis055 = RGB(30/255, 155/255, 137/255)
viridis060 = RGB(34/255, 167/255, 132/255)
viridis065 = RGB(47/255, 179/255, 123/255)
viridis070 = RGB(68/255, 190/255, 112/255)
viridis075 = RGB(94/255, 201/255, 97/255)
viridis080 = RGB(121/255, 209/255, 81/255)
viridis085 = RGB(154/255, 216/255, 60/255)
viridis090 = RGB(189/255, 222/255, 38/255)
viridis095 = RGB(223/255, 227/255, 24/255)


#= gen_col = dark_blue_col =#
#= vo_gen_col = comp_blue_col =#
gen_col = dpurp_col
vo_gen_col = lpurp_col
#= gen_col = viridis020 =#
#= vo_gen_col = viridis035 =# 
#= vo_gen_col = dark_blue_col =#
#= bif_col = dark_orange_col =# 
#= vo_bif_col = comp_orange_col =# 
bif_col = dblue_col
vo_bif_col = lblue_col 
#= bif_col = viridis065 =#
#= vo_bif_col = viridis080 =#
#= vo_bif_col = dark_orange_col =# 

marker_alpha = 0.9
fill_alpha = 0.75

gen_shape = :circle
bif_shape = :diamond
vo_gen_shape = gen_shape
vo_bif_shape = bif_shape

gen_marker_size = 3.0
bif_marker_size = 3.5
vo_gen_marker_size = gen_marker_size
vo_bif_marker_size = bif_marker_size

gen_model = L"$M_{G}$"
bif_model = L"$M_{IB}$"

gen_marker_alpha = marker_alpha
#= vo_gen_marker_alpha = marker_alpha - 0.3 =#
vo_gen_marker_alpha = marker_alpha
gen_fill_alpha = fill_alpha
#= vo_gen_fill_alpha = fill_alpha - 0.3 =#
vo_gen_fill_alpha = fill_alpha

bif_marker_alpha = marker_alpha
vo_bif_marker_alpha = marker_alpha
#= vo_bif_marker_alpha = 0.8 =#
bif_fill_alpha = fill_alpha
vo_bif_fill_alpha = fill_alpha
#= vo_bif_fill_alpha = 0.65 =#

scatter_marker_col = RGB(96/250, 96/250, 96/250)
scatter_marker_alpha = 0.8
#= identity_line_color = RGB(192/250, 192/250, 192/250) =#
identity_line_color = scatter_marker_col
identity_line_alpha = 0.18
error_bar_alpha = 0.2
error_bar_rgba = RGBA(96/250, 96/250, 96/250, error_bar_alpha)

axis_pattern = r"\\begin\{axis\}\["
axis_replace = "\\begin{axis}[clip=false, "

struct SimResults
    df::DataFrame
    suffix::AbstractString
    function SimResults(; suffix = "")
        main_tsv_path = joinpath(ProjectUtil.RESULTS_DIR, "results$(suffix).tsv")
        if Base.Filesystem.isfile(main_tsv_path * ".gz")
            results = ProjectUtil.get_data_frame_from_gz([main_tsv_path * ".gz"])
            return new(results, suffix)
        end
        sim_conditions = ["bifurcating", "generalized"]
        analysis_conditions = ["bifurcating", "generalized"]
        if length(suffix) > 0
            sim_conditions = ["generalized"]
            analysis_conditions = ["generalized"]
        end
        results = DataFrame()
        for is_fixed in [true, false]
            for sim_condition in sim_conditions
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
                            for analysis_condition in analysis_conditions
                                var_only_bools = [true, false]
                                if endswith(locus_condition, "-linked-")
                                    var_only_bools = [false]
                                elseif endswith(locus_condition, "-unlinked-")
                                    var_only_bools = [true]
                                end
                                for is_var_only in var_only_bools
                                    result_file_name = "results$(suffix)-species-9-genomes-2-$analysis_condition-tree-random.tsv.gz"
                                    if is_var_only & (locus_condition != "locus-100-unlinked-")
                                        result_file_name = "results$(suffix)-var-only-species-9-genomes-2-$analysis_condition-tree-random.tsv.gz"
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
        return new(results, suffix)
    end
    function SimResults(data_frame::DataFrame; suffix = "")
        return new(data_frame, suffix)
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
        marker_shapes = nothing,
        marker_sizes = nothing,
        y_buffer = 0.02,
        y_axis_max = -1.0,
        error_bar_alpha = 0.2,
        show_group_means = true,
        group_mean_line_alpha = 0.8,
        group_mean_line_width = 2.0,
        show_labels_on_x = false,
        x_label_size = 8,
        x_label_offset = 0.1,
        comparisons = (),
        comparison_positions = (),
        comparisons_are_paired = false,
        p_val_size = 8)::Plots.Plot
    @assert ndims(y) == 2
    @assert size(y) == size(y_lower) 
    @assert size(y) == size(y_upper) 
    @assert size(y)[2] == size(colors)[2]
    @assert size(y)[2] == size(labels)[2]
    if isnothing(marker_shapes)
        marker_shapes = fill(:circle, 1, size(y)[2])
    end
    if isnothing(marker_sizes)
        marker_sizes = fill(4.0, 1, size(y)[2])
    end
    @assert size(y)[2] == size(marker_shapes)[2]
    @assert size(y)[2] == size(marker_sizes)[2]

    nsamples = size(y)[1]
    ngroups = size(y)[2]
    x = Array{Int64}(undef, size(y)...)
    x_val = 0
    for i in eachindex(x)
        x_val += 1
        x[i] = x_val
    end
    y_extremes = (minimum(y_lower), maximum(y_upper))
    if (y_axis_max > 0.0) & (y_axis_max > y_extremes[1])
        y_extremes = (minimum(y_lower), y_axis_max)
        @assert y_axis_max >= maximum(y)
    end
    y_buf = (y_extremes[2] - y_extremes[1]) * y_buffer
    #= y_limits = [y_extremes[1] - y_buf, y_extremes[2] + y_buf] =#
    y_limits = [0.0, y_extremes[2] + y_buf]
    x_limits = [0, x_val + 1]

    error_bar_overflows_x = []
    error_bar_overflows_y = []
    for group_idx in 1:ngroups
        overflows_x = []
        overflows_y = []
        for val_idx in eachindex(y_upper[:, group_idx])
            if (y_axis_max > 0.0) & (y_upper[val_idx, group_idx] > y_axis_max)
                y_upper[val_idx, group_idx] = y_axis_max
                push!(overflows_x, x[val_idx, group_idx])
                push!(overflows_y, y_upper[val_idx, group_idx])
            end
        end
        push!(error_bar_overflows_x, overflows_x)
        push!(error_bar_overflows_y, overflows_y)
    end

    # Get empty plot
    plt = Plots.plot(
            legend = false,
            xlims = x_limits,
            ylims = y_limits)

    # Handle WSR comparison tests and P-value labels first so that error bars
    # will be placed above the labels
    if (! isempty(comparisons))
        if isempty(comparison_positions)
            comparison_positions = Tuple(0.98 for i in comparisons)
        else
            @assert length(comparisons) == length(comparison_positions)
        end
    end
    for (comp_idx, comp) in enumerate(comparisons)
        @assert length(comp) == 2
        comp_pos = comparison_positions[comp_idx]
        y_range = y_limits[2] - y_limits[1]
        y_pval = y_limits[1] + (comp_pos * y_range)
        y_bracket_ln = y_pval - (0.03 * y_range)
        y_bracket_tick = y_bracket_ln - (0.015 * y_range)
        pval_aln = :bottom
        if comp_pos < 0.5
            y_bracket_ln = y_pval + (0.035 * y_range)
            y_bracket_tick = y_bracket_ln + (0.015 * y_range)
            pval_aln = :top
        end
        if comparisons_are_paired
            comp_test = HypothesisTests.SignedRankTest(
                    y[:,comp[1]],
                    y[:,comp[2]])
        else
            comp_test = HypothesisTests.MannWhitneyUTest(
                    y[:,comp[1]],
                    y[:,comp[2]])
        end
        pval = HypothesisTests.pvalue(comp_test, tail = :both)
        pval_str = @sprintf("%.2g", pval)
        m = match(ProjectUtil.SCI_NOTATION_PATTERN, pval_str)
        if ! isnothing(m)
            e_sign = ""
            if m[:sign] == "-"
                e_sign = "-"
            end
            pval_str = "$(m[:digits]) \\times 10^{$(e_sign)$(m[:exponent])}"
        end
        #= mid_1 = Statistics.mean(x[:,comp[1]]) =#
        #= mid_2 = Statistics.mean(x[:,comp[2]]) =#
        mid_1 = Statistics.quantile(x[:,comp[1]], 0.8)
        mid_2 = Statistics.quantile(x[:,comp[2]], 0.2)
        x_position = mid_1 + ((mid_2 - mid_1) / 2)
        annotate!(plt, x_position, y_pval,
                  text(L"p = %$(pval_str)",
                       :center,
                       pval_aln,
                       p_val_size,
                       rotation = 0),
                  annotation_clip = false)
        Plots.plot!(plt,
                [mid_1, mid_1, mid_2, mid_2],
                [y_bracket_tick, y_bracket_ln, y_bracket_ln, y_bracket_tick],
                seriestype = :line,
                legend = false,
                linecolor = "black",
                linealpha = 0.3,
                linewidth = 1.0)
    end

    xtick_pos = [x[end, 1]]
    Plots.plot!(plt,
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
            markersize = 0,
            linecolor = colors[:,1],
            linealpha = error_bar_alpha)
    Plots.plot!(plt,
            x[:,1],
            y[:,1],
            seriestype = :scatter,
            legend = false,
            xlims = x_limits,
            ylims = y_limits,
            label = labels[:,1],
            markershape = marker_shapes[:,1],
            markersize = marker_sizes[:,1],
            markercolor = colors[:,1],
            markeralpha = alphas[:,1],
            markerstrokecolor = colors[:,1],
            markerstrokealpha = 0.0)
            #= markerstrokealpha = alphas[:,1]) =#
    if length(error_bar_overflows_x[1]) > 0
        Plots.plot!(plt,
                error_bar_overflows_x[1],
                error_bar_overflows_y[1],
                seriestype = :scatter,
                legend = false,
                markershape = :utriangle,
                markersize = 2.0,
                markercolor = colors[:,1],
                markeralpha = error_bar_alpha,
                markerstrokecolor = colors[:,1],
                markerstrokealpha = error_bar_alpha)
    end
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
        #= y_position = y_limits[1] =#
        y_position = y_limits[1] - (x_label_offset * (y_limits[2] - y_limits[1]))
        annotate!(plt, x_position, y_position,
                  text(labels[1,1],
                       :center,
                       :top,
                       x_label_size,
                       rotation = 0),
                  annotation_clip = false)
    end
    for i in 2:ngroups
        push!(xtick_pos, x[end, i])
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
                markersize = 0,
                linecolor = colors[:,i],
                linealpha = error_bar_alpha)
        Plots.plot!(plt,
                x[:,i],
                y[:,i],
                seriestype = :scatter,
                legend = false,
                label = labels[:,i],
                markershape = marker_shapes[:,i],
                markersize = marker_sizes[:,i],
                markercolor = colors[:,i],
                markeralpha = alphas[:,i],
                markerstrokecolor = colors[:,i],
                markerstrokealpha = 0.0)
                #= markerstrokealpha = alphas[:,i]) =#
        if length(error_bar_overflows_x[i]) > 0
            Plots.plot!(plt,
                    error_bar_overflows_x[i],
                    error_bar_overflows_y[i],
                    seriestype = :scatter,
                    legend = false,
                    markershape = :utriangle,
                    markersize = 2.0,
                    markercolor = colors[:,i],
                    markeralpha = error_bar_alpha,
                    markerstrokecolor = colors[:,i],
                    markerstrokealpha = error_bar_alpha)
        end
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
            #= y_position = y_limits[1] =#
            y_position = y_limits[1] - (x_label_offset * (y_limits[2] - y_limits[1]))
            annotate!(plt, x_position, y_position,
                      text(labels[1,i],
                           :center,
                           :top,
                           x_label_size,
                           rotation = 0),
                      annotation_clip = false)
        end
    end
    Plots.plot!(plt, xticks = xtick_pos)
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
    xlabs = xlabels
    for i in eachindex(xlabs)
        xlabs[i] = "\\textrm{\\sffamily $(xlabs[i])}"
    end
    vln = StatsPlots.violin(
            xlabs,
            values,
            legend = false,
            fillcolor = fill_colors,
            linecolor = fill_colors,
            fillalpha = fill_alphas,
    )
    if include_dots
        StatsPlots.dotplot!(vln,
                xlabs,
                values,
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
        dot_legend = false,
        left_marker_shapes = nothing,
        left_marker_sizes = nothing,
        right_marker_shapes = nothing,
        right_marker_sizes = nothing
       )::Plots.Plot
    xlabs = xlabels
    for i in eachindex(xlabs)
        xlabs[i] = "\\textrm{\\sffamily $(xlabs[i])}"
    end
    if isnothing(left_marker_shapes)
        left_marker_shapes = fill(:circle, 1, length(left_values))
    end
    if isnothing(right_marker_shapes)
        right_marker_shapes = fill(:circle, 1, length(right_values))
    end
    if isnothing(left_marker_sizes)
        left_marker_sizes = fill(4.0, 1, length(left_values))
    end
    if isnothing(right_marker_sizes)
        right_marker_sizes = fill(4.0, 1, length(right_values))
    end
    vln = StatsPlots.violin(
            xlabs,
            left_values,
            legend = legend,
            side = :left,
            fillcolor = left_fill_colors,
            fillalpha = left_fill_alphas,
            label = left_labels
    )
    StatsPlots.dotplot!(vln,
            xlabs,
            left_values,
            legend = dot_legend,
            side = :left,
            markercolor = left_marker_colors,
            markeralpha = left_marker_alphas,
            markershape = left_marker_shapes,
            #= markerstrokecolor = left_marker_colors, =#
            #= markerstrokealpha = left_marker_alphas, =#
            markersize = left_marker_sizes,
            label = left_labels
    )
    StatsPlots.violin!(vln,
            xlabs,
            right_values,
            legend = legend,
            side = :right,
            fillcolor = right_fill_colors,
            fillalpha = right_fill_alphas,
            label = right_labels
    )
    StatsPlots.dotplot!(vln,
            xlabs,
            right_values,
            legend = dot_legend,
            side = :right,
            markercolor = right_marker_colors,
            markeralpha = right_marker_alphas,
            markershape = right_marker_shapes,
            #= markerstrokecolor = right_marker_colors, =#
            #= markerstrokealpha = right_marker_alphas, =#
            markersize = right_marker_sizes,
            label = right_labels
    )
    return vln
end

function get_split_dot_plot(
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
        legend = true)::Plots.Plot
    xlabs = xlabels
    for i in eachindex(xlabs)
        xlabs[i] = "\\textrm{\\sffamily $(xlabs[i])}"
    end
    vln = StatsPlots.dotplot!(
            xlabs,
            left_values,
            legend = legend,
            side = :left,
            markercolor = left_marker_colors,
            markeralpha = left_marker_alphas,
            label = left_labels
    )
    StatsPlots.dotplot!(vln,
            xlabs,
            right_values,
            legend = legend,
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
    #= identity_line = [extremes[1] - (xy_buffer * 20.0) =#
    #=                  extremes[2] + (xy_buffer * 20.0)] =#

    #= plt = Plots.plot(identity_line, identity_line, =#
    #=         seriestype = :line, =#
    #=         legend = false, =#
    #=         xlims = axis_limits, =#
    #=         ylims = axis_limits, =#
    #=         color = identity_line_color, =#
    #=         alpha = identity_line_alpha) =#
    #= Plots.plot!(plt, =#
    plt = Plots.plot(
            x,
            y,
            yerror = (y - y_lower, y_upper - y),
            seriestype = :scatter,
            link = :all,
            legend = false,
            xlims = axis_limits,
            ylims = axis_limits,
            markercolor = scatter_marker_col,
            markeralpha = scatter_marker_alpha,
            #= markerstrokecolor = error_bar_rgba, =#
            #= markerstrokecolor = scatter_marker_col, =#
            #= markerstrokealpha = scatter_marker_alpha, =#
            #= markerstrokealpha = 0.2, =#
            #= markershape = :circle, =#
            #= markerfillcolor = :transparent, =#
            #= markerfillalpha = 0.0, =#
            markerstrokecolor = scatter_marker_col,
            markerstrokealpha = error_bar_alpha,
            markersize = 0,
            #= markerstrokealpha = scatter_marker_alpha, =#
            #= linecolor = error_bar_rgba) =#
            linecolor = scatter_marker_col,
            linealpha = error_bar_alpha)
    Plots.plot!(plt,
            x,
            y,
            seriestype = :scatter,
            link = :all,
            legend = false,
            xlims = axis_limits,
            ylims = axis_limits,
            markercolor = scatter_marker_col,
            markeralpha = scatter_marker_alpha,
            markerstrokecolor = scatter_marker_col,
            markerstrokealpha = 0.0)
            #= markerstrokealpha = scatter_marker_alpha) =#
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
                push!(highlight_x, x[i])
                push!(highlight_y, y[i])
                push!(highlight_y_lower, y_lower[i])
                push!(highlight_y_upper, y_upper[i])
            end
        end
        if length(highlight_x) > 0
            Plots.plot!(plt,
                    highlight_x,
                    highlight_y,
                    yerror = (highlight_y - highlight_y_lower, highlight_y_upper - highlight_y),
                    seriestype = :scatter,
                    markercolor = highlight_color,
                    markeralpha = 1.0,
                    markersize = 0,
                    markerstrokealpha = 1.0,
                    markerstrokecolor = highlight_color,
                    linecolor = highlight_color,
                    linealpha = 1.0)
            Plots.plot!(plt,
                    highlight_x,
                    highlight_y,
                    seriestype = :scatter,
                    markercolor = highlight_color,
                    markeralpha = 1.0,
                    #= markerstrokealpha = 1.0, =#
                    markerstrokealpha = 0.0,
                    markerstrokecolor = highlight_color)
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

function add_identity_line!(
        plots::Plots.Plot...)
    for plt in plots
        xy_lims = get_shared_xy_limits(plt)
        identity_line = [xy_lims[1]
                         xy_lims[2]]
        Plots.plot!(plt, identity_line, identity_line,
                seriestype = :line,
                legend = false,
                color = identity_line_color,
                alpha = identity_line_alpha)
                #= zorder = 0) =#
    end
    return nothing
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

function mean_squared_error(
        x::Vector{Float64},
        y::Vector{Float64})::Float64
    @assert length(x) == length(y)
    sse = 0.0
    for i in 1:length(x)
        sse += (x[i] - y[i])^2
    end
    return sse / length(x)
end

function root_mean_squared_error(
        x::Vector{Float64},
        y::Vector{Float64})::Float64
    return sqrt(mean_squared_error(x, y))
end

function insert_lines(
        in_path::String,
        out_path::String,
        after::Regex,
        lines::Vector{String};
        target::Union{String, Regex, Nothing} = nothing,
        replacement::Union{String, Nothing} = nothing
        )
    open(in_path, "r") do istream
        open(out_path, "w") do ostream
            for raw_line in eachline(istream)
                line = chomp(raw_line)
                m = match(after, line)
                if (! isnothing(target)) & (! isnothing(replacement))
                    line = replace(line, target => replacement)
                end
                write(ostream, "$(line)\n")
                if ! isnothing(m)
                    for l in lines
                        write(ostream, "$(l)\n")
                    end
                end
            end
        end
    end
    return nothing
end

function tex_to_sans(
        path::String;
        target::Union{String, Regex, Nothing} = nothing,
        replacement::Union{String, Nothing} = nothing
        )::String;
    prefix, ext = splitext(path)
    out_path = "$(prefix)-sans$(ext)"
    insert_lines(path, out_path, sans_tex_insert, sans_tex_lines,
            target = target,
            replacement = replacement)
    return out_path
end

function run_latexmk(
        path::String,
        options::String = "-pdf")
    if isnothing(Sys.which("latexmk"))
        return nothing
    end
    dir_path = dirname(path)
    file_name = basename(path)
    cmd = `latexmk $options $file_name`
    if isnothing(dir_path) | (dir_path == "") | (dir_path == ".")
        run(cmd)
        return nothing
    end
    cd(dir_path) do
        run(cmd)
    end
    return nothing
end

function tex_to_pdf(
        path::String)
    run_latexmk(path, "-C")
    run_latexmk(path, "-pdf")
    return nothing
end

function tex_clean(
        path::String)
    run_latexmk(path, "-C")
    return nothing
end

function crop_pdf(
        path::String)
    if isnothing(Sys.which("pdfcrop"))
        return nothing
    end
    prefix, ext = splitext(path)
    out_path = "$(prefix)-cropped$(ext)"
    cmd = `pdfcrop $path $out_path`
    #= cmd = `cp $path $out_path` =#
    run(cmd)
    return out_path
end

function blank_pdf(
        path::String)
    if isnothing(Sys.which("convert"))
        return nothing
    end
    prefix, ext = splitext(path)
    out_path = "$(prefix)-blank$(ext)"
    cmd = `convert "$path" -threshold -1 -alpha off "$out_path"`
    #= cmd = `cp $path $out_path` =#
    run(cmd)
    return out_path
end

function process_tex(path::String;
        target::Union{String, Regex, Nothing} = nothing,
        replacement::Union{String, Nothing} = nothing
        )
    sans_path = tex_to_sans(path, target = target, replacement = replacement)
    tex_to_pdf(sans_path)
    sans_pdf = splitext(sans_path)[1] * ".pdf"
    cropped_pdf_path = crop_pdf(sans_pdf)
    if isnothing(cropped_pdf_path)
        return sans_pdf
    end
    tex_clean(sans_path)
    return cropped_pdf_path
end

function make_legends()
    plt = Plots.plot(
            [1],
            [2],
            seriestype = :scatter,
            legend = false,
            markershape = gen_shape,
            markersize = 8.0,
            markercolor = gen_col,
            markeralpha = gen_marker_alpha,
            markerstrokecolor = gen_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [2],
            [2],
            seriestype = :scatter,
            markershape = vo_gen_shape,
            markersize = 8.0,
            markercolor = vo_gen_col,
            markeralpha = vo_gen_marker_alpha,
            markerstrokecolor = vo_gen_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [1],
            [1],
            seriestype = :scatter,
            legend = false,
            axis = nothing,
            markershape = bif_shape,
            markersize = 8.0,
            markercolor = bif_col,
            markeralpha = bif_marker_alpha,
            markerstrokecolor = bif_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [2],
            [1],
            seriestype = :scatter,
            markershape = vo_bif_shape,
            markersize = 8.0,
            markercolor = vo_bif_col,
            markeralpha = vo_bif_marker_alpha,
            markerstrokecolor = vo_bif_col,
            markerstrokealpha = 0.0)
    annotate!(plt, 1, 0.2,
              text("All sites",
                   :center,
                   :top,
                   8),
              annotation_clip = false)
    annotate!(plt, 2, 0.18,
              text("Variable-only",
                   :center,
                   :top,
                   8),
              annotation_clip = false)
    annotate!(plt, 2.1, 2,
              text("Generalized model",
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 2.1, 1,
              text("Independent-bifurcating model",
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    Plots.plot!(plt, size = (135, 90), grid = false, ticks = false)
    Plots.xaxis!(plt, false)
    Plots.yaxis!(plt, false)
    plot_path = joinpath(ProjectUtil.RESULTS_DIR, "legend.tex")
    Plots.savefig(plt, plot_path)
    process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

    plt = Plots.plot(
            [1],
            [2],
            seriestype = :scatter,
            legend = false,
            markershape = gen_shape,
            markersize = 8.0,
            markercolor = gen_col,
            markeralpha = gen_marker_alpha,
            markerstrokecolor = gen_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [2],
            [2],
            seriestype = :scatter,
            legend = false,
            axis = nothing,
            markershape = bif_shape,
            markersize = 8.0,
            markercolor = bif_col,
            markeralpha = bif_marker_alpha,
            markerstrokecolor = bif_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [1],
            [1],
            seriestype = :scatter,
            markershape = vo_gen_shape,
            markersize = 8.0,
            markercolor = vo_gen_col,
            markeralpha = vo_gen_marker_alpha,
            markerstrokecolor = vo_gen_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [2],
            [1],
            seriestype = :scatter,
            markershape = vo_bif_shape,
            markersize = 8.0,
            markercolor = vo_bif_col,
            markeralpha = vo_bif_marker_alpha,
            markerstrokecolor = vo_bif_col,
            markerstrokealpha = 0.0)
    annotate!(plt, 1, 0.2,
              text(gen_model,
                   :center,
                   :top,
                   8),
              annotation_clip = false)
    annotate!(plt, 2, 0.18,
              text(bif_model,
                   :center,
                   :top,
                   8),
              annotation_clip = false)
    annotate!(plt, 2.15, 2,
              text("All sites",
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 2.15, 1,
              text("Variable sites",
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    Plots.plot!(plt, size = (95, 90), grid = false, ticks = false)
    Plots.xaxis!(plt, false)
    Plots.yaxis!(plt, false)
    plot_path = joinpath(ProjectUtil.RESULTS_DIR, "legend-short.tex")
    Plots.savefig(plt, plot_path)
    process_tex(plot_path, target = axis_pattern, replacement = axis_replace)


    plt = Plots.plot(
            [1],
            [4],
            seriestype = :scatter,
            legend = false,
            markershape = gen_shape,
            markersize = 6.0,
            markercolor = gen_col,
            markeralpha = gen_marker_alpha,
            markerstrokecolor = gen_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [1],
            [3],
            seriestype = :scatter,
            markershape = vo_gen_shape,
            markersize = 6.0,
            markercolor = vo_gen_col,
            markeralpha = vo_gen_marker_alpha,
            markerstrokecolor = vo_gen_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [1],
            [2],
            seriestype = :scatter,
            legend = false,
            axis = nothing,
            markershape = bif_shape,
            markersize = 6.0,
            markercolor = bif_col,
            markeralpha = bif_marker_alpha,
            markerstrokecolor = bif_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [1],
            [1],
            seriestype = :scatter,
            markershape = vo_bif_shape,
            markersize = 6.0,
            markercolor = vo_bif_col,
            markeralpha = vo_bif_marker_alpha,
            markerstrokecolor = vo_bif_col,
            markerstrokealpha = 0.0)
    annotate!(plt, 0.85, 4.9,
              text(LaTeXString("\\textbf{Sites}"),
                   :right,
                   :bottom,
                   8),
              annotation_clip = false)
    annotate!(plt, 0.85, 4,
              text("All",
                   :right,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 0.85, 3,
              text("Variable",
                   :right,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 0.85, 2,
              text("All",
                   :right,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 0.85, 1,
              text("Variable",
                   :right,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 1.15, 4.9,
              text(LaTeXString("\\textbf{Model}"),
                   :left,
                   :bottom,
                   8),
              annotation_clip = false)
    annotate!(plt, 1.15, 4,
              text("Generalized",
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 1.15, 2,
              text("Independent-bifurcating",
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    Plots.plot!(plt, size = (100, 115), grid = false, ticks = false)
    Plots.xaxis!(plt, false)
    Plots.yaxis!(plt, false)
    plot_path = joinpath(ProjectUtil.RESULTS_DIR, "legend-vertical.tex")
    Plots.savefig(plt, plot_path)
    process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

    plt = Plots.plot(
            [1],
            [4],
            seriestype = :scatter,
            legend = false,
            markershape = gen_shape,
            markersize = 6.0,
            markercolor = gen_col,
            markeralpha = gen_marker_alpha,
            markerstrokecolor = gen_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [1],
            [3],
            seriestype = :scatter,
            markershape = vo_gen_shape,
            markersize = 6.0,
            markercolor = vo_gen_col,
            markeralpha = vo_gen_marker_alpha,
            markerstrokecolor = vo_gen_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [1],
            [2],
            seriestype = :scatter,
            legend = false,
            axis = nothing,
            markershape = bif_shape,
            markersize = 6.0,
            markercolor = bif_col,
            markeralpha = bif_marker_alpha,
            markerstrokecolor = bif_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [1],
            [1],
            seriestype = :scatter,
            markershape = vo_bif_shape,
            markersize = 6.0,
            markercolor = vo_bif_col,
            markeralpha = vo_bif_marker_alpha,
            markerstrokecolor = vo_bif_col,
            markerstrokealpha = 0.0)
    annotate!(plt, 1.15, 4,
              text(LaTeXString("$(gen_model), all sites"),
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 1.15, 3,
              text(LaTeXString("$(gen_model), variable sites"),
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 1.15, 2,
              text(LaTeXString("$(bif_model), all sites"),
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 1.15, 1,
              text(LaTeXString("$(bif_model), variable sites"),
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    Plots.plot!(plt, size = (100, 114), grid = false, ticks = false)
    Plots.xaxis!(plt, false)
    Plots.yaxis!(plt, false)
    plot_path = joinpath(ProjectUtil.RESULTS_DIR, "legend-vertical-short.tex")
    Plots.savefig(plt, plot_path)
    process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

    plt = Plots.plot(
            [0],
            [0],
            seriestype = :scatter,
            legend = false,
            markershape = gen_shape,
            markersize = 6.0,
            markercolor = gen_col,
            markeralpha = gen_marker_alpha,
            markerstrokecolor = gen_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [0.8],
            [0],
            seriestype = :scatter,
            markershape = vo_gen_shape,
            markersize = 6.0,
            markercolor = vo_gen_col,
            markeralpha = vo_gen_marker_alpha,
            markerstrokecolor = vo_gen_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [1.85],
            [0],
            seriestype = :scatter,
            legend = false,
            axis = nothing,
            markershape = bif_shape,
            markersize = 6.0,
            markercolor = bif_col,
            markeralpha = bif_marker_alpha,
            markerstrokecolor = bif_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [2.65],
            [0],
            seriestype = :scatter,
            markershape = vo_bif_shape,
            markersize = 6.0,
            markercolor = vo_bif_col,
            markeralpha = vo_bif_marker_alpha,
            markerstrokecolor = vo_bif_col,
            markerstrokealpha = 0.0)
    annotate!(plt, 0.03, 0,
              text(LaTeXString("$(gen_model), all sites"),
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 0.83, 0,
              text(LaTeXString("$(gen_model), variable sites"),
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 1.88, 0,
              text(LaTeXString("$(bif_model), all sites"),
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 2.68, 0,
              text(LaTeXString("$(bif_model), variable sites"),
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    Plots.plot!(plt, size = (380, 80), grid = false, ticks = false)
    Plots.xaxis!(plt, false)
    Plots.yaxis!(plt, false)
    plot_path = joinpath(ProjectUtil.RESULTS_DIR, "legend-horizontal-short.tex")
    Plots.savefig(plt, plot_path)
    process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

    plt = Plots.plot(
            [0],
            [0],
            seriestype = :scatter,
            legend = false,
            markershape = gen_shape,
            markersize = 6.0,
            markercolor = gen_col,
            markeralpha = gen_marker_alpha,
            markerstrokecolor = gen_col,
            markerstrokealpha = 0.0)
    Plots.plot!(plt,
            [0.8],
            [0],
            seriestype = :scatter,
            markershape = bif_shape,
            markersize = 6.0,
            markercolor = bif_col,
            markeralpha = bif_marker_alpha,
            markerstrokecolor = bif_col,
            markerstrokealpha = 0.0)
    annotate!(plt, 0.03, 0,
              text(LaTeXString("Generalized"),
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    annotate!(plt, 0.83, 0,
              text(LaTeXString("Bifurcating"),
                   :left,
                   :vcenter,
                   8),
              annotation_clip = false)
    Plots.plot!(plt, size = (170, 80), grid = false, ticks = false)
    Plots.xaxis!(plt, false)
    Plots.yaxis!(plt, false)
    plot_path = joinpath(ProjectUtil.RESULTS_DIR, "legend-horizontal-short-all-sites.tex")
    Plots.savefig(plt, plot_path)
    process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
    return nothing
end

function main_cli()::Cint
    if ! Base.Filesystem.ispath(ProjectUtil.RESULTS_DIR)
        Base.Filesystem.mkdir(ProjectUtil.RESULTS_DIR)
    end

    write(stdout, "Plotting backend: $(backend())\n")

    make_legends()

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

    shared_div_results = SimResults(suffix = "-true-shared-height-probs")
    polytomy_results = SimResults(suffix = "-true-polytomy-probs")

    locus_sizes = [1, 100]

    sum_stats = DataFrame()
    bools = [true, false]
    for locus_size in locus_sizes
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
                        tree_len_rmse = root_mean_squared_error(
                                get_floats(results,
                                    sim_fixed,
                                    sim_generalized,
                                    analysis_generalized,
                                    only_var_sites,
                                    :tree_length_true,
                                    locus_size),
                                get_floats(results,
                                    sim_fixed,
                                    sim_generalized,
                                    analysis_generalized,
                                    only_var_sites,
                                    :tree_length_mean,
                                    locus_size),
                                )
                        root_height_perc = get_floats(results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :root_height_true_percentile,
                                locus_size)
                        root_height_cover = Statistics.mean(
                                0.025 .<= root_height_perc .<= 0.975)
                        root_height_rmse = root_mean_squared_error(
                                get_floats(results,
                                    sim_fixed,
                                    sim_generalized,
                                    analysis_generalized,
                                    only_var_sites,
                                    :root_height_true,
                                    locus_size),
                                get_floats(results,
                                    sim_fixed,
                                    sim_generalized,
                                    analysis_generalized,
                                    only_var_sites,
                                    :root_height_mean,
                                    locus_size),
                                )
                        root_pop_size_perc = get_floats(results,
                                sim_fixed,
                                sim_generalized,
                                analysis_generalized,
                                only_var_sites,
                                :root_pop_size_true_percentile,
                                locus_size)
                        root_pop_size_cover = Statistics.mean(
                                0.025 .<= root_pop_size_perc .<= 0.975)
                        root_pop_size_rmse = root_mean_squared_error(
                                get_floats(results,
                                    sim_fixed,
                                    sim_generalized,
                                    analysis_generalized,
                                    only_var_sites,
                                    :root_pop_size_true,
                                    locus_size),
                                get_floats(results,
                                    sim_fixed,
                                    sim_generalized,
                                    analysis_generalized,
                                    only_var_sites,
                                    :root_pop_size_mean,
                                    locus_size),
                                )

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
                                rmse_tree_length = tree_len_rmse,
                                coverage_root_height = root_height_cover,
                                rmse_root_height = root_height_rmse,
                                coverage_root_pop_size = root_pop_size_cover,
                                rmse_root_pop_size = root_pop_size_rmse,
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
    end

    sum_results = SimResults(sum_stats)

    for locus_size in locus_sizes
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
    end

    for locus_size in locus_sizes
        locus_prefix = ""
        if locus_size > 1
            locus_prefix = "locus-$locus_size-"
        end

        true_shared_div_probs = get_floats(shared_div_results,
                false,
                true,
                true,
                false,
                :shared_height_prob,
                locus_size)
        vo_true_shared_div_probs = get_floats(shared_div_results,
                false,
                true,
                true,
                true,
                :shared_height_prob,
                locus_size)
        true_polytomy_probs = get_floats(polytomy_results,
                false,
                true,
                true,
                false,
                :polytomy_node_prob,
                locus_size)
        vo_true_polytomy_probs = get_floats(polytomy_results,
                false,
                true,
                true,
                true,
                :polytomy_node_prob,
                locus_size)

        v_shared_div_probs = get_split_violin_plot(
                true_shared_div_probs,
                vo_true_shared_div_probs,
                xlabels = [ "True shared divs" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_marker_shapes = gen_shape,
                left_marker_sizes = gen_marker_size - 1,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_marker_shapes = vo_gen_shape,
                right_marker_sizes = vo_gen_marker_size - 1,
                right_labels = [ "Variable sites" ],
                legend = false,
                dot_legend = false)
        Plots.plot!(v_shared_div_probs, size = (190, 220))
        Plots.ylims!(v_shared_div_probs, (-0.02, 1.02))
        Plots.ylabel!(v_shared_div_probs, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-gen-true-shared-height-probs.tex")
        Plots.savefig(v_shared_div_probs, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v_poly_probs = get_split_violin_plot(
                true_polytomy_probs,
                vo_true_polytomy_probs,
                xlabels = [ "True polytomies" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_marker_shapes = gen_shape,
                left_marker_sizes = gen_marker_size - 1,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_marker_shapes = vo_gen_shape,
                right_marker_sizes = vo_gen_marker_size - 1,
                right_labels = [ "Variable sites" ],
                legend = false,
                dot_legend = false)
        Plots.plot!(v_poly_probs, size = (190, 220))
        Plots.ylims!(v_poly_probs, (-0.02, 1.02))
        Plots.ylabel!(v_poly_probs, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-gen-true-polytomy-probs.tex")
        Plots.savefig(v_poly_probs, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        # This doesn't work because the vectors have different lengths (i.e.,
        # there are different numbers of shared divs and multifurctations
        #= v = get_split_violin_plot( =#
        #=         hcat( =#
        #=              true_shared_div_probs, =#
        #=              true_polytomy_probs), =#
        #=         hcat( =#
        #=              vo_true_shared_div_probs, =#
        #=              vo_true_polytomy_probs), =#
        #=         xlabels = [ "True shared divs" "True polytomies" ], =#
        #=         left_fill_colors = gen_col, =#
        #=         left_marker_colors = gen_col, =#
        #=         left_fill_alphas = gen_fill_alpha, =#
        #=         left_marker_alphas = gen_marker_alpha, =#
        #=         left_marker_shapes = gen_shape, =#
        #=         left_marker_sizes = gen_marker_size - 1, =#
        #=         left_labels = [ "All sites" ], =#
        #=         right_fill_colors = vo_gen_col, =#
        #=         right_marker_colors = vo_gen_col, =#
        #=         right_fill_alphas = vo_gen_fill_alpha, =#
        #=         right_marker_alphas = vo_gen_marker_alpha, =#
        #=         right_marker_shapes = vo_gen_shape, =#
        #=         right_marker_sizes = vo_gen_marker_size - 1, =#
        #=         right_labels = [ "Variable sites" ], =#
        #=         legend = false, =#
        #=         dot_legend = false) =#
        #= Plots.plot!(v, size = (350, 220)) =#
        #= Plots.ylims!(v, (-0.02, 1.02)) =#
        #= Plots.ylabel!(v, "Posterior probability") =#
        #= plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-gen-true-event-probs.tex") =#
        #= Plots.savefig(v, plot_path) =#
        #= process_tex(plot_path, target = axis_pattern, replacement = axis_replace) =#


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
        #= Plots.savefig(v, =#
        #=         joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-root-vln.pdf")) =#
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-root-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-root-all-sites-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        h = Plots.histogram(root_node_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, yticks = false, grid = :x, xticks = (0:0.5:1), tick_direction = :out)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-root-all-sites-hist.tex")
        Plots.savefig(h, plot_path)
        pdf_path = process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
        blank_pdf(pdf_path)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-789-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-789-all-sites-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        h = Plots.histogram(node_789_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, yticks = false, grid = :x, xticks = (0:0.5:1), tick_direction = :out)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-789-all-sites-hist.tex")
        Plots.savefig(h, plot_path)
        pdf_path = process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
        blank_pdf(pdf_path)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-456-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-456-all-sites-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        h = Plots.histogram(node_456_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, yticks = false, grid = :x, xticks = (0:0.5:1), tick_direction = :out)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-456-all-sites-hist.tex")
        Plots.savefig(h, plot_path)
        pdf_path = process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
        blank_pdf(pdf_path)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-12-3-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-12-3-all-sites-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        h = Plots.histogram(node_12_3_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, yticks = false, grid = :x, xticks = (0:0.5:1), tick_direction = :out)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-node-probs-12-3-all-sites-hist.tex")
        Plots.savefig(h, plot_path)
        pdf_path = process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
        blank_pdf(pdf_path)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-height-probs-12-789-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-height-probs-12-789-all-sites-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        h = Plots.histogram(height_12_789_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, yticks = false, grid = :x, xticks = (0:0.5:1), tick_direction = :out)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-height-probs-12-789-all-sites-hist.tex")
        Plots.savefig(h, plot_path)
        pdf_path = process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
        blank_pdf(pdf_path)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-height-probs-123-456-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-height-probs-123-456-all-sites-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        h = Plots.histogram(height_123_456_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, yticks = false, grid = :x, xticks = (0:0.5:1), tick_direction = :out)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-height-probs-123-456-all-sites-hist.tex")
        Plots.savefig(h, plot_path)
        pdf_path = process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
        blank_pdf(pdf_path)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-split-probs-12-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-split-probs-12-all-sites-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        h = Plots.histogram(split_12_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, yticks = false, grid = :x, xticks = (0:0.5:1), tick_direction = :out)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-split-probs-12-all-sites-hist.tex")
        Plots.savefig(h, plot_path)
        pdf_path = process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
        blank_pdf(pdf_path)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-topo-probs-vln.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        h = Plots.histogram(true_topo_probs,
                bins = prob_bins,
                normalize = :probability,
                fillcolor = gen_col,
                linecolor = gen_col,
                fillalpha = 1.0,
                #= xlabel = "Posterior probability", =#
                #= ylabel = "Frequency", =#
                legend = false)
        Plots.plot!(h, size = (140, 120), yaxis = false, yticks = false, grid = :x, xticks = (0:0.5:1), tick_direction = :out)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-topo-probs-all-sites-hist.tex")
        Plots.savefig(h, plot_path)
        pdf_path = process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
        blank_pdf(pdf_path)

        v = get_split_violin_plot(
                hcat(
                     height_123_456_probs,
                     height_12_789_probs,
                     root_node_probs,
                     node_789_probs,
                     node_456_probs,
                     node_12_3_probs,
                     split_12_probs),
                     #= true_topo_probs), =#
                hcat(
                     vo_height_123_456_probs,
                     vo_height_12_789_probs,
                     vo_root_node_probs,
                     vo_node_789_probs,
                     vo_node_456_probs,
                     vo_node_12_3_probs,
                     vo_split_12_probs),
                     #= vo_true_topo_probs), =#
                xlabels = [ "\$\\tau_2\$" "\$\\tau_1\$" "\$t_5\$" "\$t_4\$" "\$t_3\$" "\$t_2\$" "\$t_1\$" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_marker_shapes = gen_shape,
                left_marker_sizes = gen_marker_size,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_marker_shapes = vo_gen_shape,
                right_marker_sizes = vo_gen_marker_size,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (650, 220), xtickfontsize = 16)
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-true-tree-probs.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

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

        # plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-gen-true-tree-probs.tex")
        # write(stdout, "Writing to $(plot_path)\n")
        # Plots.savefig(vln_grid, plot_path)
        # process_tex(plot_path, target = axis_pattern, replacement = axis_replace)


        gen_max_456_subsplit_prob = get_floats(results, true, true, true, false, :max_456_subsplit_prob, locus_size)
        vo_gen_max_456_subsplit_prob = get_floats(results, true, true, true, true, :max_456_subsplit_prob, locus_size)
        bif_max_456_subsplit_prob = get_floats(results, true, true, false, false, :max_456_subsplit_prob, locus_size)
        vo_bif_max_456_subsplit_prob = get_floats(results, true, true, false, true, :max_456_subsplit_prob, locus_size)

        v = get_split_violin_plot(
                hcat(gen_max_456_subsplit_prob, vo_gen_max_456_subsplit_prob),
                hcat(bif_max_456_subsplit_prob, vo_bif_max_456_subsplit_prob),
                #= xlabels = ["All sites" "Variable sites"], =#
                xlabels = [ "Splitting \$t_3\$" "\\ Splitting \$t_3\$\\ " ],
                #= xlabels = [ latexstring("\\begin{tabular}{c} {\\Large \$t_3\$} \\\\ All \\\\ sites \\end{tabular}") latexstring( "\\begin{tabular}{c} {\\Large \$t_3\$} \\\\ Variable \\\\ sites \\end{tabular}") ], =#
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
        #= Plots.plot!(v, size = (250, 220), xtickfontsize = 16) =#
        Plots.plot!(v, size = (250, 220))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-max-456-subsplit-probs.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)


        gen_max_789_subsplit_prob = get_floats(results, true, true, true, false, :max_789_subsplit_prob, locus_size)
        vo_gen_max_789_subsplit_prob = get_floats(results, true, true, true, true, :max_789_subsplit_prob, locus_size)
        bif_max_789_subsplit_prob = get_floats(results, true, true, false, false, :max_789_subsplit_prob, locus_size)
        vo_bif_max_789_subsplit_prob = get_floats(results, true, true, false, true, :max_789_subsplit_prob, locus_size)

        v = get_split_violin_plot(
                hcat(gen_max_789_subsplit_prob, vo_gen_max_789_subsplit_prob),
                hcat(bif_max_789_subsplit_prob, vo_bif_max_789_subsplit_prob),
                #= xlabels = ["All sites" "Variable sites"], =#
                xlabels = [ "Splitting \$t_4\$" "\\ Splitting \$t_4\$\\ " ],
                #= xlabels = [ latexstring("\\begin{tabular}{c} {\\Large \$t_4\$} \\\\ All \\\\ sites \\end{tabular}") latexstring( "\\begin{tabular}{c} {\\Large \$t_4\$} \\\\ Variable \\\\ sites \\end{tabular}") ], =#
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
        #= Plots.plot!(v, size = (250, 220), xtickfontsize = 16) =#
        Plots.plot!(v, size = (250, 220))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-max-789-subsplit-probs.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        fixed_gen_gen_max_wrong_root_prob = get_floats(results, true, true, true, false, :max_wrong_root_prob, locus_size)
        vo_fixed_gen_gen_max_wrong_root_prob = get_floats(results, true, true, true, true, :max_wrong_root_prob, locus_size)
        fixed_gen_bif_max_wrong_root_prob = get_floats(results, true, true, false, false, :max_wrong_root_prob, locus_size)
        vo_fixed_gen_bif_max_wrong_root_prob = get_floats(results, true, true, false, true, :max_wrong_root_prob, locus_size)

        v = get_split_violin_plot(
                hcat(fixed_gen_gen_max_wrong_root_prob, vo_fixed_gen_gen_max_wrong_root_prob),
                hcat(fixed_gen_bif_max_wrong_root_prob, vo_fixed_gen_bif_max_wrong_root_prob),
                xlabels = [ "Wrong \$t_5\$" "\\ Wrong \$t_5\$\\ " ],
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
        Plots.plot!(v, size = (250, 220))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-max-wrong-root-probs.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_split_violin_plot(
                hcat(gen_max_456_subsplit_prob, vo_gen_max_456_subsplit_prob, gen_max_789_subsplit_prob, vo_gen_max_789_subsplit_prob),
                hcat(bif_max_456_subsplit_prob, vo_bif_max_456_subsplit_prob, bif_max_789_subsplit_prob, vo_bif_max_789_subsplit_prob),
                xlabels = [ "Splitting \$t_3\$" "\\ Splitting \$t_3\$\\ " "Splitting \$t_4\$" "\\ Splitting \$t_4\$\\ " ],
                left_fill_colors =   [gen_col vo_gen_col gen_col vo_gen_col],
                left_marker_colors = [gen_col vo_gen_col gen_col vo_gen_col],
                left_fill_alphas = [gen_fill_alpha vo_gen_fill_alpha gen_fill_alpha vo_gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha vo_gen_marker_alpha gen_marker_alpha vo_gen_marker_alpha],
                left_labels = "Generalized",
                right_fill_colors = [bif_col vo_bif_col bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha bif_marker_alpha vo_bif_marker_alpha],
                right_labels = "Bifurcating",
                legend = false,
                dot_legend = false)
        #= Plots.plot!(v, size = (250, 220), xtickfontsize = 16) =#
        Plots.plot!(v, size = (450, 220), xtickfontsize = 12)
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-max-poly-subsplit-probs.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_split_violin_plot(
                hcat(fixed_gen_gen_max_wrong_root_prob, vo_fixed_gen_gen_max_wrong_root_prob, gen_max_456_subsplit_prob, vo_gen_max_456_subsplit_prob, gen_max_789_subsplit_prob, vo_gen_max_789_subsplit_prob),
                hcat(fixed_gen_bif_max_wrong_root_prob, vo_fixed_gen_bif_max_wrong_root_prob, bif_max_456_subsplit_prob, vo_bif_max_456_subsplit_prob, bif_max_789_subsplit_prob, vo_bif_max_789_subsplit_prob),
                xlabels = [ "Wrong \$t_5\$" "\\ Wrong \$t_5\$\\ " "Splitting \$t_3\$" "\\ Splitting \$t_3\$\\ " "Splitting \$t_4\$" "\\ Splitting \$t_4\$\\ " ],
                left_fill_colors =   [gen_col vo_gen_col gen_col vo_gen_col gen_col vo_gen_col],
                left_marker_colors = [gen_col vo_gen_col gen_col vo_gen_col gen_col vo_gen_col],
                left_fill_alphas = [gen_fill_alpha vo_gen_fill_alpha gen_fill_alpha vo_gen_fill_alpha gen_fill_alpha vo_gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha vo_gen_marker_alpha gen_marker_alpha vo_gen_marker_alpha gen_marker_alpha vo_gen_marker_alpha],
                left_marker_shapes = [gen_shape vo_gen_shape gen_shape vo_gen_shape gen_shape vo_gen_shape],
                left_marker_sizes = [gen_marker_size vo_gen_marker_size gen_marker_size vo_gen_marker_size gen_marker_size vo_gen_marker_size],
                left_labels = "Generalized",
                right_fill_colors = [bif_col vo_bif_col bif_col vo_bif_col bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col bif_col vo_bif_col bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha bif_fill_alpha vo_bif_fill_alpha bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha bif_marker_alpha vo_bif_marker_alpha bif_marker_alpha vo_bif_marker_alpha],
                right_marker_shapes = [bif_shape vo_bif_shape bif_shape vo_bif_shape bif_shape vo_bif_shape],
                right_marker_sizes = [bif_marker_size vo_bif_marker_size bif_marker_size vo_bif_marker_size bif_marker_size vo_bif_marker_size],
                right_labels = "Bifurcating",
                legend = false,
                dot_legend = false)
        Plots.plot!(v, size = (650, 220), xtickfontsize = 12)
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-max-wrong-root-poly-subsplit-probs.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_split_violin_plot(
                hcat(gen_max_456_subsplit_prob, gen_max_789_subsplit_prob),
                hcat(bif_max_456_subsplit_prob, bif_max_789_subsplit_prob),
                #= xlabels = [ "Node 456" "Node 789" ], =#
                xlabels = [ "\$t_3\$" "\$t_4\$" ],
                left_fill_colors = [ gen_col gen_col ],
                left_marker_colors = [gen_col gen_col],
                left_fill_alphas = [gen_fill_alpha gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha gen_marker_alpha],
                left_marker_shapes = [ gen_shape gen_shape ],
                left_marker_sizes = [ gen_marker_size gen_marker_size],
                left_labels = "Generalized",
                right_fill_colors = [bif_col bif_col],
                right_marker_colors = [bif_col bif_col],
                right_fill_alphas = [bif_fill_alpha bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha bif_marker_alpha],
                right_marker_shapes = [ bif_shape bif_shape ],
                right_marker_sizes = [ bif_marker_size bif_marker_size],
                right_labels = "Bifurcating",
                legend = :top,
                dot_legend = false)
        Plots.plot!(v, size = (250, 220), xtickfontsize = 16)
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-max-poly-subsplit-probs-all-sites.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v456 = get_split_violin_plot(
                gen_max_456_subsplit_prob,
                bif_max_456_subsplit_prob,
                xlabels = [ "Splitting \$t_3\$" ],
                left_fill_colors = [ gen_col ],
                left_marker_colors = [gen_col],
                left_fill_alphas = [gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha],
                left_marker_shapes = [gen_shape],
                left_marker_sizes = [gen_marker_size],
                left_labels = "Generalized",
                right_fill_colors = [bif_col],
                right_marker_colors = [bif_col],
                right_fill_alphas = [bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha],
                right_marker_shapes = [bif_shape],
                right_marker_sizes = [bif_marker_size],
                right_labels = "Bifurcating",
                legend = :top,
                dot_legend = false)

        v789 = get_split_violin_plot(
                gen_max_789_subsplit_prob,
                bif_max_789_subsplit_prob,
                xlabels = [ "Splitting \$t_4\$" ],
                left_fill_colors = [ gen_col ],
                left_marker_colors = [gen_col],
                left_fill_alphas = [gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha],
                left_marker_shapes = [gen_shape],
                left_marker_sizes = [gen_marker_size],
                left_labels = "Generalized",
                right_fill_colors = [bif_col],
                right_marker_colors = [bif_col],
                right_fill_alphas = [bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha],
                right_marker_shapes = [bif_shape],
                right_marker_sizes = [bif_marker_size],
                right_labels = "Bifurcating",
                legend = :top,
                dot_legend = false)

        # Add to grid to link y-axes
        g = Plots.plot(
                v456,
                v789,
                layout = (1, 2),
                legend = false,
                link = :y, # :none, :x, :y, :both, :all
        )
        Plots.plot!(v456, size = (220, 280))
        Plots.ylabel!(v456, "Posterior probability")
        Plots.plot!(v789, size = (220, 280))
        Plots.ylabel!(v789, "Posterior probability")

        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-max-456-subsplit-probs-all-sites.tex")
        Plots.savefig(v456, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-max-789-subsplit-probs-all-sites.tex")
        Plots.savefig(v789, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)


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
                left_fill_colors = [ gen_col  vo_gen_col ],
                left_marker_colors = [ gen_col  vo_gen_col ],
                left_fill_alphas = [ gen_fill_alpha  vo_gen_fill_alpha ],
                left_marker_alphas = [ gen_marker_alpha  vo_gen_marker_alpha ],
                left_labels = [ "Generalized" ],
                right_fill_colors = [bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha],
                right_labels = [ "Bifurcating" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (1200, 300))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-root-node-probs.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_split_violin_plot(
                hcat(bif_gen_true_topo_probs,
                     vo_bif_gen_true_topo_probs),
                hcat(bif_bif_true_topo_probs,
                     vo_bif_bif_true_topo_probs),
                #= xlabels = [ "All sites" "Variable sites" ], =#
                xlabels = [ "True topology" "\\ True topology\\ " ],
                #= xlabels = [ latexstring("\\begin{tabular}{c} True topology \\\\ All sites \\end{tabular}") latexstring( "\\begin{tabular}{c} True topology \\\\ Variable sites \\end{tabular}") ], =#
                left_fill_colors = [ gen_col  vo_gen_col ],
                left_marker_colors = [ gen_col  vo_gen_col ],
                left_fill_alphas = [ gen_fill_alpha  vo_gen_fill_alpha ],
                left_marker_alphas = [ gen_marker_alpha  vo_gen_marker_alpha ],
                left_marker_shapes = [ gen_shape vo_gen_shape ],
                left_marker_sizes = [ gen_marker_size vo_gen_marker_size ],
                left_labels = [ "Generalized" ],
                right_fill_colors = [bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha],
                right_marker_shapes = [ bif_shape vo_bif_shape ],
                right_marker_sizes = [ bif_marker_size vo_bif_marker_size ],
                right_labels = [ "Bifurcating" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (210, 220), xtickfontsize = 8)
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-true-topo-probs.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_split_violin_plot(
                bif_gen_true_topo_probs,
                bif_bif_true_topo_probs,
                #= xlabels = [ LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model (true model) \\end{tabular}") ], =#
                xlabels = [ "True topology" ],
                left_fill_colors = [ gen_col ],
                left_marker_colors = [ gen_col ],
                left_fill_alphas = [ gen_fill_alpha ],
                left_marker_alphas = [ gen_marker_alpha ],
                left_marker_shapes = [gen_shape],
                left_marker_sizes = [gen_marker_size],
                left_labels = [ "Generalized" ],
                right_fill_colors = [bif_col],
                right_marker_colors = [bif_col],
                right_fill_alphas = [bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha],
                right_marker_shapes = [bif_shape],
                right_marker_sizes = [bif_marker_size],
                right_labels = [ "Bifurcating" ],
                legend = :top,
                dot_legend = false)
        Plots.plot!(v, size = (180, 220))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-true-topo-probs-all-sites.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_split_violin_plot(
                hcat(bif_gen_true_split_prob_mean,
                     vo_bif_gen_true_split_prob_mean),
                hcat(bif_bif_true_split_prob_mean,
                     vo_bif_bif_true_split_prob_mean),
                xlabels = [ "All sites" "Variable sites" ],
                left_fill_colors = [ gen_col  vo_gen_col ],
                left_marker_colors = [ gen_col  vo_gen_col ],
                left_fill_alphas = [ gen_fill_alpha  vo_gen_fill_alpha ],
                left_marker_alphas = [ gen_marker_alpha  vo_gen_marker_alpha ],
                left_labels = [ "Generalized" ],
                right_fill_colors = [bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha],
                right_labels = [ "Bifurcating" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (250, 220))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-true-split-prob-means.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_split_violin_plot(
                bif_gen_true_split_prob_mean,
                bif_bif_true_split_prob_mean,
                xlabels = [ latexstring("\\begin{tabular}{c} independent bifurcating \\\\ model (true model) \\end{tabular}") ],
                left_fill_colors = [ gen_col ],
                left_marker_colors = [ gen_col ],
                left_fill_alphas = [ gen_fill_alpha ],
                left_marker_alphas = [ gen_marker_alpha ],
                left_marker_shapes = [gen_shape],
                left_marker_sizes = [gen_marker_size],
                left_labels = [ "Generalized" ],
                right_fill_colors = [bif_col],
                right_marker_colors = [bif_col],
                right_fill_alphas = [bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha],
                right_marker_shapes = [bif_shape],
                right_marker_sizes = [bif_marker_size],
                right_labels = [ "Bifurcating" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (1200, 300))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-true-split-prob-means-all-sites.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_split_violin_plot(
                hcat(
                     bif_gen_nd_123_456_789_prob,
                     bif_gen_ht_456_789_prob,
                     bif_gen_nd_4_5_6_prob,
                     bif_gen_ht_45_789_prob,
                     bif_gen_nd_1_2_3_prob,
                     bif_gen_ht_12_78_prob,
                    ),
                hcat(
                     vo_bif_gen_nd_123_456_789_prob,
                     vo_bif_gen_ht_456_789_prob,
                     vo_bif_gen_nd_4_5_6_prob,
                     vo_bif_gen_ht_45_789_prob,
                     vo_bif_gen_nd_1_2_3_prob,
                     vo_bif_gen_ht_12_78_prob,
                    ),
                xlabels = [ "\$t_7 = t_8\$" "\$t_4 = t_6\$" "\$t_3 = t_4\$" "\$t_3 = t_6\$" "\$t_1 = t_2\$" "\$t_1 = t_5\$" ],
                left_fill_colors = gen_col,
                left_marker_colors = gen_col,
                left_fill_alphas = gen_fill_alpha,
                left_marker_alphas = gen_marker_alpha,
                left_marker_shapes = gen_shape,
                left_marker_sizes = gen_marker_size,
                left_labels = [ "All sites" ],
                right_fill_colors = vo_gen_col,
                right_marker_colors = vo_gen_col,
                right_fill_alphas = vo_gen_fill_alpha,
                right_marker_alphas = vo_gen_marker_alpha,
                right_marker_shapes = vo_gen_shape,
                right_marker_sizes = vo_gen_marker_size,
                right_labels = [ "Variable sites" ],
                legend = :best,
                dot_legend = false)
        Plots.plot!(v, size = (650, 200), xtickfontsize = 16)
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-gen-wrong-probs.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_violin_plot(
                hcat(bif_gen_ht_12_78_prob,
                     bif_gen_ht_45_789_prob,
                     bif_gen_ht_456_789_prob,
                     bif_gen_nd_1_2_3_prob,
                     bif_gen_nd_4_5_6_prob,
                     bif_gen_nd_123_456_789_prob),
                xlabels = [ "\$t_1 = t_5\$" "\$t_3 = t_6\$" "\$t_4 = t_6\$" "\$t_1 = t_2\$" "\$t_3 = t_4\$" "\$t_7 = t_8\$" ],
                fill_colors = gen_col,
                marker_colors = gen_col,
                fill_alphas = gen_fill_alpha,
                marker_alphas = gen_marker_alpha,
                include_dots = true)
        Plots.plot!(v, size = (800, 180))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-gen-wrong-probs-all-sites.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)


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
        mwu_test = HypothesisTests.MannWhitneyUTest(
                fixed_gen_gen_dist_mean,
                fixed_gen_bif_dist_mean)
        pval = HypothesisTests.pvalue(mwu_test, tail = :both)
        write(stdout, "WSR test, fixed gen, all sites:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(fixed_gen_gen_dist_mean - fixed_gen_bif_dist_mean))\n")
        write(stdout, "\n")
        mwu_test = HypothesisTests.MannWhitneyUTest(
                vo_fixed_gen_gen_dist_mean,
                vo_fixed_gen_bif_dist_mean)
        pval = HypothesisTests.pvalue(mwu_test, tail = :both)
        write(stdout, "WSR test, fixed gen, var only:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(vo_fixed_gen_gen_dist_mean - vo_fixed_gen_bif_dist_mean))\n")
        write(stdout, "\n")
        mwu_test = HypothesisTests.MannWhitneyUTest(
                fixed_bif_gen_dist_mean,
                fixed_bif_bif_dist_mean)
        pval = HypothesisTests.pvalue(mwu_test, tail = :both)
        write(stdout, "WSR test, fixed bif, all sites:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(fixed_bif_gen_dist_mean - fixed_bif_bif_dist_mean))\n")
        write(stdout, "\n")
        mwu_test = HypothesisTests.MannWhitneyUTest(
                vo_fixed_bif_gen_dist_mean,
                vo_fixed_bif_bif_dist_mean)
        pval = HypothesisTests.pvalue(mwu_test, tail = :both)
        write(stdout, "WSR test, fixed bif, var only:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(vo_fixed_bif_gen_dist_mean - vo_fixed_bif_bif_dist_mean))\n")
        write(stdout, "\n")

        mwu_test = HypothesisTests.MannWhitneyUTest(
                unfixed_gen_gen_dist_mean,
                unfixed_gen_bif_dist_mean)
        pval = HypothesisTests.pvalue(mwu_test, tail = :both)
        write(stdout, "WSR test, unfixed gen, all sites:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(unfixed_gen_gen_dist_mean - unfixed_gen_bif_dist_mean))\n")
        write(stdout, "\n")
        mwu_test = HypothesisTests.MannWhitneyUTest(
                vo_unfixed_gen_gen_dist_mean,
                vo_unfixed_gen_bif_dist_mean)
        pval = HypothesisTests.pvalue(mwu_test, tail = :both)
        write(stdout, "WSR test, unfixed gen, var only:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(vo_unfixed_gen_gen_dist_mean - vo_unfixed_gen_bif_dist_mean))\n")
        write(stdout, "\n")
        mwu_test = HypothesisTests.MannWhitneyUTest(
                unfixed_bif_gen_dist_mean,
                unfixed_bif_bif_dist_mean)
        pval = HypothesisTests.pvalue(mwu_test, tail = :both)
        write(stdout, "WSR test, unfixed bif, all sites:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(unfixed_bif_gen_dist_mean - unfixed_bif_bif_dist_mean))\n")
        write(stdout, "\n")
        mwu_test = HypothesisTests.MannWhitneyUTest(
                vo_unfixed_bif_gen_dist_mean,
                vo_unfixed_bif_bif_dist_mean)
        pval = HypothesisTests.pvalue(mwu_test, tail = :both)
        write(stdout, "WSR test, unfixed bif, var only:\n$pval\n")
        write(stdout, "mean diff = $(Statistics.mean(vo_unfixed_bif_gen_dist_mean - vo_unfixed_bif_bif_dist_mean))\n")
        write(stdout, "\n")

        comparison_positions = (0.995, 0.995, 0.86, 0.02)
        if length(locus_prefix) > 0
            comparison_positions = (0.998, 0.998, 0.80, 0.9)
        end
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
                #= [ LaTeXString("\\begin{tabular}{c} Generalized model \\\\ (true model) \\\\ All sites \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model \\\\ All sites \\end{tabular}") LaTeXString("\\begin{tabular}{c} Generalized model \\\\ (true model) \\\\ Only SNPs \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model \\\\ Only SNPs \\end{tabular}") ], =#
                [ gen_model bif_model gen_model bif_model ],
                [ gen_col bif_col vo_gen_col vo_bif_col ],
                [ gen_marker_alpha bif_marker_alpha vo_gen_marker_alpha vo_bif_marker_alpha ],
                marker_shapes = [ gen_shape bif_shape vo_gen_shape vo_bif_shape ],
                marker_sizes = [ gen_marker_size bif_marker_size vo_gen_marker_size vo_bif_marker_size ],
                y_buffer = 0.02,
                x_label_size = 10.0,
                x_label_offset = 0.05,
                show_labels_on_x = true,
                comparisons = ((1, 2), (3, 4), (1, 3), (2, 4)),
                comparison_positions = comparison_positions,
                comparisons_are_paired = true,
               )
        Plots.ylabel!(p, "Euclidean distance from true tree")
        Plots.plot!(p, size = (375, 290))
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-euclidean-distances.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(p, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        if length(locus_prefix) == 0
            l_100_unfixed_gen_gen_dist_mean = get_floats(results,
                    false, true, true, false,
                    :euclidean_distance_mean, 100)
            l_100_unfixed_gen_gen_dist_lower = get_floats(results,
                    false, true, true, false,
                    :euclidean_distance_eti_95_lower, 100)
            l_100_unfixed_gen_gen_dist_upper = get_floats(results,
                    false, true, true, false,
                    :euclidean_distance_eti_95_upper, 100)
            p = get_groups_by_y(
                    hcat(unfixed_gen_gen_dist_mean,
                         l_100_unfixed_gen_gen_dist_mean),
                    hcat(unfixed_gen_gen_dist_lower,
                         l_100_unfixed_gen_gen_dist_lower),
                    hcat(unfixed_gen_gen_dist_upper,
                         l_100_unfixed_gen_gen_dist_upper),
                    #= [ LaTeXString("\\begin{tabular}{c} $(gen_model) \\\\ unlinked \\end{tabular}") LaTeXString("\\begin{tabular}{c} $(gen_model) \\\\ 100bp loci \\end{tabular}") ], =#
                    [ LaTeXString("{\\normalsize $(gen_model) unlinked}") LaTeXString("{\\normalsize $(gen_model) 100bp loci}") ],
                    [ gen_col gen_col ],
                    [ gen_marker_alpha gen_marker_alpha ],
                    marker_shapes = [ gen_shape :pentagon ],
                    marker_sizes = [ gen_marker_size gen_marker_size ],
                    y_buffer = 0.02,
                    x_label_size = 8.0,
                    x_label_offset = 0.05,
                    show_labels_on_x = true,
                    comparisons = ((1, 2),),
                    comparison_positions = (0.995,),
                    comparisons_are_paired = false,
                   )
            Plots.ylabel!(p, "Euclidean distance from true tree")
            Plots.plot!(p, size = (300, 290), xtickfontsize = 8)
            plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-euclidean-distances-unlinked-v-linked.tex")
            write(stdout, "Writing to $(plot_path)\n")
            Plots.savefig(p, plot_path)
            process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
        end

        p = get_groups_by_y(
                hcat(fixed_gen_gen_dist_mean,
                     fixed_gen_bif_dist_mean),
                hcat(fixed_gen_gen_dist_lower,
                     fixed_gen_bif_dist_lower),
                hcat(fixed_gen_gen_dist_upper,
                     fixed_gen_bif_dist_upper),
                [ LaTeXString("\\begin{tabular}{c} Generalized model \\\\ (true model) \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model \\end{tabular}") ],
                [ gen_col bif_col ],
                [ gen_marker_alpha bif_marker_alpha ],
                marker_shapes = [ gen_shape bif_shape ],
                marker_sizes = [ gen_marker_size bif_marker_size ],
                y_buffer = 0.02,
                show_labels_on_x = true,
                comparisons = ((1, 2),),
                comparison_positions = (0.995,),
                comparisons_are_paired = true,
               )
        Plots.ylabel!(p, "Euclidean distance from true tree")
        Plots.plot!(p, size = (320, 290), xtickfontsize = 8)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-all-sites-euclidean-distances.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(p, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        Plots.plot!(p, legend = :top)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-all-sites-euclidean-distances-with-legend.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(p, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        comparison_positions = (0.995, 0.995, 0.09, 0.03)
        if length(locus_prefix) > 0
            comparison_positions = (0.998, 0.998, 0.83, 0.92)
        end
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
                #= [ LaTeXString("\\begin{tabular}{c} Generalized model \\\\ All sites \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model (true model) \\\\ All sites \\end{tabular}") LaTeXString("\\begin{tabular}{c} Generalized model \\\\ Only SNPs \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model (true model) \\\\ Only SNPs \\end{tabular}") ], =#
                [ gen_model bif_model gen_model bif_model ],
                [ gen_col bif_col vo_gen_col vo_bif_col ],
                [ gen_marker_alpha bif_marker_alpha vo_gen_marker_alpha vo_bif_marker_alpha ],
                marker_shapes = [ gen_shape bif_shape vo_gen_shape vo_bif_shape ],
                marker_sizes = [ gen_marker_size bif_marker_size vo_gen_marker_size vo_bif_marker_size ],
                y_buffer = 0.02,
                x_label_size = 10.0,
                x_label_offset = 0.05,
                show_labels_on_x = true,
                comparisons = ((1, 2), (3, 4), (1, 3), (2, 4)),
                comparison_positions = comparison_positions,
                comparisons_are_paired = true,
               )
        Plots.ylabel!(p, "Euclidean distance from true tree")
        #= Plots.plot!(p, size = (600, 350)) =#
        Plots.plot!(p, size = (375, 290))
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-euclidean-distances.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(p, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        Plots.plot!(p, size = (600, 290))
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-euclidean-distances-large.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(p, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        p = get_groups_by_y(
                hcat(fixed_bif_gen_dist_mean,
                     fixed_bif_bif_dist_mean),
                hcat(fixed_bif_gen_dist_lower,
                     fixed_bif_bif_dist_lower),
                hcat(fixed_bif_gen_dist_upper,
                     fixed_bif_bif_dist_upper),
                [ LaTeXString("\\begin{tabular}{c} Generalized \\\\ model \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model (true model) \\end{tabular}") ],
                [ gen_col bif_col ],
                [ gen_marker_alpha bif_marker_alpha ],
                marker_shapes = [ gen_shape bif_shape ],
                marker_sizes = [ gen_marker_size bif_marker_size ],
                y_buffer = 0.02,
                show_labels_on_x = true,
                comparisons = ((1, 2),),
                comparison_positions = (0.995,),
                comparisons_are_paired = true,
               )
        Plots.ylabel!(p, "Euclidean distance from true tree")
        Plots.plot!(p, size = (320, 290), xtickfontsize = 8)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-bif-all-sites-euclidean-distances.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(p, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)


        y_axis_max = -1.0
        if length(locus_prefix) == 0
            y_axis_max = 0.03
        end
        comparison_positions = (0.995, 0.995)
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
                #= [ LaTeXString("\\begin{tabular}{c} Generalized model \\\\ (true model) \\\\ All sites \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model \\\\ All sites \\end{tabular}") LaTeXString("\\begin{tabular}{c} Generalized model \\\\ (true model) \\\\ Only SNPs \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model \\\\ Only SNPs \\end{tabular}") ], =#
                [ gen_model bif_model gen_model bif_model ],
                [ gen_col bif_col vo_gen_col vo_bif_col ],
                [ gen_marker_alpha bif_marker_alpha vo_gen_marker_alpha vo_bif_marker_alpha ],
                marker_shapes = [ gen_shape bif_shape vo_gen_shape vo_bif_shape ],
                marker_sizes = [ gen_marker_size bif_marker_size vo_gen_marker_size vo_bif_marker_size ],
                y_buffer = 0.02,
                y_axis_max = y_axis_max,
                x_label_size = 10.0,
                x_label_offset = 0.05,
                show_labels_on_x = true,
                comparisons = ((1, 2), (3, 4)),
                comparison_positions = comparison_positions,
                comparisons_are_paired = true,
               )
        Plots.ylabel!(p, "Euclidean distance from true tree")
        #= Plots.plot!(p, size = (600, 340)) =#
        Plots.plot!(p, size = (420, 290))
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-euclidean-distances.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(p, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        p = get_groups_by_y(
                hcat(unfixed_gen_gen_dist_mean,
                     unfixed_gen_bif_dist_mean),
                hcat(unfixed_gen_gen_dist_lower,
                     unfixed_gen_bif_dist_lower),
                hcat(unfixed_gen_gen_dist_upper,
                     unfixed_gen_bif_dist_upper),
                [ LaTeXString("\\begin{tabular}{c} Generalized model \\\\ (true model) \\\\ All sites \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model \\\\ All sites \\end{tabular}") ],
                [ gen_col bif_col ],
                [ gen_marker_alpha bif_marker_alpha ],
                marker_shapes = [ gen_shape bif_shape ],
                marker_sizes = [ gen_marker_size bif_marker_size ],
                y_buffer = 0.02,
                show_labels_on_x = true)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-all-sites-euclidean-distances.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(p, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        p = get_groups_by_y(
                hcat(vo_unfixed_gen_gen_dist_mean,
                     vo_unfixed_gen_bif_dist_mean),
                hcat(vo_unfixed_gen_gen_dist_lower,
                     vo_unfixed_gen_bif_dist_lower),
                hcat(vo_unfixed_gen_gen_dist_upper,
                     vo_unfixed_gen_bif_dist_upper),
                [ LaTeXString("\\begin{tabular}{c} Generalized model \\\\ (true model) \\\\ Only SNPs \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model \\\\ Only SNPs \\end{tabular}") ],
                [ vo_gen_col vo_bif_col ],
                [ vo_gen_marker_alpha vo_bif_marker_alpha ],
                y_buffer = 0.02,
                show_labels_on_x = true)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-var-only-euclidean-distances.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(p, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        comparison_positions = (0.995, 0.995)
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
                #= [ LaTeXString("\\begin{tabular}{c} Generalized model \\\\ (true model) \\\\ All sites \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model \\\\ All sites \\end{tabular}") LaTeXString("\\begin{tabular}{c} Generalized model \\\\ (true model) \\\\ Only SNPs \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model \\\\ Only SNPs \\end{tabular}") ], =#
                [ gen_model bif_model gen_model bif_model ],
                [ gen_col bif_col vo_gen_col vo_bif_col ],
                [ gen_marker_alpha bif_marker_alpha vo_gen_marker_alpha vo_bif_marker_alpha ],
                marker_shapes = [ gen_shape bif_shape vo_gen_shape vo_bif_shape ],
                marker_sizes = [ gen_marker_size bif_marker_size vo_gen_marker_size vo_bif_marker_size ],
                y_buffer = 0.02,
                x_label_size = 10.0,
                x_label_offset = 0.05,
                show_labels_on_x = true,
                comparisons = ((1, 2), (3, 4)),
                comparison_positions = comparison_positions
               )
        Plots.ylabel!(p, "Euclidean distance from true tree")
        #= Plots.plot!(p, size = (600, 340)) =#
        Plots.plot!(p, size = (420, 290))
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-bif-euclidean-distances.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(p, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        p = get_groups_by_y(
                hcat(unfixed_bif_gen_dist_mean,
                     unfixed_bif_bif_dist_mean),
                hcat(unfixed_bif_gen_dist_lower,
                     unfixed_bif_bif_dist_lower),
                hcat(unfixed_bif_gen_dist_upper,
                     unfixed_bif_bif_dist_upper),
                [ LaTeXString("\\begin{tabular}{c} Generalized model \\\\ All sites \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model (true model) \\\\ All sites \\end{tabular}") ],
                [ gen_col bif_col ],
                [ gen_marker_alpha bif_marker_alpha ],
                marker_shapes = [ gen_shape bif_shape ],
                marker_sizes = [ gen_marker_size bif_marker_size ],
                y_buffer = 0.02,
                show_labels_on_x = true)
        Plots.ylabel!(p, "Euclidean distance from true tree")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-bif-all-sites-euclidean-distances.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(p, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        p = get_groups_by_y(
                hcat(vo_unfixed_bif_gen_dist_mean,
                     vo_unfixed_bif_bif_dist_mean),
                hcat(vo_unfixed_bif_gen_dist_lower,
                     vo_unfixed_bif_bif_dist_lower),
                hcat(vo_unfixed_bif_gen_dist_upper,
                     vo_unfixed_bif_bif_dist_upper),
                [ LaTeXString("\\begin{tabular}{c} Generalized model \\\\ Only SNPs \\end{tabular}") LaTeXString("\\begin{tabular}{c} Independent bifurcating \\\\ model (true model) \\\\ Only SNPs \\end{tabular}") ],
                [ vo_gen_col vo_bif_col ],
                [ vo_gen_marker_alpha vo_bif_marker_alpha ],
                y_buffer = 0.02,
                show_labels_on_x = true)
        Plots.ylabel!(p, "Euclidean distance from true tree")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-bif-var-only-euclidean-distances.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(p, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)


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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed_gen_asdsf.tex")
        Plots.savefig(vln_fixed_gen_asdsf, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_split_violin_plot(
                fixed_gen_gen_asdsf,
                fixed_gen_bif_asdsf,
                xlabels = [""],
                left_fill_colors = [gen_col],
                left_marker_colors = [gen_col],
                left_fill_alphas = [gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha],
                left_marker_shapes = [gen_shape],
                left_marker_sizes = [gen_marker_size],
                left_labels = "Generalized",
                right_fill_colors = [bif_col],
                right_marker_colors = [bif_col],
                right_fill_alphas = [bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha],
                right_marker_shapes = [bif_shape],
                right_marker_sizes = [bif_marker_size],
                right_labels = "Bifurcating",
                legend = :top,
                dot_legend = false)
        Plots.plot!(v, size = (220, 280))
        Plots.ylabel!(v, "Ave std dev of split freq (ASDSF)")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-asdsf-all-sites.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_split_dot_plot(
                fixed_gen_gen_asdsf,
                fixed_gen_bif_asdsf,
                xlabels = [""],
                left_fill_colors = [gen_col],
                left_marker_colors = [gen_col],
                left_fill_alphas = [gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha],
                left_labels = "Generalized",
                right_fill_colors = [bif_col],
                right_marker_colors = [bif_col],
                right_fill_alphas = [bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha],
                right_labels = "Bifurcating",
                legend = :top)
        Plots.plot!(v, size = (220, 280))
        Plots.ylabel!(v, "Ave std dev of split freq (ASDSF)")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-gen-asdsf-all-sites-dot.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        vln_fixed_bif_asdsf = get_split_violin_plot(
                hcat(fixed_bif_gen_asdsf, vo_fixed_bif_gen_asdsf),
                hcat(fixed_bif_bif_asdsf, vo_fixed_bif_bif_asdsf),
                xlabels = ["All sites" "Variable sites"],
                left_fill_colors = [gen_col vo_gen_col],
                left_marker_colors = [gen_col vo_gen_col],
                left_fill_alphas = [gen_fill_alpha vo_gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha vo_gen_marker_alpha],
                left_marker_shapes = [ gen_shape vo_gen_shape ],
                left_marker_sizes = [ gen_marker_size vo_gen_marker_size],
                left_labels = "Generalized",
                right_fill_colors = [bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha],
                right_marker_shapes = [ bif_shape vo_bif_shape ],
                right_marker_sizes = [ bif_marker_size vo_bif_marker_size],
                right_labels = "Bifurcating",
                legend = false,
                dot_legend = false)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed_bif_asdsf.tex")
        Plots.savefig(vln_fixed_bif_asdsf, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)fixed-asdsf.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(vln_fixed_asdsf_grid, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)


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
                left_marker_shapes = [ gen_shape vo_gen_shape ],
                left_marker_sizes = [ gen_marker_size vo_gen_marker_size],
                left_labels = "Generalized",
                right_fill_colors = [bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha],
                right_marker_shapes = [ bif_shape vo_bif_shape ],
                right_marker_sizes = [ bif_marker_size vo_bif_marker_size],
                right_labels = "Bifurcating",
                legend = false,
                dot_legend = false)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed_gen_asdsf.tex")
        Plots.savefig(vln_unfixed_gen_asdsf, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        vln_unfixed_bif_asdsf = get_split_violin_plot(
                hcat(unfixed_bif_gen_asdsf, vo_unfixed_bif_gen_asdsf),
                hcat(unfixed_bif_bif_asdsf, vo_unfixed_bif_bif_asdsf),
                xlabels = ["All sites" "Variable sites"],
                left_fill_colors = [gen_col vo_gen_col],
                left_marker_colors = [gen_col vo_gen_col],
                left_fill_alphas = [gen_fill_alpha vo_gen_fill_alpha],
                left_marker_alphas = [gen_marker_alpha vo_gen_marker_alpha],
                left_marker_shapes = [ gen_shape vo_gen_shape ],
                left_marker_sizes = [ gen_marker_size vo_gen_marker_size],
                left_labels = "Generalized",
                right_fill_colors = [bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha],
                right_marker_shapes = [ bif_shape vo_bif_shape ],
                right_marker_sizes = [ bif_marker_size vo_bif_marker_size],
                right_labels = "Bifurcating",
                legend = false,
                dot_legend = false)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed_bif_asdsf.tex")
        Plots.savefig(vln_unfixed_bif_asdsf, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-asdsf.tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(vln_unfixed_asdsf_grid, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)


        # Plot topology-related probs

        unfixed_gen_gen_true_topo_probs = get_floats(results, false, true, true, false, :topo_true_prob, locus_size)
        unfixed_gen_gen_vo_true_topo_probs = get_floats(results, false, true, true, true, :topo_true_prob, locus_size)
        unfixed_gen_gen_root_node_probs = get_floats(results, false, true, true, false, :root_node_true_prob, locus_size)
        unfixed_gen_gen_vo_root_node_probs = get_floats(results, false, true, true, true, :root_node_true_prob, locus_size)
        unfixed_gen_gen_true_split_prob_means = get_floats(results, false, true, true, false, :true_split_prob_mean, locus_size)
        unfixed_gen_gen_vo_true_split_prob_means = get_floats(results, false, true, true, true, :true_split_prob_mean, locus_size)
        unfixed_gen_gen_true_node_prob_means = get_floats(results, false, true, true, false, :true_node_prob_mean, locus_size)
        unfixed_gen_gen_vo_true_node_prob_means = get_floats(results, false, true, true, true, :true_node_prob_mean, locus_size)
        unfixed_gen_gen_true_height_prob_means = get_floats(results, false, true, true, false, :true_height_prob_mean, locus_size)
        unfixed_gen_gen_vo_true_height_prob_means = get_floats(results, false, true, true, true, :true_height_prob_mean, locus_size)

        v = get_floats(results, false, true, true, false, :true_shared_height_prob_mean, locus_size)
        true_shared_height_prob_means = v[.!any.(isnan, v)]
        v = get_floats(results, false, true, true, true, :true_shared_height_prob_mean, locus_size)
        vo_true_shared_height_prob_means = v[.!any.(isnan, v)]

        v = get_split_violin_plot(
                true_shared_height_prob_means,
                vo_true_shared_height_prob_means,
                xlabels = [ "Shared divergence" ],
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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-gen-mean-true-shared-height-probs.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_violin_plot(
                true_shared_height_prob_means,
                xlabels = [ "Shared divergence" ],
                fill_colors = gen_col,
                marker_colors = gen_col,
                fill_alphas = gen_fill_alpha,
                marker_alphas = gen_marker_alpha,
                include_dots = true)
        Plots.plot!(v, size = (200, 200))
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-gen-mean-true-shared-height-probs-all-sites.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_split_violin_plot(
                hcat(unfixed_gen_gen_true_topo_probs,
                     unfixed_gen_gen_root_node_probs,
                     unfixed_gen_gen_true_split_prob_means,
                     unfixed_gen_gen_true_node_prob_means,
                     unfixed_gen_gen_true_height_prob_means),
                hcat(unfixed_gen_gen_vo_true_topo_probs,
                     unfixed_gen_gen_vo_root_node_probs,
                     unfixed_gen_gen_vo_true_split_prob_means,
                     unfixed_gen_gen_vo_true_node_prob_means,
                     unfixed_gen_gen_vo_true_height_prob_means),
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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-gen-true-tree-probs.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)


        unfixed_bif_bif_true_topo_probs = get_floats(results, false, false, false, false, :topo_true_prob, locus_size)
        unfixed_bif_bif_vo_true_topo_probs = get_floats(results, false, false, false, true, :topo_true_prob, locus_size)
        unfixed_bif_bif_root_node_probs = get_floats(results, false, false, false, false, :root_node_true_prob, locus_size)
        unfixed_bif_bif_vo_root_node_probs = get_floats(results, false, false, false, true, :root_node_true_prob, locus_size)
        unfixed_bif_bif_true_split_prob_means = get_floats(results, false, false, false, false, :true_split_prob_mean, locus_size)
        unfixed_bif_bif_vo_true_split_prob_means = get_floats(results, false, false, false, true, :true_split_prob_mean, locus_size)
        unfixed_bif_bif_true_node_prob_means = get_floats(results, false, false, false, false, :true_node_prob_mean, locus_size)
        unfixed_bif_bif_vo_true_node_prob_means = get_floats(results, false, false, false, true, :true_node_prob_mean, locus_size)
        unfixed_bif_bif_true_height_prob_means = get_floats(results, false, false, false, false, :true_height_prob_mean, locus_size)
        unfixed_bif_bif_vo_true_height_prob_means = get_floats(results, false, false, false, true, :true_height_prob_mean, locus_size)

        v = get_split_violin_plot(
                hcat(unfixed_bif_bif_true_topo_probs,
                     unfixed_bif_bif_root_node_probs,
                     unfixed_bif_bif_true_split_prob_means,
                     unfixed_bif_bif_true_node_prob_means,
                     unfixed_bif_bif_true_height_prob_means),
                hcat(unfixed_bif_bif_vo_true_topo_probs,
                     unfixed_bif_bif_vo_root_node_probs,
                     unfixed_bif_bif_vo_true_split_prob_means,
                     unfixed_bif_bif_vo_true_node_prob_means,
                     unfixed_bif_bif_vo_true_height_prob_means),
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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-bif-bif-true-tree-probs.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)


        unfixed_gen_bif_true_topo_probs = get_floats(results, false, true, false, false, :topo_true_prob, locus_size)
        unfixed_gen_bif_vo_true_topo_probs = get_floats(results, false, true, false, true, :topo_true_prob, locus_size)
        unfixed_gen_bif_root_node_probs = get_floats(results, false, true, false, false, :root_node_true_prob, locus_size)
        unfixed_gen_bif_vo_root_node_probs = get_floats(results, false, true, false, true, :root_node_true_prob, locus_size)
        unfixed_gen_bif_true_split_prob_means = get_floats(results, false, true, false, false, :true_split_prob_mean, locus_size)
        unfixed_gen_bif_vo_true_split_prob_means = get_floats(results, false, true, false, true, :true_split_prob_mean, locus_size)
        unfixed_gen_bif_true_node_prob_means = get_floats(results, false, true, false, false, :true_node_prob_mean, locus_size)
        unfixed_gen_bif_vo_true_node_prob_means = get_floats(results, false, true, false, true, :true_node_prob_mean, locus_size)
        unfixed_gen_bif_true_height_prob_means = get_floats(results, false, true, false, false, :true_height_prob_mean, locus_size)
        unfixed_gen_bif_vo_true_height_prob_means = get_floats(results, false, true, false, true, :true_height_prob_mean, locus_size)

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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-bif-true-shared-height-probs.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_split_violin_plot(
                hcat(unfixed_gen_bif_true_topo_probs,
                     unfixed_gen_bif_root_node_probs,
                     unfixed_gen_bif_true_split_prob_means,
                     unfixed_gen_bif_true_node_prob_means,
                     unfixed_gen_bif_true_height_prob_means),
                hcat(unfixed_gen_bif_vo_true_topo_probs,
                     unfixed_gen_bif_vo_root_node_probs,
                     unfixed_gen_bif_vo_true_split_prob_means,
                     unfixed_gen_bif_vo_true_node_prob_means,
                     unfixed_gen_bif_vo_true_height_prob_means),
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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-gen-bif-true-tree-probs.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)


        unfixed_bif_gen_true_topo_probs = get_floats(results, false, false, true, false, :topo_true_prob, locus_size)
        unfixed_bif_gen_vo_true_topo_probs = get_floats(results, false, false, true, true, :topo_true_prob, locus_size)
        unfixed_bif_gen_root_node_probs = get_floats(results, false, false, true, false, :root_node_true_prob, locus_size)
        unfixed_bif_gen_vo_root_node_probs = get_floats(results, false, false, true, true, :root_node_true_prob, locus_size)
        unfixed_bif_gen_true_split_prob_means = get_floats(results, false, false, true, false, :true_split_prob_mean, locus_size)
        unfixed_bif_gen_vo_true_split_prob_means = get_floats(results, false, false, true, true, :true_split_prob_mean, locus_size)
        unfixed_bif_gen_true_node_prob_means = get_floats(results, false, false, true, false, :true_node_prob_mean, locus_size)
        unfixed_bif_gen_vo_true_node_prob_means = get_floats(results, false, false, true, true, :true_node_prob_mean, locus_size)
        unfixed_bif_gen_true_height_prob_means = get_floats(results, false, false, true, false, :true_height_prob_mean, locus_size)
        unfixed_bif_gen_vo_true_height_prob_means = get_floats(results, false, false, true, true, :true_height_prob_mean, locus_size)

        v = get_split_violin_plot(
                hcat(unfixed_bif_gen_true_topo_probs,
                     unfixed_bif_gen_root_node_probs,
                     unfixed_bif_gen_true_split_prob_means,
                     unfixed_bif_gen_true_node_prob_means,
                     unfixed_bif_gen_true_height_prob_means),
                hcat(unfixed_bif_gen_vo_true_topo_probs,
                     unfixed_bif_gen_vo_root_node_probs,
                     unfixed_bif_gen_vo_true_split_prob_means,
                     unfixed_bif_gen_vo_true_node_prob_means,
                     unfixed_bif_gen_vo_true_height_prob_means),
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
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-bif-gen-true-tree-probs.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

        v = get_split_violin_plot(
                hcat(unfixed_bif_gen_true_split_prob_means,
                     unfixed_bif_gen_vo_true_split_prob_means),
                hcat(unfixed_bif_bif_true_split_prob_means,
                     unfixed_bif_bif_vo_true_split_prob_means),
                xlabels = [ "True split mean" " True split mean " ],
                left_fill_colors = [ gen_col  vo_gen_col ],
                left_marker_colors = [ gen_col  vo_gen_col ],
                left_fill_alphas = [ gen_fill_alpha  vo_gen_fill_alpha ],
                left_marker_alphas = [ gen_marker_alpha  vo_gen_marker_alpha ],
                left_marker_shapes = [ gen_shape vo_gen_shape ],
                left_marker_sizes = [ gen_marker_size vo_gen_marker_size ],
                left_labels = [ "Generalized" ],
                right_fill_colors = [bif_col vo_bif_col],
                right_marker_colors = [bif_col vo_bif_col],
                right_fill_alphas = [bif_fill_alpha vo_bif_fill_alpha],
                right_marker_alphas = [bif_marker_alpha vo_bif_marker_alpha],
                right_marker_shapes = [ bif_shape vo_bif_shape ],
                right_marker_sizes = [ bif_marker_size vo_bif_marker_size ],
                right_labels = [ "Bifurcating" ],
                legend = false,
                dot_legend = false)
        Plots.plot!(v, size = (300, 220), xtickfontsize = 8)
        Plots.ylims!(v, (-0.02, 1.02))
        Plots.ylabel!(v, "Posterior probability")
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "$(locus_prefix)unfixed-bif-mean-split-probs.tex")
        Plots.savefig(v, plot_path)
        process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

    end

    parameters_to_plot = ["tree_length", "root_height", "root_pop_size"]
    parameter_symbols = ["\\textrm{\\sffamily TL}", "\\tau_{n(\\tau)}", "N_e"]
    for param_index in 1:length(parameters_to_plot)
        parameter = parameters_to_plot[param_index]
        parameter_symbol = parameter_symbols[param_index]
        true_v_est_plots = []
        coverage_stats = Vector{Float64}()
        rmse_stats = Vector{Float64}()
        plot_labels = Vector{String}()
        for locus_size in locus_sizes
            locus_prefix = ""
            if locus_size > 1
                locus_prefix = "locus-$locus_size-"
            end
            for sim_fixed in [false]
                for sim_generalized in [true, false]
                    for analysis_generalized in [true, false]
                        for only_var_sites in [false, true]
                            plot_title = "$(sim_fixed ? "fixed" : "free")"
                            plot_title *= "$(sim_generalized ? "-gen" : "-bif")"
                            plot_title *= "$(analysis_generalized ? "-gen" : "-bif")"
                            plot_title *= "$(only_var_sites ? "-var" : "-all")"
                            plot_label = "$(locus_prefix)$(plot_title)"
                            scatter_plot = get_true_v_est_plot(
                                    results,
                                    sim_fixed,
                                    sim_generalized,
                                    analysis_generalized,
                                    only_var_sites,
                                    parameter,
                                    locus_size,
                                    false,
                                    0.05,
                                    brooks_gelman_1998_recommended_psrf,
                                    ess_threshold,
                                    "red")
                            # Plots.title!(scatter_plot, plot_title) 
                            #= Plots.plot!(scatter_plot, margin = 0px) =#
                            coverage = get_floats(
                                    sum_results,
                                    sim_fixed,
                                    sim_generalized,
                                    analysis_generalized,
                                    only_var_sites,
                                    "coverage_$(parameter)",
                                    locus_size)[1]
                            rmse = get_floats(
                                    sum_results,
                                    sim_fixed,
                                    sim_generalized,
                                    analysis_generalized,
                                    only_var_sites,
                                    "rmse_$(parameter)",
                                    locus_size)[1]
                            push!(true_v_est_plots, scatter_plot)
                            push!(coverage_stats, coverage)
                            push!(rmse_stats, rmse)
                            push!(plot_labels, plot_label)
                        end
                    end
                end
            end
        end

        xy_lims = share_xy_limits!(true_v_est_plots...)
        add_identity_line!(true_v_est_plots...)
        coverage_position = relative_xy(true_v_est_plots[1], 0.03, 0.98)
        rmse_position = relative_xy(true_v_est_plots[1], 0.03, 0.91)

        for i in 1:length(true_v_est_plots)
            if xy_lims[2] < 0.01
                Plots.plot!(true_v_est_plots[i], ticks = optimize_ticks(xy_lims..., k_min = 2, k_max = 4)[1])
            end
            ci_str = "\\textrm{\\sffamily CI}"
            rmse_str = "\\textrm{\\sffamily RMSE}"
            cov_str = L"$p(%$(parameter_symbol) \in %$(ci_str)) = %$(coverage_stats[i])$"
            rmse_val_str = @sprintf("%5.3g", rmse_stats[i])
            m = match(ProjectUtil.SCI_NOTATION_PATTERN, rmse_val_str)
            if ! isnothing(m)
                e_sign = ""
                if m[:sign] == "-"
                    e_sign = "-"
                end
                rmse_val_str = "$(m[:digits]) \\times 10^{$(e_sign)$(m[:exponent])}"
            end
            rmse_str = L"$ %$(rmse_str) = %$(rmse_val_str)$"
            Plots.plot!(true_v_est_plots[i], size = (290, 250))
            #= Plots.plot!(true_v_est_plots[i], =#
            #=         bottom_margin = -2.5mm, =#
            #=         left_margin = -1mm, =#
            #=         top_margin = -1mm, =#
            #=         right_marin = -2mm) =#
            annotate!(true_v_est_plots[i],
                    coverage_position...,
                    #= text(L"$p(\textrm{TL} \in \textrm{CI}) = %$tree_len_coverage$", =#
                    #= text("\$\\textrm{\\sffamily\\it p}(\\textrm{\\sffamily TL} \\in \\textrm{\\sffamily CI}) = $(tree_len_coverage)\$", =#
                    #= text(L"$\textrm{\sffamily\it p}(\textrm{\sffamily TL} \in \textrm{\sffamily CI}) = %$(tree_len_coverage)$", =#
                    #= text("p(TL in CI) = $(tree_len_coverage[i])", =#
                    #= text(L"p(\textrm{TL} \in \textrm{CI}) = %$(tree_len_coverage[i])", =#
                    #= text(L"p(\textrm{TL} \in \textrm{CI}) = %$(tree_len_coverage[i])", =#
                    text(cov_str,
                            :left,
                            :top,
                            8),
                    annotation_clip = false)
            annotate!(true_v_est_plots[i],
                    rmse_position...,
                    text(rmse_str,
                            :left,
                            :top,
                            8),
                    annotation_clip = false)
            plot_path = joinpath(ProjectUtil.RESULTS_DIR,
                    "$(plot_labels[i])-$(parameter).tex")
            Plots.savefig(true_v_est_plots[i], plot_path)
            process_tex(plot_path, target = axis_pattern, replacement = axis_replace)

            x_label = "True $(replace(parameter, "_" => " "))"
            y_label = "Estimated $(replace(parameter, "_" => " "))"
            if parameter == "root_pop_size"
                x_label = "True population size"
                y_label = "Estimated population size"
            end
            Plots.plot!(true_v_est_plots[i], xlabel = x_label, ylabel = y_label)
            plot_path = joinpath(ProjectUtil.RESULTS_DIR,
                    "$(plot_labels[i])-$(parameter)-axis-labels.tex")
            Plots.savefig(true_v_est_plots[i], plot_path)
            process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
            Plots.plot!(true_v_est_plots[i], xlabel = "", ylabel = "")
        end

        plot_grid = Plots.plot(
                true_v_est_plots...,
                layout = (4, 4),
                legend = false,
                size = (700, 1200),
                link = :all,
                xlabel = "True $(parameter)",
                ylabel = "Posterior mean $(parameter)",
                )
        share_xy_axes!(plot_grid)
        plot_path = joinpath(ProjectUtil.RESULTS_DIR, "unfixed-$(parameter).tex")
        write(stdout, "Writing to $(plot_path)\n")
        Plots.savefig(plot_grid, plot_path)
        # process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
    end

    return 0
end

main_cli()
