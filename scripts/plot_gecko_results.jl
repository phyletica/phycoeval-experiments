#! /usr/bin/env julia

using Plots
using StatsPlots

pgfplotsx()
default(tex_output_standalone = true)

include("project_util.jl")


sans_tex_lines = [
        "\\usepackage[cm]{sfmath}",
        "\\usepackage[T1]{fontenc}",
        "\\renewcommand{\\familydefault}{\\sfdefault}",
        ]
sans_tex_insert = r"^\\usepackage\{pgfplots\}$"

dark_blue_col = RGB(2/255, 75/255, 120/255)

gen_col = dark_blue_col

axis_pattern = r"\\begin\{axis\}\["
axis_replace = "\\begin{axis}[clip=false, "

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

function parse_number_of_divs(
        state_log_paths::String...;
        skip::Int = 0)::Vector{Int}
    df = ProjectUtil.get_data_frame_from_gz(collect(state_log_paths);
            skip = skip)
    return df[:, :number_of_heights]
end

function parse_number_of_tips(
        state_log_path::String)::Int
    df = ProjectUtil.get_data_frame_from_gz([state_log_path])
    num_tips = 0
    for col_name in names(df)
        if startswith(col_name, "pop_size_") & (col_name != "pop_size_root")
            num_tips += 1
        end
    end
    return num_tips
end

function get_number_of_div_probs(
        state_log_paths::String...;
        skip::Int = 0)::Vector{Float64}
    num_divs = Vector{Int}()
    num_paths = 0
    num_tips = -1
    for p in state_log_paths
        write(stdout, "Parsing numbers of divs from $(p)...\n")
        ndivs = parse_number_of_divs(p; skip = skip)
        write(stdout, "Parsed $(length(ndivs)) samples after ignoring first $(skip).\n")
        append!(num_divs, ndivs)
        num_paths += 1
        if num_tips < 0
            num_tips = parse_number_of_tips(p)
        end
    end

    nsamples = length(num_divs)
    write(stdout, "\nParsed a total of $(nsamples) samples from $(num_paths) files.\n")
    write(stdout, "Number of tips: $(num_tips)\n")

    @assert num_tips > 0
    max_num_divs = num_tips - 1
    num_div_counts = zeros(Int, max_num_divs)
    for n in num_divs
        num_div_counts[n] += 1
    end
    num_div_freqs = zeros(Float64, length(num_div_counts))
    for i in 1:length(num_div_counts)
        num_div_freqs[i] = num_div_counts[i] / nsamples
    end

    return num_div_freqs
end


function main_cli()::Cint
    write(stdout, "Plotting backend: $(backend())\n")

    taxa = ["gekko", "cyrt"]
    burnin = 101
    for taxon in taxa
        state_log_path_iter = ProjectUtil.flat_file_path_iter(
                ProjectUtil.GEK_OUT_DIR,
                r"^run-\d+-threads-\d+-" * taxon * r"-nopoly-state-run-1.log.gz$")
        num_div_freqs = get_number_of_div_probs(state_log_path_iter...;
                skip = burnin)
        max_num_divs = length(num_div_freqs)

        have_prior_samples = false
        prior_state_log_path_iter = ProjectUtil.flat_file_path_iter(
                ProjectUtil.GEK_OUT_DIR,
                r"^nodata-run-\d+-threads-\d+-" * taxon * r"-nopoly-long-state-run-1.log.gz$")
        if isready(prior_state_log_path_iter)
            have_prior_samples = true
            prior_num_div_freqs = get_number_of_div_probs(
                    prior_state_log_path_iter...;
                    skip = burnin)
            @assert length(prior_num_div_freqs) == length(num_div_freqs)

            x_ticks = collect(1:max_num_divs)
            x_tick_labels = string.(1:max_num_divs)
            for i in 1:max_num_divs
                if i % 2 == 0
                    x_tick_labels[i] = ""
                end
            end
            groups = repeat(["Prior", "Posterior"], inner = max_num_divs)

            p = StatsPlots.groupedbar(vcat(x_ticks, x_ticks),
                    hcat(prior_num_div_freqs, num_div_freqs),
                    group = groups,
                    xlabel = "Number of divergences",
                    ylabel = "Posterior probability",
                    xticks = (1:max_num_divs, x_tick_labels),
                    linewidth = 0.0,
                    linealpha = 0.0,
                    legend = :topright
                    #= framestyle = :box =#
                   )
            Plots.plot!(p, size = (600, 280))
            plot_path = joinpath(ProjectUtil.GEK_OUT_DIR, "number-of-divs-$(taxon).tex")
            Plots.savefig(p, plot_path)
            pdf_path = process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
            continue
        end

        x_tick_labels = string.(1:max_num_divs)
        for i in 1:max_num_divs
            if i % 2 == 0
                x_tick_labels[i] = ""
            end
        end
        p = Plots.bar(num_div_freqs,
                      fillcolor = gen_col,
                      linecolor = gen_col,
                      fillalpha = 1.0,
                      linealpha = 0.0,
                      xlabel = "Number of divergences",
                      ylabel = "Posterior probability",
                      #= xticks = (1:max_num_divs, string.(1:max_num_divs)), =#
                      xticks = (1:max_num_divs, x_tick_labels),
                      legend = false)
        Plots.plot!(p, size = (400, 280))
        plot_path = joinpath(ProjectUtil.GEK_OUT_DIR, "number-of-divs-$(taxon).tex")
        Plots.savefig(p, plot_path)
        pdf_path = process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
    end
    return 0
end

main_cli()
