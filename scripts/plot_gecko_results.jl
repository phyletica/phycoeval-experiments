#! /usr/bin/env julia

using YAML
using Plots

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


function main_cli()::Cint
    write(stdout, "Plotting backend: $(backend())\n")

    taxa = ["gekko", "cyrt"]
    for taxon in taxa
        yml_path = joinpath(ProjectUtil.GEK_OUT_DIR, "posterior-summary-$(taxon)-nopoly.yml")
        d = YAML.load_file(yml_path)
        max_num_divs = length(d["leaf_label_map"]) - 1
        num_div_freqs = zeros(Float64, max_num_divs)
        for num_height_entry in d["numbers_of_heights"]
            num_heights::Int64 = num_height_entry["number_of_heights"]
            freq::Float64 = num_height_entry["frequency"]
            num_div_freqs[num_heights] = freq
        end

        p = Plots.bar(num_div_freqs,
                      fillcolor = gen_col,
                      linecolor = gen_col,
                      fillalpha = 1.0,
                      linealpha = 0.0,
                      xlabel = "Number of divergences",
                      ylabel = "Posterior probability",
                      #= yticks = (1:max_num_divs, string.(1:max_num_divs)), =#
                      legend = false)
        Plots.plot!(p, size = (400, 280))
        plot_path = joinpath(ProjectUtil.GEK_OUT_DIR, "number-of-divs-$(taxon).tex")
        Plots.savefig(p, plot_path)
        pdf_path = process_tex(plot_path, target = axis_pattern, replacement = axis_replace)
    end
    return 0
end

main_cli()
