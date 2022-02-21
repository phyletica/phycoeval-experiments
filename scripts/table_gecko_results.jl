#! /usr/bin/env julia

using Printf

include("project_util.jl")


struct SharedDiv
    index::Int
    freq::Float64
    mean::Float64
    median::Float64
    eti_95_lower::Float64
    eti_95_upper::Float64
    hpdi_95_lower::Float64
    hpdi_95_upper::Float64
    num_nodes::Int
end

function parse_shared_divergences(newick_str)
    height_indices = Set()
    shared_divs = Dict{Int, RegexMatch}()
    shared_div_counts = Dict{Int, Int}()
    for m in eachmatch(ProjectUtil.INTERNAL_NODE_ANNOTATION_PATTERN, newick_str)
        height_index = parse(Int, m[:height_index])
        if height_index in height_indices
            if height_index in keys(shared_divs)
                shared_div_counts[height_index] += 1
                @assert m[:index_height_mean] == shared_divs[height_index][:index_height_mean]
            else
                shared_divs[height_index] = m
                shared_div_counts[height_index] = 2
            end
        end
        push!(height_indices, height_index)
    end
    @assert keys(shared_divs) == keys(shared_div_counts)
    max_index = maximum(height_indices)
    for i in 1:max_index
        @assert i in height_indices
    end
    shared_divergences = []
    for height_index in sort([k for k in keys(shared_divs)], rev = true)
        m = shared_divs[height_index]
        sd = SharedDiv(parse(Int, m[:height_index]),
                       parse(Float64, m[:index_freq]),
                       parse(Float64, m[:index_height_mean]),
                       parse(Float64, m[:index_height_median]),
                       parse(Float64, m[:index_height_eti_95_lower]),
                       parse(Float64, m[:index_height_eti_95_upper]),
                       parse(Float64, m[:index_height_hpdi_95_lower]),
                       parse(Float64, m[:index_height_hpdi_95_upper]),
                       shared_div_counts[height_index])
        push!(shared_divergences, sd)
    end
    return shared_divergences
end

function main_cli()::Cint
    sep = "\t"
    hline = ""
    line_end = "\n"
    for a in ARGS
        if a == "--tex"
            sep = " & "
            hline = "\\hline\n"
            line_end = " \\\\\n"
            break
        end
    end
    nex_tree_paths = [
        joinpath(ProjectUtil.GEK_OUT_DIR, "scaled-map-tree-cyrt-nopoly.nex"),
        joinpath(ProjectUtil.GEK_OUT_DIR, "scaled-map-tree-gekko-nopoly.nex"),
       ]
    genus_labels = ["C", "G"]

    write(stdout, "Shared divergence label$(sep)Number of nodes$(sep)Posterior prob.$(sep)Mean time$(sep)Time 95% HPDI$(line_end)")
    for genus_index in 1:length(nex_tree_paths)
        write(stdout, "$hline")
        tree_path = nex_tree_paths[genus_index]
        genus_label = genus_labels[genus_index]
        open(tree_path, "r") do istream
            for line in eachline(istream)
                m = match(ProjectUtil.NEXUS_TREE_LINE_PATTERN, line)
                if ! isnothing(m)
                    newick_str = m[:newick_str]
                    shared_divs = parse_shared_divergences(newick_str)
                    for i in 1:length(shared_divs)
                        sd = shared_divs[i]
                        label = "$(genus_labels[genus_index])$i"
                        pp = @sprintf("%.3f", sd.freq)
                        mean = @sprintf("%.3f", sd.mean)
                        hpdi = @sprintf("%.3f--%.3f", sd.hpdi_95_lower, sd.hpdi_95_upper)
                        write(stdout, "$label$(sep)$(sd.num_nodes)$(sep)$pp$(sep)$mean$(sep)$hpdi$line_end")
                    end
                end
            end
        end
    end
    write(stdout, "$hline")
    return 0
end


main_cli()
