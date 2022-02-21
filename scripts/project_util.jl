#! /usr/bin/env julia

module ProjectUtil

using CSV
using DataFrames
using Statistics
using CodecZlib
using Mmap
using Random
using Printf

# Project paths
SCRIPT_DIR = Base.Filesystem.dirname(@__FILE__)
PROJECT_DIR = Base.Filesystem.dirname(SCRIPT_DIR)
CONFIG_DIR = Base.Filesystem.joinpath(PROJECT_DIR, "configs")
BIN_DIR = Base.Filesystem.joinpath(PROJECT_DIR, "bin")
DATA_DIR = Base.Filesystem.joinpath(PROJECT_DIR, "data")
SIM_DIR = Base.Filesystem.joinpath(PROJECT_DIR, "simulations")
SIM_SCRIPT_DIR = Base.Filesystem.joinpath(SCRIPT_DIR, "simphycoeval-scripts")
RESULTS_DIR = Base.Filesystem.joinpath(PROJECT_DIR, "results")
GEK_OUT_DIR = Base.Filesystem.joinpath(PROJECT_DIR, "gekkonid-output")

function get_project_dir()::AbstractString
    return PROJECT_DIR
end

# Project regular expressions
SIMPHYCOEVAL_CONFIG_NAME_PATTERN_STR = (
        raw"(?<var_only>var-only-)?" *
        raw"(?<config_prefix>\S*simphycoeval)" *
        raw"-sim-(?<sim_num>\d+)" *
        raw"-(?<config_name>\S+)" *
        raw"-config.yml"
       )

SIMPHYCOEVAL_CONFIG_NAME_PATTERN = Regex(join([
        raw"^",
        SIMPHYCOEVAL_CONFIG_NAME_PATTERN_STR,
        raw"$"
       ]))

SCI_NOTATION_STR = (
        raw"(?<digits>[\d.]+)[Ee](?<sign>[+-])?[0]*(?<exponent>\d+)"
        )
SCI_NOTATION_PATTERN = Regex(join([
        raw"^",
        SCI_NOTATION_STR,
        raw"$"
       ]))

SIM_STATE_LOG_PATTERN_STR = (
        raw"run-(?<run_num>\d+)-" *
        raw"(?<var_only>var-only-)?" *
        raw"(?<config_prefix>\S*simphycoeval)" *
        raw"-sim-(?<sim_num>\d+)" *
        raw"-(?<config_name>\S+)" *
        raw"-config-state-run-(?P<dummy_run_num>\d+).log"
       )
SIM_STATE_LOG_PATTERN = Regex(join([
        raw"^",
        SIM_STATE_LOG_PATTERN_STR,
        raw"$"
       ]))

SIM_TREE_LOG_PATTERN_STR = (
        raw"run-(?<run_num>\d+)-" *
        raw"(?<var_only>var-only-)?" *
        raw"(?<config_prefix>\S*simphycoeval)" *
        raw"-sim-(?<sim_num>\d+)" *
        raw"-(?<config_name>\S+)" *
        raw"-config-trees-run-(?P<dummy_run_num>\d+).nex"
       )
SIM_TREE_LOG_PATTERN = Regex(join([
        raw"^",
        SIM_TREE_LOG_PATTERN_STR,
        raw"$"
       ]))

SIM_TRUE_TREE_FILE_PATTERN_STR = (
        raw"(?<config_prefix>\S*simphycoeval)" *
        raw"-sim-(?<sim_num>\d+)" *
        raw"-true-tree.phy"
       )
SIM_TRUE_TREE_FILE_PATTERN = Regex(join([
        raw"^",
        SIM_TRUE_TREE_FILE_PATTERN_STR,
        raw"$"
       ]))

RUNTIME_PATTERN = Regex(raw"^\s*Runtime:\s+(?<runtime>\d+)\s+seconds\.\s*$")

BATCH_DIR_PATTERN_STR = raw"batch-(?<batch_num>\d+)"
BATCH_DIR_PATTERN = Regex(
        raw"^" * BATCH_DIR_PATTERN_STR * raw"$")
BATCH_DIR_ENDING_PATTERN = Regex(
            raw"^.*" * BATCH_DIR_PATTERN_STR * raw"(" * Base.Filesystem.path_separator * raw")?$")

NEXUS_TREE_LINE_STR = raw"\s*TREE\s+(?<tree_label>\S+)\s+=\s+\[&R]\s(?<newick_str>\(.*;)\s*"
NEXUS_TREE_LINE_PATTERN = Regex(join([
        raw"^",
        NEXUS_TREE_LINE_STR,
        raw"$"
       ]))

FLOAT_STR = (
        raw"\d+[.]?[Ee]?[+-]?\d*"
        )
INTERNAL_NODE_ANNOTATION_STR = (
        raw"\[&height_index=(?<height_index>\d+)," *
        raw"index_freq=(?<index_freq>" * FLOAT_STR * ")," *
        raw"[^]]*" *
        raw"index_height_mean=(?<index_height_mean>" * FLOAT_STR * ")," *
        raw"index_height_median=(?<index_height_median>" * FLOAT_STR * ")," *
        raw"[^]]*" *
        raw"index_height_eti_95={(?<index_height_eti_95_lower>" * FLOAT_STR * ")," * "(?<index_height_eti_95_upper>" * FLOAT_STR * ")}," *
        raw"index_height_hpdi_95={(?<index_height_hpdi_95_lower>" * FLOAT_STR * ")," * "(?<index_height_hpdi_95_upper>" * FLOAT_STR * ")}," *
        raw"[^]]*\]"
       )
INTERNAL_NODE_ANNOTATION_PATTERN = Regex(INTERNAL_NODE_ANNOTATION_STR)


function pretty_sci_not(n::Float64)::String
    n_str = @sprintf("%.2g", n)
    m = match(SCI_NOTATION_PATTERN, n_str)
    if ! isnothing(m)
        e_sign = ""
        if m[:sign] == "-"
            e_sign = "-"
        end
        n_str = "$(m[:digits]) \\times 10^{$(e_sign)$(m[:exponent])}"
    end
    return n_str
end

function get_pbs_header(pbs_script_path::AbstractString;
        exe_name::AbstractString = "phycoeval",
        exe_var_name::AbstractString = "exe_path")::AbstractString
    script_dir = Base.Filesystem.dirname(Base.Filesystem.abspath(pbs_script_path))
    relative_project_dir = Base.Filesystem.relpath(PROJECT_DIR, script_dir)
    return  """#! /bin/bash
               
               set -e
               
               if [ -n "\$PBS_JOBNAME" ]
               then
                   if [ -f "\${PBS_O_HOME}/.bashrc" ]
                   then
                       source "\${PBS_O_HOME}/.bashrc"
                   fi
                   cd "$script_dir"
               else
                   cd "\$( dirname "\${BASH_SOURCE[0]}" )"
               fi

               project_dir="$relative_project_dir"
               $exe_var_name="\${project_dir}/bin/$exe_name"
               
               if [ ! -x "\$$exe_var_name" ]
               then
                   echo "ERROR: No executable '\${$exe_var_name}'."
                   echo "       You probably need to run the project setup script."
                   exit 1
               fi
               
               source "\${project_dir}/modules-to-load.sh" >/dev/null 2>&1 || echo "    No modules loaded"
               
               """
end

file_path_iter(directory::AbstractString,
               regex_pattern::Regex
              ) = Channel(ctype = AbstractString) do c
    for (dir_path, dir_names, file_names) in Base.Filesystem.walkdir(directory)
        for f_name in file_names
            m = match(regex_pattern, f_name)
            if ! isnothing(m)
                path = Base.Filesystem.joinpath(dir_path, f_name)
                push!(c, path)
            end
        end
    end
end

dir_path_iter(directory::AbstractString,
               regex_pattern::Regex
              ) = Channel(ctype = AbstractString) do c
    for (dir_path, dir_names, file_names) in Base.Filesystem.walkdir(directory)
        for d_name in dir_names
            m = match(regex_pattern, d_name)
            if ! isnothing(m)
                path = Base.Filesystem.joinpath(dir_path, d_name)
                push!(c, path)
            end
        end
    end
end

flat_file_path_iter(directory::AbstractString,
                    regex_pattern::Regex
                   ) = Channel(ctype = AbstractString) do c
    for file_name in Base.Filesystem.readdir(directory)
        m = match(regex_pattern, file_name)
        if ! isnothing(m)
            path = Base.Filesystem.joinpath(directory, file_name)
            push!(c, path)
        end
    end
end

function simphycoeval_config_iter(
        sim_directory::AbstractString = Nothing
       )::Channel{AbstractString}
    if isnothing(sim_directory)
        sim_directory = SIM_DIR
    end
    return file_path_iter(sim_directory, SIMPHYCOEVAL_CONFIG_NAME_PATTERN)
end

function batch_dir_iter(directory::AbstractString = Nothing
              )::Channel{AbstractString}
    if isnothing(directory)
        directory = SIM_DIR
    end
    return dir_path_iter(directory, BATCH_DIR_ENDING_PATTERN)
end

function get_data_frame(paths::Vector{String};
                        skip::Int = 0)::DataFrame
    try
        return DataFrame(mapreduce(
                x -> CSV.File(x, header = 1, skipto = skip + 2),
                vcat,
                paths))
    catch e
        write(stderr, "ProjectUtil::get_data_frame:: Error parsing paths:\n")
        for pth in paths
            write(stderr, "  $pth\n")
        end
        throw(e)
    end
end

function get_data_frame_from_gz(paths::Vector{String};
                        skip::Int = 0)::DataFrame
    try
        return DataFrame(mapreduce(
                x -> CSV.File(transcode(GzipDecompressor, Mmap.mmap(x)),
                              header = 1,
                              skipto = skip + 2),
                vcat,
                paths))
    catch e
        write(stderr, "ProjectUtil::get_data_frame_from_gz:: Error parsing paths:\n")
        for pth in paths
            write(stderr, "  $pth\n")
        end
        throw(e)
    end
end

function parse_runtime(path::String)::String
    open(path, "r") do istream
        for line in eachline(istream)
            m = match(RUNTIME_PATTERN, chomp(line))
            if ! isnothing(m)
                return m[:runtime]
            end
        end
        error("Could not find runtime in '$path'")
    end
end

function get_batch_id(rng::MersenneTwister,
        num_id_digits::Int = 9)::AbstractString
    batch_num_str = Random.randstring(rng, '0':'9',
            num_id_digits)
    return batch_num_str
end


"""
Calculate Monte Carlo standard error.

Adapted from 'mcse' function of the 'mcmcse' R package (Gnu GPL version 2).
See:

Flegal, J. M. and Jones, G. L. (2010) Batch means and spectral variance
estimators in Markov chain Monte Carlo. The Annals of Statistics,
38:1034--1070.
"""
struct MonteCarloSummary
    mean::Float64
    variance::Float64
    std_error::Float64
    function MonteCarloSummary(samples::Vector{T}) where T <: Number
        n = length(samples)
        b = floor(Int, sqrt(n))
        a = floor(Int, n / b)
        y = Vector{Float64}(undef, a)
        for k in 1:a
            count = 0
            z_sum = 0.0
            for i in (((k-1)*b)+1):(k*b)
                z_sum += samples[i]
                count += 1
            end
            y[k] = z_sum / count
        end
        s_sum = sum(samples)
        mu_hat = s_sum / n
        sq_diffs = Vector{Float64}(undef, length(y))
        for i in 1:length(y)
            sq_diffs[i] = (y[i] - mu_hat)^2
        end
        sum_sq_diffs = sum(sq_diffs)
        var_hat = b * sum_sq_diffs / (a - 1.0)
        se = sqrt(var_hat / n)
        new(mu_hat, Statistics.var(samples), se)
    end
end


"""
    effective_sample_size([r1, r2, ...])

Estimate effective sample size of MCMC sample.

Adapted from 'ess' function of the 'mcmcse' R package (Gnu GPL version 2).
See:

    Gong, Lei, and James M. Flegal. A practical sequential stopping rule for
    high-dimensional MCMC and its application to spatial-temporal Bayesian
    models. arXiv:1403.5536v1 [stat.CO].

# Examples
```jldoctest
x = [0.7977061294666541, 0.9150350307910423, 0.7209626707423714,
        0.5954848559944081, 0.18032194756853182, 0.210042410144069,
        0.3673333965443635, 0.8740467791825761, 0.6874289295702046,
        0.22144353794416716, 0.3233467553676893, 0.10398479380458114,
        0.5243615565040305, 0.5877894894599294, 0.42089823773318724,
        0.6266108731616019, 0.3343859686141625, 0.512551474670303,
        0.6446230257104236, 0.36282234951752024, 0.6228723575494212,
        0.7568718761184856, 0.3718316658814024, 0.6861537858829704,
        0.1257109245390987, 0.6412426639048084, 0.48211219814972295,
        0.593973829940721, 0.4036132973697879, 0.42477867300229544,
        0.31213513805943194, 0.7963559245316685, 0.6941826857155579,
        0.6805456463190873, 0.49143482763009017, 0.6290575158052324,
        0.08661756315679514, 0.7995156973771527, 0.27539069568104,
        0.3139293111140057, 0.32288336271807183, 0.2612070751385418,
        0.4545704301079062, 0.6359171147861155, 0.3737093467417866,
        0.9232642159501455, 0.8271543021690014, 0.34958286197540656,
        0.34815266170044323, 0.6056828909353177, 0.5011441473017468,
        0.8184372611091862, 0.06710536859043326, 0.019983484122365947,
        0.3176095570458911, 0.9800154385339, 0.5319803418547973,
        0.2523950819849151, 0.04169284733227552, 0.5240020836881362,
        0.040929832798068166, 0.5024077861662805, 0.7176655502585366,
        0.6306537858831496, 0.5774716670659389, 0.9104292864296849,
        0.35302437929192343, 0.8624334312505447, 0.6990861575487167,
        0.8394941343135478, 0.5795304077084198, 0.12535068024747653,
        0.7025132099214821, 0.177220279120623, 0.9070732428670005,
        0.7666417009611808, 0.08750652002252135, 0.9948532901833365,
        0.44265582277400917, 0.10322490371849158, 0.5094288068541217,
        0.13640416841602576, 0.20328541281100587, 0.7289261198868512,
        0.8040861608469766, 0.9670617517210303, 0.23243617749946088,
        0.25068739997092004, 0.2742590187495584, 0.307652725552081,
        0.8997811130977051, 0.35615376615317706, 0.0211059298791072,
        0.03263965076194353, 0.4416542975034954, 0.5586675733736068,
        0.21167935845287156, 0.47810451475326077, 0.7395889690656308,
        0.24135469373818985]
# expected results calculated with R package mcmcse: Monte Carlo
# Standard Errors for MCMC Version 1.2-1 created on 2016-03-24.
ess1 = effective_sample_size(x, false)
ess2 = effective_sample_size(x, true)
isapprox(ess1, 154.581627605617) & isapprox(ess2, 100.0)

# output

true
```
"""
function effective_sample_size(
        samples::Vector{T},
        limit_to_number_of_samples::Bool = true
       )::Float64 where T <: Number
    sample_size = length(samples)
    if sample_size == 0
        return 0.0
    end
    mcs = MonteCarloSummary(samples)
    if mcs.std_error == 0.0
        return 0.0
    end
    ess = ((sample_size * mcs.variance) / ((mcs.std_error^2) * sample_size))
    if (ess > sample_size) & limit_to_number_of_samples
        return sample_size * 1.0
    end
    return ess
end


"""
    effective_sample_size([[r1, ....], ...])

Estimate effective sample size of MCMC sample.

Adapted from 'ess' function of the 'mcmcse' R package (Gnu GPL version 2).
See:

    Gong, Lei, and James M. Flegal. A practical sequential stopping rule for
    high-dimensional MCMC and its application to spatial-temporal Bayesian
    models. arXiv:1403.5536v1 [stat.CO].

# Examples
```jldoctest
x = [[0.7977061294666541, 0.9150350307910423, 0.7209626707423714,
        0.5954848559944081, 0.18032194756853182, 0.210042410144069,
        0.3673333965443635, 0.8740467791825761, 0.6874289295702046,
        0.22144353794416716, 0.3233467553676893, 0.10398479380458114,
        0.5243615565040305, 0.5877894894599294, 0.42089823773318724,
        0.6266108731616019, 0.3343859686141625, 0.512551474670303,
        0.6446230257104236, 0.36282234951752024, 0.6228723575494212,
        0.7568718761184856, 0.3718316658814024, 0.6861537858829704,
        0.1257109245390987, 0.6412426639048084, 0.48211219814972295],
        [0.593973829940721, 0.4036132973697879, 0.42477867300229544,
        0.31213513805943194, 0.7963559245316685, 0.6941826857155579,
        0.6805456463190873, 0.49143482763009017, 0.6290575158052324,
        0.08661756315679514, 0.7995156973771527, 0.27539069568104,
        0.3139293111140057, 0.32288336271807183, 0.2612070751385418,
        0.4545704301079062, 0.6359171147861155, 0.3737093467417866,
        0.9232642159501455, 0.8271543021690014, 0.34958286197540656,
        0.34815266170044323, 0.6056828909353177, 0.5011441473017468,
        0.8184372611091862, 0.06710536859043326, 0.019983484122365947,
        0.3176095570458911, 0.9800154385339, 0.5319803418547973,
        0.2523950819849151, 0.04169284733227552, 0.5240020836881362,
        0.040929832798068166, 0.5024077861662805, 0.7176655502585366],
        [0.6306537858831496, 0.5774716670659389, 0.9104292864296849,
        0.35302437929192343, 0.8624334312505447, 0.6990861575487167,
        0.8394941343135478, 0.5795304077084198, 0.12535068024747653,
        0.7025132099214821, 0.177220279120623, 0.9070732428670005,
        0.7666417009611808, 0.08750652002252135, 0.9948532901833365,
        0.44265582277400917, 0.10322490371849158, 0.5094288068541217,
        0.13640416841602576, 0.20328541281100587, 0.7289261198868512,
        0.8040861608469766, 0.9670617517210303, 0.23243617749946088,
        0.25068739997092004, 0.2742590187495584, 0.307652725552081,
        0.8997811130977051, 0.35615376615317706, 0.0211059298791072,
        0.03263965076194353, 0.4416542975034954, 0.5586675733736068,
        0.21167935845287156, 0.47810451475326077, 0.7395889690656308,
        0.24135469373818985]]
# expected results calculated with R package mcmcse: Monte Carlo
# Standard Errors for MCMC Version 1.2-1 created on 2016-03-24.
ess1 = effective_sample_size(x, false)
ess2 = effective_sample_size(x, true)
isapprox(ess1, 154.581627605617) & isapprox(ess2, 100.0)

# output

true
```
"""
function effective_sample_size(
        samples::Vector{Vector{T}},
        limit_to_number_of_samples::Bool = true
       )::Float64 where T <: Number
    s::Vector{Real} = vcat((v for v in samples)...)
    return effective_sample_size(s, limit_to_number_of_samples)
end


"""
    effective_sample_size(data_frame, "column_label")

Estimate effective sample size of MCMC sample.

Adapted from 'ess' function of the 'mcmcse' R package (Gnu GPL version 2).
See:

    Gong, Lei, and James M. Flegal. A practical sequential stopping rule for
    high-dimensional MCMC and its application to spatial-temporal Bayesian
    models. arXiv:1403.5536v1 [stat.CO].

# Examples
```jldoctest
df = DataFrame(A = [0.7977061294666541, 0.9150350307910423, 0.7209626707423714,
        0.5954848559944081, 0.18032194756853182, 0.210042410144069,
        0.3673333965443635, 0.8740467791825761, 0.6874289295702046,
        0.22144353794416716, 0.3233467553676893, 0.10398479380458114,
        0.5243615565040305, 0.5877894894599294, 0.42089823773318724,
        0.6266108731616019, 0.3343859686141625, 0.512551474670303,
        0.6446230257104236, 0.36282234951752024, 0.6228723575494212,
        0.7568718761184856, 0.3718316658814024, 0.6861537858829704,
        0.1257109245390987, 0.6412426639048084, 0.48211219814972295,
        0.593973829940721, 0.4036132973697879, 0.42477867300229544,
        0.31213513805943194, 0.7963559245316685, 0.6941826857155579,
        0.6805456463190873, 0.49143482763009017, 0.6290575158052324,
        0.08661756315679514, 0.7995156973771527, 0.27539069568104,
        0.3139293111140057, 0.32288336271807183, 0.2612070751385418,
        0.4545704301079062, 0.6359171147861155, 0.3737093467417866,
        0.9232642159501455, 0.8271543021690014, 0.34958286197540656,
        0.34815266170044323, 0.6056828909353177, 0.5011441473017468,
        0.8184372611091862, 0.06710536859043326, 0.019983484122365947,
        0.3176095570458911, 0.9800154385339, 0.5319803418547973,
        0.2523950819849151, 0.04169284733227552, 0.5240020836881362,
        0.040929832798068166, 0.5024077861662805, 0.7176655502585366,
        0.6306537858831496, 0.5774716670659389, 0.9104292864296849,
        0.35302437929192343, 0.8624334312505447, 0.6990861575487167,
        0.8394941343135478, 0.5795304077084198, 0.12535068024747653,
        0.7025132099214821, 0.177220279120623, 0.9070732428670005,
        0.7666417009611808, 0.08750652002252135, 0.9948532901833365,
        0.44265582277400917, 0.10322490371849158, 0.5094288068541217,
        0.13640416841602576, 0.20328541281100587, 0.7289261198868512,
        0.8040861608469766, 0.9670617517210303, 0.23243617749946088,
        0.25068739997092004, 0.2742590187495584, 0.307652725552081,
        0.8997811130977051, 0.35615376615317706, 0.0211059298791072,
        0.03263965076194353, 0.4416542975034954, 0.5586675733736068,
        0.21167935845287156, 0.47810451475326077, 0.7395889690656308,
        0.24135469373818985])
# expected results calculated with R package mcmcse: Monte Carlo
# Standard Errors for MCMC Version 1.2-1 created on 2016-03-24.
ess1 = effective_sample_size(df, "A", false)
ess2 = effective_sample_size(df, "A", true)
isapprox(ess1, 154.581627605617) & isapprox(ess2, 100.0)

# output

true
```
"""
function effective_sample_size(
        data_frame::DataFrame,
        column_label::String,
        limit_to_number_of_samples::Bool = true)::Float64
    return effective_sample_size(
            data_frame[!, column_label],
            limit_to_number_of_samples)
end


"""
    effective_sample_size([df1, ...], "column_label")

Estimate effective sample size of MCMC sample.

Adapted from 'ess' function of the 'mcmcse' R package (Gnu GPL version 2).
See:

    Gong, Lei, and James M. Flegal. A practical sequential stopping rule for
    high-dimensional MCMC and its application to spatial-temporal Bayesian
    models. arXiv:1403.5536v1 [stat.CO].

# Examples
```jldoctest
df1 = DataFrame(A = [0.7977061294666541, 0.9150350307910423,
        0.7209626707423714, 0.5954848559944081, 0.18032194756853182,
        0.210042410144069, 0.3673333965443635, 0.8740467791825761,
        0.6874289295702046, 0.22144353794416716, 0.3233467553676893,
        0.10398479380458114, 0.5243615565040305, 0.5877894894599294,
        0.42089823773318724, 0.6266108731616019, 0.3343859686141625,
        0.512551474670303, 0.6446230257104236, 0.36282234951752024,
        0.6228723575494212, 0.7568718761184856, 0.3718316658814024,
        0.6861537858829704, 0.1257109245390987])
df2 = DataFrame(A = [0.6412426639048084, 0.48211219814972295,
        0.593973829940721, 0.4036132973697879, 0.42477867300229544,
        0.31213513805943194, 0.7963559245316685, 0.6941826857155579,
        0.6805456463190873, 0.49143482763009017, 0.6290575158052324,
        0.08661756315679514, 0.7995156973771527, 0.27539069568104,
        0.3139293111140057, 0.32288336271807183, 0.2612070751385418,
        0.4545704301079062, 0.6359171147861155, 0.3737093467417866,
        0.9232642159501455, 0.8271543021690014, 0.34958286197540656,
        0.34815266170044323, 0.6056828909353177])
df3 = DataFrame(A = [0.5011441473017468, 0.8184372611091862,
        0.06710536859043326, 0.019983484122365947, 0.3176095570458911,
        0.9800154385339, 0.5319803418547973, 0.2523950819849151,
        0.04169284733227552, 0.5240020836881362, 0.040929832798068166,
        0.5024077861662805, 0.7176655502585366, 0.6306537858831496,
        0.5774716670659389, 0.9104292864296849, 0.35302437929192343,
        0.8624334312505447, 0.6990861575487167, 0.8394941343135478,
        0.5795304077084198, 0.12535068024747653, 0.7025132099214821,
        0.177220279120623, 0.9070732428670005])
df4 = DataFrame(A = [0.7666417009611808, 0.08750652002252135,
        0.9948532901833365, 0.44265582277400917, 0.10322490371849158,
        0.5094288068541217, 0.13640416841602576, 0.20328541281100587,
        0.7289261198868512, 0.8040861608469766, 0.9670617517210303,
        0.23243617749946088, 0.25068739997092004, 0.2742590187495584,
        0.307652725552081, 0.8997811130977051, 0.35615376615317706,
        0.0211059298791072, 0.03263965076194353, 0.4416542975034954,
        0.5586675733736068, 0.21167935845287156, 0.47810451475326077,
        0.7395889690656308, 0.24135469373818985])
dfs = [df1, df2, df3, df4]
# expected results calculated with R package mcmcse: Monte Carlo
# Standard Errors for MCMC Version 1.2-1 created on 2016-03-24.
ess1 = effective_sample_size(dfs, "A", false)
ess2 = effective_sample_size(dfs, "A", true)
isapprox(ess1, 154.581627605617) & isapprox(ess2, 100.0)

# output

true
```
"""
function effective_sample_size(
        data_frames::Vector{DataFrame},
        column_label::String,
        limit_to_number_of_samples::Bool = true)::Float64
    s::Vector{Real} = vcat((d[!, column_label] for d in data_frames)...)
    return effective_sample_size(s, limit_to_number_of_samples)
end



"""
    potential_scale_reduction_factor([[r11, ...], [r21, ...], ...])

Calculate the potential scale reduction factor.

Returns the square root of Equation 1.1 in:

Brooks, Stephen P. and Andrew Gelman. 1998. General Methods for Monitoring
Convergence of Iterative Simulations. Journal of Computational and
Graphical Statistics, Volume7, Number 4, Pages 434-455.

# Examples
```jldoctest
chains = [[1.1, 1.3, 1.2, 1.6, 1.5],
          [1.2, 1.7, 1.5, 1.9, 1.6]]
psrf = potential_scale_reduction_factor(chains)
# expectation calculated with commit aa83c8cc8584ba2d
# of pymc.diagnostics.gelman_rubin
# <https://github.com/pymc-devs/pymc/blob/master/pymc/diagnostics.py>
e_pymc = 1.2591483413222384
isapprox(psrf, e_pymc)

# output

true
```

```jldoctest
chains = [[1.1, 1.3, 1.2, 1.6, 1.5],
          [1.1, 1.3, 1.2, 1.6, 1.5]]
psrf = potential_scale_reduction_factor(chains)
# expectation calculated with commit aa83c8cc8584ba2d
# of pymc.diagnostics.gelman_rubin
# <https://github.com/pymc-devs/pymc/blob/master/pymc/diagnostics.py>
e_pymc = 0.89442719099991586
isapprox(psrf, e_pymc)

# output

true
```
"""
function potential_scale_reduction_factor(
        chains::Vector{Vector{T}}
       )::Float64 where T <: Number
    nchains = length(chains)
    @assert nchains > 1
    nsamples = length(chains[1])
    means = Vector{Float64}(undef, nchains)
    variances = Vector{Float64}(undef, nchains)
    for i in 1:nchains
        means[i] = Statistics.mean(chains[i])
        variances[i] = Statistics.var(chains[i])
    end
    within_chain_var = Statistics.mean(variances)
    between_chain_var = Statistics.var(means)
    pooled_var_term1 = (1.0 - (1.0 / nsamples)) * within_chain_var
    pooled_var = pooled_var_term1 + between_chain_var
    pooled_posterior_var = pooled_var + (between_chain_var / nchains)
    if within_chain_var == 0.0
        return Inf
    end
    return sqrt(pooled_posterior_var / within_chain_var)
end

"""
    potential_scale_reduction_factor([df1, ...], "column_label")

Calculate the potential scale reduction factor.

Returns the square root of Equation 1.1 in:

Brooks, Stephen P. and Andrew Gelman. 1998. General Methods for Monitoring
Convergence of Iterative Simulations. Journal of Computational and
Graphical Statistics, Volume7, Number 4, Pages 434-455.

# Examples
```jldoctest
df1 = DataFrame(A = [1.1, 1.3, 1.2, 1.6, 1.5])
df2 = DataFrame(A = [1.2, 1.7, 1.5, 1.9, 1.6])
psrf = potential_scale_reduction_factor([df1, df2], "A")
# expectation calculated with commit aa83c8cc8584ba2d
# of pymc.diagnostics.gelman_rubin
# <https://github.com/pymc-devs/pymc/blob/master/pymc/diagnostics.py>
e_pymc = 1.2591483413222384
isapprox(psrf, e_pymc)

# output

true
```

```jldoctest
df1 = DataFrame(A = [1.1, 1.3, 1.2, 1.6, 1.5])
df2 = DataFrame(A = [1.1, 1.3, 1.2, 1.6, 1.5])
psrf = potential_scale_reduction_factor([df1, df2], "A")
# expectation calculated with commit aa83c8cc8584ba2d
# of pymc.diagnostics.gelman_rubin
# <https://github.com/pymc-devs/pymc/blob/master/pymc/diagnostics.py>
e_pymc = 0.89442719099991586
isapprox(psrf, e_pymc)

# output

true
```
"""
function potential_scale_reduction_factor(
        data_frames::Vector{DataFrame},
        column_label::String)::Float64
    chains = [d[!, column_label] for d in data_frames]
    return potential_scale_reduction_factor(chains)
end

export effective_sample_size
export potential_scale_reduction_factor

end # ProjectUtil module
