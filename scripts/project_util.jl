#! /usr/bin/env julia

module ProjectUtil

# Project paths
SCRIPT_DIR = Base.Filesystem.dirname(@__FILE__)
PROJECT_DIR = Base.Filesystem.dirname(SCRIPT_DIR)
CONFIG_DIR = Base.Filesystem.joinpath(PROJECT_DIR, "configs")
BIN_DIR = Base.Filesystem.joinpath(PROJECT_DIR, "bin")
DATA_DIR = Base.Filesystem.joinpath(PROJECT_DIR, "data")
SIM_DIR = Base.Filesystem.joinpath(PROJECT_DIR, "simulations")
SIM_SCRIPT_DIR = Base.Filesystem.joinpath(SCRIPT_DIR, "simphycoeval-scripts")

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

BATCH_DIR_PATTERN_STR = raw"batch-(?<batch_num>\d+)"
BATCH_DIR_PATTERN = Regex(
        raw"^" * BATCH_DIR_PATTERN_STR * raw"$")
BATCH_DIR_ENDING_PATTERN = Regex(
            raw"^.*" * BATCH_DIR_PATTERN_STR * raw"(" * Base.Filesystem.path_separator * raw")?$")

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
               fi
               
               cd $script_dir
               
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
    return file_path_iter(directory, BATCH_DIR_ENDING_PATTERN)
end

end # ProjectUtil module
