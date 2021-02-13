#! /usr/bin/env julia

using YAML

function main_cli()::Cint
    try
        path = ARGS[1]
        d = YAML.load_file(path)
        h = d["splits"]["root"]["height_mean"]
        write(stdout, "$h\n")
    catch e
    end

    return 0
end

main_cli()
