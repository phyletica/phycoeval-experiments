#! /usr/bin/env julia

using YAML
using Printf

include("project_util.jl")


function main_cli()::Cint
    write(stdout, "Summary & \\textit{Cyrtodactylus} & \\textit{Gekko} \\\\\n")
    write(stdout, "\\hline\n")

    d = YAML.load_file("../gekkonid-output/posterior-summary-cyrt-nopoly.yml")
    cyrt_asdsf = d["summary_of_split_freq_std_deviations"]["mean"]
    cyrt_root_age_psrf = d["splits"]["root"]["height_psrf"]
    cyrt_root_age_ess = d["splits"]["root"]["height_ess"]
    cyrt_pop_size_psrf = d["splits"]["root"]["pop_size_psrf"]
    cyrt_pop_size_ess = d["splits"]["root"]["pop_size_ess"]
    cyrt_n_samples = d["summary_of_tree_sources"]["total_number_of_trees_sampled"]
    cyrt_n_chains = length(d["summary_of_tree_sources"]["sources"])
    cyrt_n_samples_per_chain = d["summary_of_tree_sources"]["sources"][1]["number_of_trees_sampled"]
    cyrt_burnin_per_chain = d["summary_of_tree_sources"]["sources"][1]["number_of_trees_skipped"]

    d = YAML.load_file("../gekkonid-output/posterior-summary-gekko-nopoly.yml")
    gekko_asdsf = d["summary_of_split_freq_std_deviations"]["mean"]
    gekko_root_age_psrf = d["splits"]["root"]["height_psrf"]
    gekko_root_age_ess = d["splits"]["root"]["height_ess"]
    gekko_pop_size_psrf = d["splits"]["root"]["pop_size_psrf"]
    gekko_pop_size_ess = d["splits"]["root"]["pop_size_ess"]
    gekko_n_samples = d["summary_of_tree_sources"]["total_number_of_trees_sampled"]
    gekko_n_chains = length(d["summary_of_tree_sources"]["sources"])
    gekko_n_samples_per_chain = d["summary_of_tree_sources"]["sources"][1]["number_of_trees_sampled"]
    gekko_burnin_per_chain = d["summary_of_tree_sources"]["sources"][1]["number_of_trees_skipped"]

    cyrt_asdsf = ProjectUtil.pretty_sci_not(cyrt_asdsf)
    cyrt_root_age_psrf = @sprintf("%.7f", cyrt_root_age_psrf)
    cyrt_root_age_ess = @sprintf("%.2f", cyrt_root_age_ess)
    cyrt_pop_size_psrf = @sprintf("%.7f", cyrt_pop_size_psrf)
    cyrt_pop_size_ess = @sprintf("%.2f", cyrt_pop_size_ess)

    gekko_asdsf = ProjectUtil.pretty_sci_not(gekko_asdsf)
    gekko_root_age_psrf = @sprintf("%.7f", gekko_root_age_psrf)
    gekko_root_age_ess = @sprintf("%.2f", gekko_root_age_ess)
    gekko_pop_size_psrf = @sprintf("%.7f", gekko_pop_size_psrf)
    gekko_pop_size_ess = @sprintf("%.2f", gekko_pop_size_ess)

    write(stdout, "Number of chains & $cyrt_n_chains & $gekko_n_chains \\\\\n")
    write(stdout, "Samples skipped per chain & $cyrt_burnin_per_chain & $gekko_burnin_per_chain \\\\\n")
    write(stdout, "Samples retained per chain & $cyrt_n_samples_per_chain & $gekko_n_samples_per_chain \\\\\n")
    write(stdout, "Total samples & $cyrt_n_samples & $gekko_n_samples \\\\\n")
    write(stdout, "ASDSF & $cyrt_asdsf & $gekko_asdsf \\\\\n")
    write(stdout, "Root age PSRF & $cyrt_root_age_psrf & $gekko_root_age_psrf \\\\\n")
    write(stdout, "Root age ESS & $cyrt_root_age_ess & $gekko_root_age_ess \\\\\n")
    write(stdout, "Population size PSRF & $cyrt_pop_size_psrf & $gekko_pop_size_psrf \\\\\n")
    write(stdout, "Population size ESS & $cyrt_pop_size_ess & $gekko_pop_size_ess \\\\\n")
    write(stdout, "\\hline\n")

    return 0
end

main_cli()
