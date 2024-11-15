import os
import sys
import yaml
#from workflow.scripts.utility_functions import extract_outputs_from_rules, replace_galaxy_indices


#with open("config/figures.yaml", 'r') as f:
#    figures = yaml.safe_load(f)
#
#figure_files, figure_files2 = (["outputs/publication/figures/" +
#    fig for fig in figures[key]] for key in ['figure1', 'figure2'])
#figures = figure_files + figure_files2

#suffix = '_meerkat_gbt_spec.pdf'
#position = next((i for i, filename in enumerate(figures) if filename.endswith(suffix)), -1)
## figures to be made before this rule is executed
#required_fig = figures[:position]

galaxy_indices = ["hcg16", "hcg30", "hcg31", "hcg90", "hcg91", "hcg97"]

# Retrieve outputs for each galaxy index
#sofia_outs = extract_outputs_from_rules("workflow/rules/run_sofia.smk")

#sofia_outputs = replace_galaxy_indices(sofia_outs, galaxy_indices)
#sofia_products = extract_outputs_from_rules("workflow/rules/plot_sofia_products.smk")

rule plot_gbt_meerkat_spectra:
    input:
        cube=(lambda wildcards: os.path.join("outputs", "data", wildcards.gal_idx) + "/" 
              + os.path.basename(config['sofia_params'][wildcards.gal_idx]['datacube_pbc'][:-5]+'_vopt.fits')),
    output:
        "outputs/data/{gal_idx}/{gal_idx}_line60_masked.pb_corr_vopt_times_beam_factor_integrated_spec.txt",
        "outputs/publication/figures/{gal_idx}_meerkat_gbt_spec.pdf"
    priority: 98
    log:
       "outputs/sofia/logs/{gal_idx}/{gal_idx}_meerkat_gbt_spec.log" 
    shell:
        ("python workflow/scripts/plot_gbt_meerkat_spectra.py"
         " --cube {input.cube}")
