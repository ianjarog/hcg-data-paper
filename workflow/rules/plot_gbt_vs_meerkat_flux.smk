#import yaml
import os
import sys
#sys.path.insert(0, os.path.abspath('data/workflow'))
#from utility_functions import extract_outputs_from_rules, replace_galaxy_indices

#with open("config/figures.yaml", 'r') as f:
#    figures = yaml.safe_load(f)
#
#figure_files, figure_files2 = (["outputs/publication/figures/" +
#    fig for fig in figures[key]] for key in ['figure1', 'figure2'])
#figures = figure_files + figure_files2
#
#suffix = 's_meerkat_gbt_flux.pdf'
#position = next((i for i, filename in enumerate(figures) if filename.endswith(suffix)), -1)
## figures to be made before this rule is executed
#required_fig = figures[:position]

#gbt_meerkat_spectra_out = extract_outputs_from_rules("workflow/rules/plot_gbt_meerkat_spectra.smk")
#gbt_meerkat_spectra_out = replace_galaxy_indices(gbt_meerkat_spectra_out, galaxy_indices)

galaxy_indices = ["hcg16", "hcg30", "hcg31", "hcg90", "hcg91", "hcg97"]

rule plot_gbt_vs_meerkat_flux:
    input:
        expand("outputs/publication/figures/{gal_idx}_meerkat_gbt_spec.pdf", gal_idx=galaxy_indices),
        expand("outputs/data/{gal_idx}/{gal_idx}_line60_masked.pb_corr_vopt_times_beam_factor_integrated_spec.txt", gal_idx=galaxy_indices),
        cubes = [os.path.join("results", "data", idx) + "/" 
                 + os.path.basename(config['sofia_params'][idx]['datacube_pbc'][:-5] 
                 + '_vopt.fits') for idx in config['sofia_params'].keys()]
        #gbt_meerkat_spectra_out
        #pdfs = required_fig 
    output:
        "outputs/publication/figures/hcgs_meerkat_gbt_flux.pdf"
    priority:97
    shell:
        """
        
        for cube in {input.cubes}; do
            python workflow/scripts/plot_gbt_vs_meerkat_flux.py --cube $cube
        done  
        """
