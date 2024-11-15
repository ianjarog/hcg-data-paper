import os
#import yaml
 
# Load figure names from the YAML file
#with open("config/figures.yaml", 'r') as f:
#    figures = yaml.safe_load(f)
#
#figure_files, figure_files2 = (["outputs/publication/figures/" +
#    fig for fig in figures[key]] for key in ['figure1', 'figure2'])
#figures = figure_files + figure_files2
#suffix = 'segmented_pv_view_mom0.pdf'
#position = next((i for i, filename in enumerate(figures) if filename.endswith(suffix)), -1)
## figures to be made before this rule is executed
#required_fig = figures[:position]

galaxy_indices = ["hcg16", "hcg30", "hcg31", "hcg90", "hcg91", "hcg97"]

rule plot_segmented_pvdiag:
    input:
        expand("outputs/publication/figures/{gal_idx}_mom0_pbc_zoom.pdf", gal_idx=galaxy_indices),
        expand("outputs/publication/figures/{gal_idx}_mom0_pbc_large.pdf", gal_idx=galaxy_indices),
        expand("outputs/publication/figures/{gal_idx}_coldens_pbc_zoom.pdf", gal_idx=galaxy_indices),
        expand("outputs/publication/figures/{gal_idx}_coldens_pbc_large.pdf", gal_idx=galaxy_indices),
        expand("outputs/publication/figures/{gal_idx}_mom1_pbc_zoom.pdf", gal_idx=galaxy_indices),
        expand("outputs/publication/figures/{gal_idx}_mom1_pbc_large.pdf", gal_idx=galaxy_indices),
        mom0 = "outputs/sofia_pbc/hcg16/hcg16_line60_masked.pb_corr_vopt_mom0_nan_zoom.fits",
        meerkat_cube = "outputs/data/hcg16/hcg16_line60_masked.pb_corr_vopt.fits",
        vla_cube = "data/hcg16/vla_data/HCG16_CD_rob2_MS.fits"
        #pdfs = required_fig
    output:
        "outputs/publication/figures/hcg16_segmented_pv_view_mom0.pdf",
        "outputs/publication/figures/hcg16_segmented_pv.pdf",
        "outputs/publication/figures/hcg16_segmented_pv_vla.pdf",
    priority: 92
    shell:
        """
        
        python workflow/scripts/plot_segmented_pvdiag.py -m {input.mom0} -c {input.meerkat_cube} -v {input.vla_cube} 
        """
