import os
#import yaml
import subprocess 
 
# Load figure names from the YAML file
#with open("config/figures.yaml", 'r') as f:
#    figures = yaml.safe_load(f)
#
#figure_files, figure_files2 = (["outputs/publication/figures/" +
#    fig for fig in figures[key]] for key in ['figure1', 'figure2'])
#figures = figure_files + figure_files2
suffix = 'hcg16_optical_mom0.pdf'
#position = next((i for i, filename in enumerate(figures) if filename.endswith(suffix)), -1)
# figures to be made before this rule is executed
#required_fig = figures[:position]
def get_path(gal):
    mom0_path = os.path.join("outputs", "sofia_pbc", gal)
    file_name = os.path.basename(config['sofia_params'][gal]['datacube_pbc'][:-5] 
                                        + f'_vopt.fits')
    mom0_fitsname = file_name[:-5] + '_mom0.fits' 
    mom0_filename =  os.path.join(mom0_path, mom0_fitsname)
    return mom0_filename


rule plot_optical_mom0:
    input:
        mom0 = [get_path(gal_id) for gal_id in config['sofia_params'].keys()], 
        #pdfs = required_fig
    output:
        "outputs/publication/figures/hcg16_optical_mom0.pdf",
        "outputs/publication/figures/hcg30_optical_mom0.pdf",
        "outputs/publication/figures/hcg31_optical_mom0.pdf",
        "outputs/publication/figures/hcg90_optical_mom0.pdf",
        "outputs/publication/figures/hcg91_optical_mom0.pdf",
        "outputs/publication/figures/hcg97_optical_mom0.pdf"
    priority: 91
    shell:
        """

        for mommap in {input.mom0}; do
            python workflow/scripts/plot_optical_mom0.py -f $mommap
        done
        """
