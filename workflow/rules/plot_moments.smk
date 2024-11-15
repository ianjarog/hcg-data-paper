import os
#import yaml
import subprocess 
 
## Load figure names from the YAML file
#with open("config/figures.yaml", 'r') as f:
#    figures = yaml.safe_load(f)
#
#figure_files, figure_files2 = (["outputs/publication/figures/" +
#    fig for fig in figures[key]] for key in ['figure1', 'figure2'])
#figures = figure_files + figure_files2
#suffix = '_mom0_pbc_zoom.pdf'
#position = next((i for i, filename in enumerate(figures) if filename.endswith(suffix)), -1)
## figures to be made before this rule is executed
#required_fig = figures[:position]

def get_mom_path(wildcards, mom_type):
    base_path = os.path.join("outputs", "sofia_pbc", wildcards.gal_idx)
    file_name = os.path.basename(config['sofia_params'][wildcards.gal_idx]['datacube_pbc'][:-5] 
                                        + f'_vopt_mom{mom_type}.fits')
    return os.path.join(base_path, file_name)


rule plot_moments:
    input:
        mom0 = lambda wildcards: get_mom_path(wildcards, 0),
        mom1 = lambda wildcards: get_mom_path(wildcards, 1),
        
    output:
        "outputs/publication/figures/{gal_idx}_mom0_pbc_zoom.pdf",
        "outputs/publication/figures/{gal_idx}_mom0_pbc_large.pdf",
        "outputs/publication/figures/{gal_idx}_coldens_pbc_zoom.pdf",
        "outputs/publication/figures/{gal_idx}_coldens_pbc_large.pdf",
        "outputs/publication/figures/{gal_idx}_mom1_pbc_zoom.pdf",
        "outputs/publication/figures/{gal_idx}_mom1_pbc_large.pdf",
    priority: 93
    run:
        plot_mom0 = (
                    f"python workflow/scripts/plot_moments.py "
                    f"--fitsfile {input.mom0} "
                   )
        plot_mom1 = (
                    f"python workflow/scripts/plot_moments.py "
                    f"--fitsfile {input.mom1}"
                   )
        subprocess.run(plot_mom0, shell=True, check=True)
        subprocess.run(plot_mom1, shell=True, check=True)
