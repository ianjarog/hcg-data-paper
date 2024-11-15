import os
#import yaml
#  
## Load figure names from the YAML file
#with open("config/figures.yaml", 'r') as f:
#    figures = yaml.safe_load(f)
#
#figure_files, figure_files2 = (["outputs/publication/figures/" +
#    fig for fig in figures[key]] for key in ['figure1', 'figure2'])
#figures = figure_files + figure_files2
#
#suffix = '_global_profile.pdf'
#position = next((i for i, filename in enumerate(figures) if filename.endswith(suffix)), -1)
## figures to be made before this rule is executed
#required_fig = figures[:position]

rule plot_global_prof:
    input:
        cube=(lambda wildcards: os.path.join("outputs", "data", wildcards.gal_idx) + "/" 
              + os.path.basename(config['sofia_params'][wildcards.gal_idx]['datacube_pbc'][:-5]+'_vopt.fits')),
        vlaspec=(lambda wildcards: os.path.join("data", wildcards.gal_idx) + "/vla_data/" 
              + wildcards.gal_idx + '_spec.txt'),
        mask=(lambda wildcards: os.path.join("outputs", "sofia_pbc", wildcards.gal_idx) + "/" 
              + os.path.basename(config['sofia_params'][wildcards.gal_idx]['datacube_pbc'][:-5]+'_vopt_mask.fits')),
        #pdfs = required_fig
    output:
        "outputs/publication/figures/{gal_idx}_global_profile.pdf"
    priority: 94
    shell:
        ("python workflow/scripts/plot_global_prof.py"
         " --cube {input.cube}"
         " --mask {input.mask}"
         " --vlaspec {input.vlaspec}")
