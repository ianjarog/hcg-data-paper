import yaml

# Load figure names from the YAML file
#with open("config/figures.yaml", 'r') as f:
#    figures = yaml.safe_load(f)
#
#figure_files, figure_files2 = (["results/publication/figures/" +
#    fig for fig in figures[key]] for key in ['figure1', 'figure2'])
#figures = figure_files + figure_files2
#suffix = '_noise_specaxis.pdf'
#position = next((i for i, filename in enumerate(figures) if filename.endswith(suffix)), -1)
## figures to be made before this rule is executed
#required_fig = figures[:position]

rule plot_noise_specaxis:
    input:
        cube=(lambda wildcards: os.path.join("outputs", "sofia", wildcards.gal_idx) + "/" 
              + os.path.basename(config['sofia_params'][wildcards.gal_idx]['datacube'][:-5]+'_vopt_noise.fits')),
        #pdfs=required_fig,
    output:
        "outputs/publication/figures/{gal_idx}_noise_specaxis.pdf"
    priority: 95
    shell:
        ("python workflow/scripts/plot_noise_specaxis.py"
         " --cube {input.cube}")
