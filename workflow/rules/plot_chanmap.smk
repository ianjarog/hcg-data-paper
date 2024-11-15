import os
import yaml
import subprocess 
 
#Load figure names from the YAML file
with open("config/figures.yaml", 'r') as f:
    figures = yaml.safe_load(f)

figure_files, figure_files2 = (["outputs/publication/figures/" +
    fig for fig in figures[key]] for key in ['figure1', 'figure2'])
figures = figure_files + figure_files2
suffix = 'hcg16_chanmap170_pbc_zoom.pdf'
position = next((i for i, filename in enumerate(figures) if filename.endswith(suffix)), -1)
## figures to be made before this rule is executed
#required_fig = figures[:position]

def get_path(gal):
    cube_path = os.path.join("outputs", "data", gal)
    noise_path = os.path.join("outputs", "sofia_pbc", gal)
    file_name = os.path.basename(config['sofia_params'][gal]['datacube_pbc'][:-5] 
                                        + f'_vopt.fits')
    noise_cube = file_name[:-5] + '_noise.fits' 
    cube_filename =  os.path.join(cube_path, file_name)
    noise_filename =  os.path.join(noise_path, noise_cube)
    return (cube_filename, noise_filename)


rule plot_chanmap:
    input:
        cube = [get_path(id)[0] for id in config['sofia_params'].keys()],
        noise = [get_path(id)[1] for id in config['sofia_params'].keys()],
        #pdfs = required_fig
    output:
        output_figs = figures[position:]
    priority: 90
    run:
        for cube_path, noise_path in zip(input.cube, input.noise):
            plot_channels = (
                f"python workflow/scripts/plot_chanmap.py "
                f"-f {cube_path} -n {noise_path}"
            )  
            subprocess.run(plot_channels, shell=True, check=True)
