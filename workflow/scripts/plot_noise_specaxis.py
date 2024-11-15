import os
import yaml
import numpy
import argparse
import pyspeckit
import figure_properties
from astropy.io import fits
import matplotlib.pyplot as plt
from analysis_tools.functions import delheader
from utility_functions import extract_medians_from_fits

def get_args():
    '''This function parses and returns arguments passed in'''
    description = 'Select dataset to plot'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-c', '--cube', dest='cube', \
                        help='Input data cube')
    args = parser.parse_args()
    return args

xlims = {"hcg16": [2800, 5000], 
         "hcg30": [3300, 5700],
         "hcg31": [3200, 4800],
         "hcg91": [4000, 11000],
         "hcg97": [5000, 8500],
         "hcg90": [1500, 4000]}
 
ylims = {"hcg16": [0.32, 0.35], 
         "hcg30": [0.32, 0.335],
         "hcg31": [0.325, 0.35],
         "hcg91": [0.31, 0.36],
         "hcg90": [0.305, 0.325],
         "hcg97": [0.31, 0.36]}


fig = plt.figure(figsize=(10.5,10.5))
custom_font = figure_properties.FontManager(
        ['/opt/fits-tools/src/analysis_tools/'], 'tex gyre heros')
fig_prop = open("config/figure_properties.yaml")
fig_prop = yaml.safe_load(fig_prop)

def plot_noise(cube):
    cube_name = os.path.basename(cube)
    gal = cube_name[:5]
    data_cube = pyspeckit.Cube(cube)
    velocities = data_cube.xarr.value / 1000
    noise = extract_medians_from_fits(cube) * 1000
    noise_val = {gal: float(numpy.median(noise))}
    with open(cube[:-5] + '_median_noise.yaml', 'w') as val:
        yaml.dump(noise_val, val, default_flow_style=False)
    print("MEDIAN NOISE VALUE OF ", gal.upper(), numpy.median(noise))

    fig = plt.figure(figsize=(fig_prop['general']['figsize'][0],
                              fig_prop['general']['figsize'][1]))
    ax = fig.add_axes([fig_prop['general']['x0'],
                       fig_prop['general']['y0'],
                       fig_prop['general']['w'],
                       fig_prop['general']['h']])

    ax.plot(velocities, noise, '-', color='#468499', lw=1.5)
    ax.set_ylim([ylims[gal][0], ylims[gal][1]])
    ax.set_xlim([xlims[gal][0], xlims[gal][1]])
    ax.hlines(numpy.median(noise), xmin=velocities.min()-250, 
        xmax=velocities.max()+250, color="#803367", lw=2, ls='--', 
    label = ("Median: " 
             + str(round(numpy.median(noise), 2))  
             + " $\mathrm{mJy~beam^{-1}}$"))
    plt.legend(frameon=False, fontsize=fig_prop['general']['axtext_fontsize'], loc='upper left')
    figprop = figure_properties.Figureprop('')
    figprop.axprop(ax, xlabel=r"$\mathrm{Velocity~(km~s^{-1})}$",
        ylabel=r"$\mathrm{Noise~(mJy~beam^{-1})}$", is_aplpy = False)
    plt.savefig("outputs/publication/figures/" + gal + '_noise_specaxis.pdf', 
        bbox_inches="tight", dpi=fig_prop['general']['dpi'])

if __name__ == '__main__':
    custom_font.setup_fonts()
    args = get_args()
    plot_noise(args.cube) 
#python workflow/scripts/plot_noise_specaxis.py -c outputs/sofia/hcg97/hcg97_line60_masked.image_vopt_noise.fits
