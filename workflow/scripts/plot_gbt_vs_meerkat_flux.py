import os
import json
import yaml
import numpy 
import argparse
import figure_properties
import matplotlib.pyplot as plt

def get_args():
    '''This function parses and returns arguments passed in'''
    description = 'Select dataset to plot'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-c', '--cube', dest='cube', \
                        help='Input data cube')
    args = parser.parse_args()
    return args

custom_font = figure_properties.FontManager(
        ['/opt/fits-tools/src/analysis_tools/'], 'tex gyre heros')

fig_prop = open("config/figure_properties.yaml")
fig_prop = yaml.safe_load(fig_prop)
# Discard noise, only calculates the sum of flux within this range
velocity_range = {'hcg16': [3500, 4300], 
                  'hcg30': [3600, 5600],
                  'hcg31': [3700, 4400],
                  'hcg91': [6400, 8000],
                  'hcg97': [5600, 7600],
                  'hcg90': [2000, 3000]}

colors = {'hcg16': '#468499', 
          'hcg30': 'black', 
          'hcg31': '#7b4ba2', 
          'hcg91': 'red',
          'hcg97': 'orange',
          'hcg90': 'purple'}

shapes = {'hcg16': 'o', 
          'hcg30': 's', 
          'hcg31': 'D', 
          'hcg91': 'X',
          'hcg97': '<',
          'hcg90': 'p'}

line_of_unity = numpy.arange(0,30,1)

#Files to save the total flux of the HCGs in Jy*Km/s
gbt_spectra = "outputs/data/hcgs/gbt_data/gbt_spectra.json"
meerkat_spectra = "outputs/data/hcgs/meerkat_spectra.json"

gbt_dict = dict()
meerkat_dict = dict()

directory = os.path.dirname(gbt_spectra)
if not os.path.exists(directory):
    os.makedirs(directory)

def write_json(spectra, input_dict, telescope):
    """ Write total flux from integrated spectra into a json file
   
    Parameters
    ----------
    spectra: str 
        name of json file to write the data into
    input_dict: dict
        Dictionary to be written to the json file
    telescope: str
        Name of the observing instrument
    """
    if not os.path.exists(spectra):
        input_dict["comment"] = telescope + " integrated spectra in Jy*km/s"
        with open(spectra, 'w') as f:
            json.dump(input_dict, f, indent=4)
    else:
        with open(spectra, 'r') as f:
            data = json.load(f)
        data.update(input_dict)
        with open(spectra, 'w') as f:
            json.dump(data, f, indent=4)

def read_spectra(cube):
    """ Calculate total flux in Jy * km/s from integrated spectra
   
    Parameters
    ----------
    cube: str
        Input data cube
    """
    gal = os.path.basename(cube)[:5]
    meerkat_flux =  os.path.join('outputs/data', gal,
                    gal + '_line60_masked.pb_corr_vopt_times_beam_factor_integrated_spec.txt')
    gbt_flux =  os.path.join('outputs/data', gal, 'gbt_data',
                  gal + '_flux_smoothed.txt')
    # gbt flux in Jansky
    gbt_j = numpy.loadtxt(gbt_flux)[:,1]
    vel_gbt = numpy.loadtxt(gbt_flux)[:,0]
    minvel = velocity_range[gal][0] 
    maxvel = velocity_range[gal][1] 
    selected_range = (vel_gbt >= minvel) & (vel_gbt <= maxvel)
    selected_gbt_j = gbt_j[selected_range]
    selected_vel_gbt = vel_gbt[selected_range]
    res_gbt = abs(numpy.mean(numpy.diff(selected_vel_gbt)))
    # meerkat flux in Jansky
    meerkat_j = numpy.loadtxt(meerkat_flux)[:,1]
    vel_meerkat = numpy.loadtxt(meerkat_flux)[:,0]
    selected_range2 = (vel_meerkat >= minvel) & (vel_meerkat <= maxvel)
    selected_meerkat_j = meerkat_j[selected_range2]
    selected_vel_meerkat = vel_meerkat[selected_range2]
    res_meerkat = abs(numpy.mean(numpy.diff(selected_vel_meerkat))) 
    # GBT flux in Jansky * km/s 
    sum_gbt = numpy.sum(selected_gbt_j) * res_gbt
    # MeerKAT flux in Jansky * km/s 
    sum_meerkat = numpy.sum(selected_meerkat_j) * res_meerkat
    gbt_dict[gal] = sum_gbt
    meerkat_dict[gal] = sum_meerkat
    # Write  the total flux to a json file
    write_json(gbt_spectra, gbt_dict, telescope="GBT")
    write_json(meerkat_spectra, meerkat_dict, telescope="MeerKAT")
  

class Totalflux:
    def __init__(self, gbt_totflux, meerkat_totflux):
        self.gbt_totflux = gbt_totflux
        self.meerkat_totflux = meerkat_totflux
        with open(self.meerkat_totflux, 'r') as f:
            self.data_meerkat = json.load(f)
        with open(self.gbt_totflux, 'r') as f:
            self.data_gbt = json.load(f)

    def check_flux(self):
        """Check if the input file contains data for all HCGs 
        
        Parameters
        ---------
        input_file: str
            Json file that contains the total flux for all HCGs
        """
        if len(self.data_meerkat.keys()) == len(velocity_range.keys()) + 1:
            return True
        else: 
            return False 

    def plot_spectra(self): 
        """Plot GBT vs MeerKAT total flux 
        
        Parameters
        ---------
        sum_gbt: str
            Total flux from GBT in Jy * km/s
        sum_meerkat: str
            Total flux from MeerKAT in Jy * km/s   
        """ 
        fig = plt.figure(figsize=(fig_prop['general']['figsize'][0],
                                  fig_prop['general']['figsize'][1]))

        ax = fig.add_axes([fig_prop['general']['x0'],
                           fig_prop['general']['y0'],
                           fig_prop['general']['w'],
                           fig_prop['general']['h']])

        for gal in self.data_meerkat.keys():
            if gal != 'comment':
                gal_label = gal.upper()[:3] + ' ' + gal[3:]
                print("MeerKAT_flux/GBT_flux = ", self.data_meerkat[gal] , self.data_gbt[gal] , self.data_meerkat[gal] / self.data_gbt[gal], " for ", gal)
                ax.plot(self.data_gbt[gal], self.data_meerkat[gal], shapes[gal], 
                    color=colors[gal], ms=fig_prop['scatter_plot']['marker_size'], 
                    label=gal_label)
        plt.legend(fontsize=fig_prop['general']['axtext_fontsize'])    
        ax.plot(line_of_unity, line_of_unity, '-', color='black', 
                 lw=fig_prop['line_plot']['line_width'])
        figprop = figure_properties.Figureprop('')
        figprop.axprop(ax, xlabel=r"$\mathrm{GBT~flux~(Jy~km~s^{-1})}$",
            ylabel=r"$\mathrm{MeerKAT~flux~(Jy~km~s^{-1})}$", is_aplpy = False)
        ouput_figure = os.path.join('outputs/publication/figures', 
                                    'hcgs_meerkat_gbt_flux.pdf')
        plt.savefig(ouput_figure, bbox_inches="tight", dpi=fig_prop['general']['dpi'])

if __name__ == '__main__':
    print("Plotting GBT vs MeerKAT total flux")
    args = get_args()
    custom_font.setup_fonts()
    read_spectra(args.cube)
    total_flux = Totalflux(gbt_spectra, meerkat_spectra)
    if total_flux.check_flux():
        total_flux.plot_spectra() 
