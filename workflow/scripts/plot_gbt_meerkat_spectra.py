import os
import json
import yaml
import numpy 
import argparse
import pyspeckit
import subprocess
import figure_properties
from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as unit
from datetime import datetime
import matplotlib.pyplot as plt
from astropy import constants as const
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from analysis_tools.functions import radec2deg
from astropy.convolution import convolve, Box1DKernel
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

# Maximum velocity range in kms to be considered for the GBT spectra
max_vels = {"hcg16": 5240, "hcg30": 5900, "hcg31": 5305, 
            "hcg90": 3913, "hcg91": 8438, "hcg97": 7800} 

x_lims = {"hcg16": [3000, 5000], "hcg30": [3600, 5600], 
          "hcg31": [3200, 4800], "hcg91": [6000, 8400],
          "hcg97": [5600, 7600], "hcg90": [1500, 4000]}

y_lims = {"hcg16": [-0.02, 0.12], "hcg30": [-0.003, 0.006], 
          "hcg31": [-0.05, 0.25], "hcg91": [-0.02, 0.06],
          "hcg90": [-0.01, 0.013], "hcg97": [-0.005, 0.015]}

xmajloc = {"hcg16": 500, "hcg30": 400, "hcg31": 400, 
           "hcg91": 400, "hcg90": 500, "hcg97": 400}
ymajloc = {"hcg16": 0.05, "hcg30": 0.003, "hcg31": 0.05, 
           "hcg91": 0.02, "hcg90": 0.01, "hcg97": 0.005}

gbt_pointings = open("data/hcgs/gbp_pointing_center.json")
gbt_center = json.load(gbt_pointings)

hcgs_gal = {
    "hcg16": {
        'HCG16a': {
            'ra' : 32.35312,
            'dec': -10.13622,
            'morth': 'SBab',
            'vel': 4152,
            'delta_y': -0.015,
            'delta_x': -50
        },  
        'HCG16b': {
            'ra' : 32.3361,
            'dec': -10.13309,
            'morth': 'Sab', 
            'vel': 3977,
            'delta_y': -0.015,
            'delta_x': -50
        },
        'HCG16c': {
            'ra' : 32.41071,
            'dec': -10.14637,
            'morth': 'Im',
            'vel': 3851,
            'delta_y': -0.015,
            'delta_x': -50
        },  
        'HCG16d': {
            'ra' : 32.42872,
            'dec': -10.18394,
            'morth': 'Im',
            'vel': 3847,
            'delta_y': -0.03,
            'delta_x': -65
        }   
    },
    "hcg30": {
        'HCG30a': {
            'ra' : 69.07744,
            'dec': -2.83129,
            'morth': 'SBa',
            'vel': 4697,
            'delta_y': -0.001,
            'delta_x': -50
        },
        'HCG30b': {
            'ra' : 69.12619,
            'dec': -2.86656,
            'morth': 'Sa',
            'vel': 4625,
            'delta_y': -0.001,
            'delta_x': -50
        },  
        'HCG30c': {
            'ra' : 69.097,
            'dec': -2.79985,
            'morth': 'SBbc',
            'vel': 4508,
            'delta_y': -0.001,
            'delta_x': -50
        },
        'HCG30d': {
            'ra' : 69.15276,
            'dec': -2.84302,
            'morth': 'S0', 
            'vel': 4666,
            'delta_y': -0.002,
            'delta_x': -50
        }
    },
    "hcg31": { 
        'HCG31g': {
            'ra' : 75.433384,
            'dec': -4.28875,
            'morth': 'cI',
            'vel': 4011,
            'delta_y': -0.03,
            'delta_x': -30
        },
        'HCG31q': {
            'ra' : 75.40974,
            'dec': -4.22245,
            'morth': 'cI',
            'vel': 4090,
            'delta_y': -0.05,
            'delta_x': -30
        },  
        'HCG31a': {
            'ra' : 75.41146,
            'dec': -4.25946,
            'morth': 'Sdm',
            'vel': 4042,
            'delta_y': -0.05,
            'delta_x': -30
        },
        'HCG31b': {
            'ra' : 75.39756,
            'dec': -4.26401,
            'morth': 'Sm', 
            'vel': 4171,
            'delta_y': -0.03,
            'delta_x': -30
        },
        'HCG31c': {
            'ra' : 75.40751,
            'dec': -4.2577,
            'morth': 'Im', 
            'vel': 4068,
            'delta_y': -0.03,
            'delta_x': -30
        }
    },
    'hcg90': { 
        'HCG90a': {
            'ra' : 330.50889,
            'dec': -31.87014,
            'morth': 'Sa',
            'vel': 2575,
            'delta_y': -0.003,
            'delta_x': -40
        },
        'HCG90b': {
            'ra' : 330.53631,
            'dec': -31.99068,
            'morth': 'E0',
            'vel': 2525,
            'delta_y': -0.003,
            'delta_x': -40
        },
        'HCG90c': {
            'ra' : 330.51423,
            'dec': -31.97451,
            'morth': 'E0',
            'vel': 2696,
            'delta_y': -0.003,
            'delta_x': -40
        },
        'HCG90d': {
            'ra' : 330.52602,
            'dec': -31.99423,
            'morth': 'Im',
            'vel': 2778,
            'delta_y': -0.003,
            'delta_x': -40
        }
    },
    'hcg91': { 
        'HCG91a': {
            'ra' : 332.28174,
            'dec': -27.80984,
            'morth': 'SBc',
            'vel': 7151,
            'delta_y': -0.01,
            'delta_x': -50
        },
        'HCG91b': {
            'ra' : 332.3183,
            'dec': -27.73134,
            'morth': 'Sc',
            'vel': 7196,
            'delta_y': -0.03,
            'delta_x': -50
        },
        'HCG91c': {
            'ra' : 332.30881,
            'dec': -27.78241,
            'morth': 'Sc',
            'vel': 7319,
            'delta_y': -0.01,
            'delta_x': -50
        },
        'HCG91d': {
            'ra' : 332.28567,
            'dec': -27.80086,
            'morth': 'SB0',
            'vel': 7195,
            'delta_y': -0.017,
            'delta_x': -50
        }
    },
    'hcg97': { 
        'HCG97a': {
            'ra' : 356.84591,
            'dec': -2.30096,
            'morth': 'E5',
            'vel': 6910,
            'delta_y': -0.003,
            'delta_x': -50
        },
        'HCG97b': {
            'ra' : 356.90753,
            'dec': -2.31716,
            'morth': 'Sc',
            'vel': 6940,
            'delta_y': -0.003,
            'delta_x': -30
        },
        'HCG97c': {
            'ra' : 356.84884,
            'dec': -2.35143,
            'morth': 'Sa',
            'vel': 5995,
            'delta_y': -0.003,
            'delta_x': -50
        },
        'HCG97d': {
            'ra' : 356.82867,
            'dec': -2.31329,
            'morth': 'E1',
            'vel': 6239,
            'delta_y': -0.003,
            'delta_x': -50
        },
        'HCG97e': { 
            'ra' : 356.83253,
            'dec': -2.28102,
            'morth': 'S0a',
            'vel': 6579,
            'delta_y': -0.003,
            'delta_x': -50
        }
    }
}

def gbt_flux(f_input, f_output):
    """ Converting the GBT spectrum in brightness temperature in K 
       to Jy and the frequency axis to velocity axis 
    
    Parameters
    ----------
    f_input: str
        Input spectrum in K and frequency
    f_out: str
        outout spectrum in Jy and velocity 
    """
    spectra = numpy.loadtxt(f_input)
    gain_fac = 1.65 # antenna gain factor
    speed_of_light = const.c.to('km/s').value
    #HI line rest frequency
    f0 = 1.420405751 #GHz
    GBT_vel = speed_of_light * (f0 - spectra[:,0])/spectra[:,0]
    GBT_spec = spectra[:,1] / gain_fac
    numpy.savetxt(f_output, numpy.c_[GBT_vel, GBT_spec])
    return f_output

def smooth_spectrum(f_spectrum, new_velres, gal):
    """ Smooth spectrum using boxcar kernel 
 
    Parameters
    ---------
    f_spectrum: str or arrays
        Spectrum to be smoothed
    new_velres: float
        New velocity resolution in km/s
    gal: str
        
    """
    if isinstance(f_spectrum, str):
        f = numpy.loadtxt(f_spectrum)
        f_spectrum_smo = f_spectrum[:-4] + "_smoothed.txt"
    else:
        f = f_spectrum
    x_or = f[:,0]
    y_or = f[:,1]
    x_max = max_vels[gal]
    if gal == 'hcg97':
        x_min = 5700 # flux below this value has bad baseline
    else:
        x_min = x_or.min()
    filtered_data = f[(f[:,0] >= x_min) & (f[:,0] <= x_max)]
    x = filtered_data[:,0]
    y = filtered_data[:,1]
    velocity_resolution = numpy.mean(numpy.diff(x))  # Assuming x is evenly spaced
    kernel_size = int(numpy.round(new_velres / velocity_resolution))
    if kernel_size % 2 == 0:
        kernel_size += 1
    box_kernel = Box1DKernel(abs(kernel_size))
    
    y_smoothed = convolve(y, box_kernel)
    new_x = numpy.arange(x_min, x_max, new_velres)
    interpolate_flux = interp1d(x, y_smoothed, kind='linear', fill_value="extrapolate")
    new_y = interpolate_flux(new_x)
    if isinstance(f_spectrum, str):
        numpy.savetxt(f_spectrum_smo, numpy.c_[new_x, new_y])
        return f_spectrum_smo
    else:
        return numpy.c_[new_x, new_y]

class Compareflux:
    def __init__(self, cube):
        self.cube = cube
        self.hdu = fits.open(cube)[0] 
        self.hdr = self.hdu.header
        self.bmin = self.hdr["BMIN"] * 3600
        self.bmaj = self.hdr["BMAJ"] * 3600
        self.pixel_size = abs(self.hdr['CDELT1']) * 3600
        self.image_size = (self.hdr['NAXIS2'], self.hdr['NAXIS1'])
        self.gal = os.path.basename(cube)[:5] 
        self.outputfits = cube[:-5] + '_beam_factor.fits'
        self.gbt_ra, self.gbt_dec = gbt_center[self.gal]["ra"], gbt_center[self.gal]["dec"] 
    def aips_patgn(self, fwhm=9.1, telescope='GBT'):
        """
        This function mimics the PATGN routine in AIPS (OPCODE='GAUSS'). 
        This return a fits file where the intensity at each point is  
        given by a two-dimensional normal distribution with a specified
        r.m.s. as described in AIPS. 
    
        Parameters
        ----------
        fwhm: str
            Full width at half maximum of the beam
        telescope: str
            Name of the telesope
    
        Returns
        -------
            outputfits: str
                Name of the fits file to save the intensity values
        """ 
        # Constants
        C3 = 0
        C4 = 1
        C5 = (fwhm * 60 / 2.35482)
        # Create an array to hold the pixel values
        y, x = numpy.indices(self.image_size)
        gbt_radec = radec2deg(self.gbt_ra, self.gbt_dec)
        U = skycoord_to_pixel(SkyCoord(ra=gbt_radec["ra"]*unit.deg, dec=gbt_radec["dec"]*unit.deg,
                                       frame='icrs', equinox='J2000'), WCS(self.hdr))
        x0, y0 = int(U[0]), int(U[1])
        R = numpy.sqrt((x - x0)**2 + (y - y0)**2) * self.pixel_size
    
        # Calculate the beam response f(R)
        f_R = C3 + (C4 - C3) * numpy.exp(-((0.707 / C5)**2 * R**2))
        f_R = f_R.astype('float32')
        # Create a FITS file
        hdu = fits.PrimaryHDU(f_R)
     
        hdr = hdu.header
    
        hdr.set('SIMPLE', True, ' NO ')
        hdr.set('NAXIS', 2, 'Number of array dimensions')
        hdr.set('NAXIS1', self.image_size[1], ' NO ')
        hdr.set('NAXIS2', self.image_size[0], ' NO ')
        hdr.set('EXTEND', True, 'Tables following main image')
        hdr.set('BLOCKED', True, 'Tape may be blocked')
        hdr.set('OBJECT', 'Unknown', 'Source name')
        hdr.set('TELESCOP', telescope, 'Telescope name')
        hdr.set('INSTRUME', 'Unknown', 'Instrument used for observation')
        hdr.set('OBSERVER', 'Roger Ianjamasimanana', 'Observer name')
        hdr.set('DATE-OBS', datetime.now().strftime('%Y-%m-%d'), 'Observation start date')
        hdr.set('DATE-MAP', datetime.now().strftime('%Y-%m-%d'), 'Date of last processing')
        hdr.set('BSCALE', 1, 'REAL = TAPE * BSCALE + BZERO ')
        hdr.set('BZERO', 0, ' NO ')
        hdr.set('BUNIT', 'UNDEFINED', 'Data units')
        hdr.set('EQUINOX', 2000.0, 'Epoch of RA DEC')
        hdr.set('DATAMAX', 1, 'Maximum data value')
        hdr.set('DATAMIN', 0, 'Minimum data value')
        hdr.set('CTYPE1', 'RA---SIN', 'Coordinate type for axis 1')
        hdr.set('CRVAL1', 0, 'Coordinate value at reference point on axis 1')
        hdr.set('CDELT1', self.hdr['CDELT1'], 'Coordinate increment along axis 1')
        hdr.set('CRPIX1', self.hdr['CRPIX1'], 'Reference pixel on axis 1')
        hdr.set('CROTA1', 0, 'Rotation angle of axis 1')
        hdr.set('CTYPE2', 'DEC--SIN', 'Coordinate type for axis 2')
        hdr.set('CRVAL2', 0, 'Coordinate value at reference point on axis 2')
        hdr.set('CDELT2', self.hdr['CDELT2'], 'Coordinate increment along axis 2')
        hdr.set('CRPIX2', self.hdr['CRPIX2'], 'Reference pixel on axis 2')
        hdr.set('CROTA2', 0, 'Rotation angle of axis 2')
        if os.path.isfile(self.outputfits):
            subprocess.run(['rm', self.outputfits])
        else:
            pass
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(self.outputfits, overwrite=True)
        print(telescope + " BEAM response created {0}.".format(self.outputfits))
        return self.outputfits

    def multiply_cube(self, factor_file):
        # Read the FITS data cube
        hdus = fits.open(self.cube)
        cube_data = numpy.squeeze(hdus[0].data)
        cube_header = hdus[0].header.copy() 
        # Read the 2D factor FITS file
        hdus_factor = fits.open(factor_file)
        factor_data = numpy.squeeze(hdus_factor[0].data)
    
        # Ensure that the shape of the factor matches the cube's shape
        factor_3d = numpy.repeat(factor_data[numpy.newaxis, :, :], cube_data.shape[0], axis=0)
    
        # Multiply each channel of the cube by the factor
        result_data = cube_data * factor_3d
        output_file = self.cube[:-5] + '_times_beam_factor.fits'
        # Save the result to a new FITS file
        fits.writeto(output_file, numpy.squeeze(result_data), cube_header, overwrite=True)
     
        print("Multiplication completed. Result saved to {0}".format(output_file))
        return output_file

def get_args():
    '''This function parses and returns arguments passed in'''
    description = 'Select dataset to plot'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-c', '--cube', dest='cube', \
                        help='Input data cube')
    args = parser.parse_args()
    return args

custom_font = figure_properties.FontManager(['/opt/fits-tools/src/analysis_tools/'], 'tex gyre heros')
fig_prop = open("config/figure_properties.yaml")
fig_prop = yaml.safe_load(fig_prop)

class Plotspec:
    def __init__(self, args):
        self.args = args 
        #self.font = font
    def plot_spectra(self):
        compareflux = Compareflux(self.args.cube)
        #calculate gbt beam response assuming gaussian
        beam_response = compareflux.aips_patgn()
        #multiply meerkat cube with the gbt beam response
        meerkat_mul_gbt = compareflux.multiply_cube(beam_response)
        hdu = fits.open(meerkat_mul_gbt)[0]
        img = numpy.squeeze(hdu.data)
        sumspec = numpy.nansum(numpy.nansum(img, axis=2), axis=1)
        cellsize = compareflux.pixel_size 
        bmaj = compareflux.bmin  
        bmin = compareflux.bmaj 
        beam_fac = 2*numpy.pi/((numpy.sqrt(8*numpy.log(2)))**2)
        sumspec_jy = sumspec/(beam_fac*bmaj*bmin/(cellsize*cellsize))
        data_cube = pyspeckit.Cube(meerkat_mul_gbt)
        velocities = data_cube.xarr.value / 1000 # in km/s        
        fig = plt.figure(figsize=(12, 5))  
        # Smooth meerKAT spectrum
        common_resolution = 20
        gal = self.args.cube.split("/")[-1][:5]
        meerkat_spec = numpy.c_[velocities, sumspec_jy]
        meerkat_jansky_smoothed = smooth_spectrum(meerkat_spec, common_resolution, gal)
        meerkat_spectrum_x = meerkat_jansky_smoothed[:,0]
        meerkat_spectrum_y = meerkat_jansky_smoothed[:,1]
        plt.plot(meerkat_spectrum_x, meerkat_spectrum_y, "#b40056", label="MeerKAT")
        # Save integrated spectrum of MeerkAT within GBT beam
        numpy.savetxt(meerkat_mul_gbt[:-5]+"_integrated_spec.txt", 
                numpy.c_[meerkat_spectrum_x, meerkat_spectrum_y])
        #common_resolution = abs(velocities[0] - velocities[1]) 
        # smooth gbt spectrum and plot
        output_dir = os.path.join('outputs', 'data', gal, 'gbt_data')
        input_dir = os.path.join('data', gal, 'gbt_data')
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        f_kelvin = os.path.join(input_dir, f'{gal}.ascii')
        f_jansky = output_dir + '/' + os.path.basename(f_kelvin)[:-6] + '_flux.txt'
        gbt_jansky = gbt_flux(f_kelvin, f_jansky) # convert intensity in K to Jy
        # smooth gbt spectrum to the MeerKAT resolution
        gbt_jansky_smoothed = smooth_spectrum(gbt_jansky, common_resolution, gal)
        gbt_jansky_smoothed = numpy.loadtxt(gbt_jansky_smoothed)
        gbt_spectrum_x = gbt_jansky_smoothed[:,0]
        gbt_spectrum_y = gbt_jansky_smoothed[:,1]
        plt.plot(gbt_spectrum_x, gbt_spectrum_y, "black", label="GBT")
        plt.xlim([x_lims[gal][0], x_lims[gal][1]])
        plt.ylim([y_lims[gal][0], y_lims[gal][1]])
        plt.legend(fontsize=fig_prop['general']['axtext_fontsize']-5, loc='upper right')
        # Add major and minor ticks
        ax = plt.gca() # get current axes
        ax.axhline(linestyle='dashed')
        ax.text(0.05, 0.78, gal.upper()[:3] + ' ' +  gal[3:], 
                fontsize=fig_prop['general']['axtext_fontsize']-5, 
                transform=ax.transAxes)
        figprop = figure_properties.Figureprop('')
        figprop.axprop(ax, xlabel=r"$\mathrm{Velocity~(km~s^{-1})}$",
            ylabel=r"$\mathrm{Flux~(Jy)}$", is_aplpy = False, 
            majorlocator=True, 
            xmajloc = xmajloc[gal], ymajloc= ymajloc[gal])
        for k in range(len(hcgs_gal[gal])):
            delta_y = hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["delta_y"]
            delta_x = hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["delta_x"]
            ax.vlines(hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["vel"], 
                    ymin=y_lims[gal][0], ymax=y_lims[gal][1], color="black", ls='dotted')
            ax.text(hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["vel"] + delta_x, 
                    y_lims[gal][1] + delta_y, list(hcgs_gal[gal].keys())[k][-1], color="black", 
                    fontsize=fig_prop['general']['axtext_fontsize']-5)
        figure_path = 'outputs/publication/figures/'
        plt.savefig(figure_path + gal + '_meerkat_gbt_spec.pdf', 
                bbox_inches="tight", dpi=fig_prop['general']['dpi'])


if __name__ == '__main__':
    print("Plotting GBT vs MeerKAT spectra")
    args = get_args()
    custom_font.setup_fonts()
    plotspec = Plotspec(args)
    plotspec.plot_spectra()
