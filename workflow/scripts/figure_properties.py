import yaml
import matplotlib as mpl
from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

class FontManager:
    """Replacing default matplotlib font to tex gyre hereos"""
    def __init__(self, font_dirs=['/opt/fits-tools/src/analysis_tools/'], 
                 font='tex gyre heros'):
        self.font_dirs = font_dirs
        self.font = font
        self.setup_fonts()

    def setup_fonts(self):
        for font_dir in self.font_dirs:
            font_files = font_manager.findSystemFonts(fontpaths=[font_dir])
            for font_file in font_files:
                font_manager.fontManager.addfont(font_file)
        mpl.rcParams['font.sans-serif'] = self.font
        mpl.rc('mathtext', fontset='custom', it=self.font + ':italic')

fig_prop = open("config/figure_properties.yaml")
fig_prop = yaml.safe_load(fig_prop)

class Figureprop:
    def __init__(self, fitsfig, 
            xlabelpad=fig_prop['general']['xlabelpad'],
            ylabelpad=fig_prop['general']['ylabelpad'],
            mat_xlabelpad=fig_prop['general']['mat_xlabelpad'],
            mat_ylabelpad=fig_prop['general']['mat_ylabelpad']):
        self.fitsfig = fitsfig
        self.xlabelpad = xlabelpad
        self.ylabelpad = ylabelpad
        self.mat_xlabelpad = mat_xlabelpad
        self.mat_ylabelpad = mat_ylabelpad
        self.ticklabelsize = fig_prop['general']['ticklabelsize']
        self.axislabelsize = fig_prop['general']['xylabel_fontsize']
   
    def aplpy_figprop(self):
        self.fitsfig.tick_labels.set_font(size = self.ticklabelsize)
        self.fitsfig.axis_labels.set_font(size = self.axislabelsize)
    
    def axprop(self, ax, xlabel='RA (J2000)',
            ylabel='Dec (J2000)', axnum=1, is_aplpy = True,
            axnumshow=True, secaxtickhide=False, majorlocator=False,
            xmajloc=3, ymajloc=4, show_yticks=True, show_xticks=True):
        if is_aplpy: # if ax is from aplpy
            ax.coords[axnum].set_ticklabel_visible(axnumshow)
            ax.coords[1].set_ticklabel_visible(show_yticks)
            ax.coords[0].set_ticklabel_visible(show_xticks)
            ax.coords[0].set_ticklabel(exclude_overlapping=True)
            ax.coords[axnum].set_ticks_visible(axnumshow)
            ax.set_xlabel(xlabel, labelpad=self.xlabelpad)
            ax.set_ylabel(ylabel, labelpad=self.ylabelpad)
        else:
            ax.set_xlabel(xlabel, labelpad=self.mat_xlabelpad, 
                    fontsize=self.axislabelsize)
            ax.set_ylabel(ylabel, labelpad=self.mat_ylabelpad, 
                    fontsize=self.axislabelsize)
            ax.tick_params(axis='both', which='major', 
                           labelsize=self.ticklabelsize)
        secax = ax.secondary_xaxis('top')
        secaxr = ax.secondary_yaxis('right')
        for axes in [ax, secax, secaxr]:
            axes.minorticks_on()
            axes.tick_params(which='both', direction='in',
                             length=fig_prop['axes']['major_length'], 
                             width=fig_prop['axes']['tick_width'], 
                             pad=fig_prop['axes']['tick_pad'])
            axes.tick_params(which='minor', 
                    length=fig_prop['axes']['minor_length'])
        secax.tick_params(labeltop=False)
        secaxr.tick_params(labelright=False)
        secaxr.tick_params(labelleft=False)
        if secaxtickhide:
            secax.set_xticks([])
            secaxr.set_yticks([])
        if majorlocator:
            ax.xaxis.set_major_locator(MultipleLocator(xmajloc))
            ax.yaxis.set_major_locator(MultipleLocator(ymajloc))
            secax.xaxis.set_major_locator(MultipleLocator(xmajloc))
            secaxr.yaxis.set_major_locator(MultipleLocator(ymajloc))

    def aplpy_colorbar(self, loc='right', ylabel='',
                       labelsize=fig_prop['general']['xylabel_fontsize'], 
                       labelpad=fig_prop['axes']['cbar_labelpad'], cbarloc='right'):
        self.fitsfig.add_colorbar()
        self.fitsfig.colorbar.set_location(loc)
        self.fitsfig.colorbar.show()
        self.fitsfig.colorbar.set_pad(0.0)
        self.fitsfig.colorbar.set_axis_label_pad(labelpad)
        self.fitsfig.colorbar.set_axis_label_text(ylabel)
        self.fitsfig.colorbar.set_axis_label_font(size=labelsize)
        self.fitsfig.colorbar.set_width(fig_prop['axes']['cbar_width'])
        self.fitsfig.colorbar.set_location(cbarloc)
        cbar = self.fitsfig.colorbar._colorbar
        # Set the tick label size
        cbar.ax.tick_params(labelsize=fig_prop['general']['ticklabelsize'])
        ax = plt.gca()
        ax.tick_params(direction='in', length=5, width=1, 
                pad=fig_prop['axes']['cbar_tickpad'])
        ax.tick_params(which='minor', length=1)

    @staticmethod
    def ra_dec_corner(fits_file, edge_offset_percent=10):
        """
        Calculate the RA and Dec coordinates near the lower left corner of a FITS image,
        with an offset that is a percentage of the image dimensions.

        Parameters
        ----------
        fits_file : str
            The path to the FITS file.
        edge_offset_percent : int or float, optional
            The offset from the edge as a percentage of the image dimensions.

        Returns
        -------
        ra : float
            The right ascension of the calculated position.
        dec : float
            The declination of the calculated position.
        """
        with fits.open(fits_file) as hdul:
            hdr = hdul[0].header
            wcs = WCS(hdr)
            nx, ny = hdr['NAXIS1'], hdr['NAXIS2']
            x_pixel = nx * (edge_offset_percent / 100)
            y_pixel = ny * (edge_offset_percent / 100)
            ra, dec = wcs.wcs_pix2world(x_pixel, y_pixel, 1)
            return ra, dec
