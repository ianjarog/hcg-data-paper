import os
import yaml
import glob
import aplpy
import argparse
import numpy as np
import figure_properties
import matplotlib as mpl
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import colors
import matplotlib.pyplot as plt
from PIL import Image, ImageOps
from utility_functions import nhitomom0
from analysis_tools.functions import delheader
from analysis_tools.functions import radec2deg
from utility_functions import get_deep_optical
from astropy.visualization import (SinhStretch, ImageNormalize)


font_size=40

class Fitshandler:
    @staticmethod
    def normalize_and_save(input_file, output_file, 
            stretch=SinhStretch(), fill_value=0):
        """
        Normalize the data in a FITS file and save it to a new file.

        This function opens a FITS file, 
        normalizes its data using the specified
        stretching function, and then saves the 
        normalized data to a new FITS file.
        If the data is a masked array, masked values 
        are replaced with a specified fill value.

        Parameters
        ----------
        input_file : str
            The path to the input FITS file.
        output_file : str
            The path where the normalized FITS file will be saved.
        stretch : astropy.visualization.Stretch instance, optional
            The stretching function to be applied 
            to the data for normalization.
            Default is `SinhStretch()`.
        fill_value : float or int, optional
            The value to replace masked values with, 
            if the data is a masked array.
            Default is 0.

        Returns
        -------
        None
        """  
        with fits.open(input_file) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            norm = ImageNormalize(data, stretch=stretch)
            normalized_data = norm(data)
            if isinstance(normalized_data, np.ma.MaskedArray):
                normalized_data = normalized_data.filled(fill_value)
            hdu = fits.PrimaryHDU(data=normalized_data, header=header)
            hdu.writeto(output_file, overwrite=True)

    @staticmethod
    def get_percentile_vmin_vmax(file_name, 
            lower_percentile, upper_percentile):
        """
        Calculate and return the lower and upper 
        percentile values of the data in a FITS file.

        This function opens a FITS file, 
        extracts the data from its first HDU, and computes
        the lower and upper percentile values specified by the user. 
        It handles the data by flattening the array and 
        removing NaN values before performing the percentile calculation.

        Parameters
        ----------
        file_name : str
            The path to the FITS file from which data will be read.
        lower_percentile : float
            The lower percentile value for the data range. 
            It should be a number between 0 and 100.
        upper_percentile : float
            The upper percentile value for the data range. 
            It should be a number between 0 and 100,
            and greater than the `lower_percentile`.

        Returns
        -------
        vmin : float
            The value at the lower percentile in the data.
        vmax : float
            The value at the upper percentile in the data.

        Notes
        -----
        This function is primarily used for 
        determining suitable scaling limits (vmin, vmax)
        for visualization of FITS image data.
        """
        with fits.open(file_name) as hdul:
            data = hdul[0].data
            data_flat = data.flatten()
            data_flat = data_flat[~np.isnan(data_flat)]
            vmin = np.percentile(data_flat, lower_percentile)
            vmax = np.percentile(data_flat, upper_percentile)
            return vmin, vmax

    @staticmethod
    def dividefits(file_name, outname, value=1000):
        """ Divides values in fits file by value

        Parameters
        ---------
        file_name: str
            Input fits file
        outname: str
            Output file name
        value: float
            The value to divide the fits file with. 
        """
        hdu = fits.open(file_name)[0]
        img = hdu.data / value
        hdr = hdu.header
        fits.writeto(outname, img, hdr, overwrite=True)
        return outname

    @staticmethod
    def get_beam(fits_file):
        """Get beam info from a fits file 
        
        Parameters
        ---------
        fits_file: str
            input fits file to get beam information from

        Returns
        -------
            beaminfo: dict
                beam information, bmaj, bmin, bpa
        """
        hdu = fits.open(fits_file)[0]
        hdr = hdu.header
        try:
            bmaj = hdr["BMAJ"]
            bmin = hdr["BMIN"]
            bpa = hdr["BPA"]
            return {"bmin": bmin, "bmaj": bmaj, "bpa": bpa}
        except KeyError as e:
            raise KeyError(f"Missing {e.args[0]} keyword in FITS header.")


custom_font = figure_properties.FontManager(
        ['/opt/fits-tools/src/analysis_tools/'], 'tex gyre heros')
fig_prop = open("config/figure_properties.yaml")
fig_prop = yaml.safe_load(fig_prop)

data_r = 'normalized_r.fits'

class Plotter:
    def __init__(self, figure_name, filetoplot, label=''):
        """
        filetoplot : str
            The path to the FITS file from which contours 
            are extracted and overlaid on the RGB image.
        """
        #self.font = font
        self.filetoplot = filetoplot
        self.beaminfo = Fitshandler.get_beam(self.filetoplot)
        self.bmin, self.bmaj, self.bpa = (self.beaminfo["bmin"], 
                                          self.beaminfo["bmaj"], 
                                          self.beaminfo["bpa"])
        self.figure_name = figure_name
        self.figure_path = 'outputs/publication/figures/' 
        self.label = label

    
    def plot_momentzero(self, r_filter, contlev):
        """
        Generate and save a false-color image and contour plots from FITS files.

        This function normalizes the data from three 
        FITS files (representing i, g, and z filters),
        creates an RGB cube, and generates a false-color image. 
        It then overlays contour plots from another 
        FITS file onto this image. The resulting plot is saved as a PDF file.

        Parameters
        ----------
        r_filter : str
            The path to the FITS file for the r filter.
        contlev: array
            contour levels 

        Notes
        -----
        The function depends on the APLpy library for FITS file 
        handling and matplotlib for plotting.
        The contours' color, levels, and other properties are hard-coded 
        and may need adjustment for different datasets.

        Outputs
        -------
        False-color image with overlaid contours and
        other annotations.
        """
        fig_size = plt.figure(figsize=(fig_prop['general']['figsize'][0],
                                       fig_prop['general']['figsize'][1]))    
        # Normalization
        Fitshandler.normalize_and_save(r_filter, data_r)
        vmin_r, vmax_r = Fitshandler.get_percentile_vmin_vmax(data_r, 1, 99)
        fig = aplpy.FITSFigure(r_filter, figure=fig_size, 
                               subplot=[fig_prop['general']['x0'],
                               fig_prop['general']['y0'],
                               fig_prop['general']['w'],
                               fig_prop['general']['h']])
        figprop = figure_properties.Figureprop(fig)
        figprop.aplpy_figprop()
        fig_size.canvas.draw()
        ax = plt.gca()
        figprop.axprop(ax, secaxtickhide=True)
        ax.text(0.35, 0.90, self.label, transform=ax.transAxes,
                 color="black", fontsize=font_size, 
                 fontfamily='sans-serif')
        #fig.show_grayscale(invert="True")
        contcol = ['#2698ff', 'red']
        fig.show_contour(self.filetoplot, colors=contcol[1], 
            levels=contlev, linewidths=1.5)
        fig.show_grayscale(invert=True, vmin=vmin_r, vmax=vmax_r)
        beampos = figprop.ra_dec_corner(r_filter)
        x_word, y_word = beampos
        fig.show_ellipses(x_word, y_word, self.bmin, 
            self.bmaj, angle=self.bpa, edgecolor='black')
        plt.tight_layout()
        output_figure = self.figure_path + self.figure_name
        plt.savefig(output_figure, bbox_inches="tight", dpi=fig_prop['general']['dpi'])
        os.system("rm cube.fits false_color.png")
        os.system("rm false_color_inverted.png")

    def plot_sources(self, colmap, beamcolor, colorbar_label, 
             contlev=None, contcol=None,
             textlabelcolor="black", text_xpos=0.35):
        """
        Plot SoFiA products such as signal-to-noise ratio map and moment maps. 

        Parameters
        ----------
        figure_name: str
            Output figure name 
        colmap: matplotlib.colors.LinearSegmentedColormap 
        Color map 
        contlev: array 
            Contour levels
        contcol: str
            Contour color
        beamcolor: str
            Edge color of the beam
        textlabelcolor: str
            Color for text inside plot
        Notes
        -----
        The function depends on the APLpy library 
        for FITS file handling and matplotlib for plotting.
        The contours' color, levels, and other properties are hard-coded 
        and may need adjustment for different datasets.

        Outputs
        -------
        False-color image with overlaid contours and
        other annotations.
        """
        figure_size = (fig_prop['general']['figsize'][0], 
                       fig_prop['general']['figsize'][1])
        fig_size = plt.figure(figsize=figure_size)    
        vmin, vmax = Fitshandler.get_percentile_vmin_vmax(self.filetoplot, 1, 99)
        # Initialize an APLpy figure using one of the FITS files (for WCS info)
        fig = aplpy.FITSFigure(self.filetoplot, figure=fig_size, 
                               subplot=[fig_prop['general']['x0'],
                                        fig_prop['general']['y0'],
                                        fig_prop['general']['w'],
                                        fig_prop['general']['h']])
        # Display the RGB image
        fig.show_colorscale(vmin=vmin, vmax=vmax, 
                    aspect='auto', cmap=colmap, stretch='linear')
        figprop = figure_properties.Figureprop(fig)
        figprop.aplpy_figprop()
        fig_size.canvas.draw()
        ax = plt.gca()  
        ax.text(text_xpos,0.90, self.label, transform=ax.transAxes, 
                color=textlabelcolor, 
                fontsize=font_size, 
                fontfamily='sans-serif')
        figprop.axprop(ax, secaxtickhide=True)
        figprop.aplpy_colorbar(ylabel=colorbar_label, cbarloc="top")
        # Overlay the contour
        if contlev is not None:
            fig.show_contour(self.filetoplot, colors=contcol, 
                            levels=contlev, linewidths=1.5)
        beampos = figprop.ra_dec_corner(self.filetoplot)
        x_word = beampos[0]
        y_word = beampos[1]
        # Show beam at the lower left corner of the image
        fig.show_ellipses(x_word, y_word, self.bmin, self.bmaj, 
                        angle=self.bpa, edgecolor=beamcolor)
        plt.tight_layout()
        output_figure = self.figure_path + self.figure_name
        plt.savefig(output_figure, bbox_inches="tight", 
                    dpi=fig_prop['general']['dpi'])

    def plotspec(self, linecol='#468499'):
        """
           Plotting integrated spectrum
        """
        fig = plt.figure(figsize=(fig_prop['general']['figsize'][0], 
                                  fig_prop['general']['figsize'][1]))
        ax = fig.add_axes([fig_prop['general']['x0'], 
                           fig_prop['general']['y0'], 
                           fig_prop['general']['w'], 
                           fig_prop['general']['h']])
        fspec = np.loadtxt(self.filetoplot)
        x = fspec[:,1] / 1000
        y = fspec[:,2]
        ax.plot(x, y, '-', color=linecol, 
            lw=fig_prop['line_plot']['line_width'])
        plt.text(0.3, 0.90, self.label, 
                 transform=ax.transAxes, color="black", 
                 fontsize=font_size)
        figprop = figure_properties.Figureprop('', mat_ylabelpad=65)
        figprop.axprop(ax, xlabel=r"$\mathrm{Velocity~(km~s^{-1})}$",
            ylabel=r"$\mathrm{Flux~(Jy)}$", is_aplpy = False)
        output_figure = self.figure_path + self.figure_name
        plt.savefig(output_figure, bbox_inches="tight", 
                dpi=fig_prop['general']['dpi'])

    def plotpvd(self, pv_label, side="", pos1=0.08,
                pos2=0.9,label="", textcol="black", pvdrms_in=None):
        """ Plot PV diagrams
        inputs:
        f: data cube
        hdr: header
        img: 3D array
        figname: name of output file name
        side: A, B or C
        dirc: directory to save files
        """
        f = self.filetoplot.replace("mom0", pv_label) 
        hdu = fits.open(f)[0]
        hdr = hdu.header
        img = hdu.data
        hdr2 = hdr.copy()
        hdr2["CDELT2"] = (hdr["CDELT2"] ) / 1000
        hdr2["CRVAL2"] = (hdr["CRVAL2"] ) / 1000
        hdr2["CDELT1"] = hdr["CDELT1"] * 60
        hdr2["CRVAL1"] = hdr["CRVAL1"] * 60
        f2 = f[:-5]+"_arcsec.fits"
        fits.writeto(f2, img, hdr2, overwrite=True)
        font = "tex gyre heros"
        mpl.rcParams['font.sans-serif'] = font #font.get_name()
        mpl.rc('mathtext', fontset='custom', it=font + ':italic')
        mpl.rc('font', size=30)
        fig = aplpy.FITSFigure(f2, dimensions=[0,1])
        fig.set_theme('publication')
        fig.ticks.set_tick_direction('in')
        fig.axis_labels.set_xpad(2)
        fig.axis_labels.set_ypad(2)
        fig.tick_labels.set_font(size = 30)
        fig.axis_labels.set_font(size = 30)
        colors_noise = plt.cm.gray(np.linspace(0, 1, 256))
        colors_galaxy = plt.cm.afmhot(np.linspace(1, 0.4, 256))
        all_colors = np.vstack((colors_noise, colors_galaxy))
        pvd_map = colors.LinearSegmentedColormap.from_list('pvd_map', all_colors)
        pvd = img
        if pvdrms_in:
            pvd_rms = pvdrms_in
        else:
            pvd_rms = 1.4826 * np.nanmedian(np.abs(pvd[pvd < 0])) # assuming emission is all noise 
        sigma_factor = 3
        divnorm = colors.TwoSlopeNorm(vmin=-sigma_factor*pvd_rms, vcenter=+sigma_factor*pvd_rms, vmax=15*pvd_rms)
        plt.gca().invert_yaxis()
        print(pvd_rms, "PVDRMS")
        ax4 = plt.gca()
        if pv_label == "pv":
            ax4.set_xlim(0.5,62)
        else:
            ax4.set_xlim(2,60)
        ax4.imshow(pvd, cmap=pvd_map, aspect='auto', norm=divnorm)
        if np.nanmax(pvd) > sigma_factor*pvd_rms:
            ax4.contour(pvd, colors=['b', ], levels=sigma_factor**np.arange(1, 10)*pvd_rms)
            # Plot negative contours
        if np.nanmin(pvd) < -sigma_factor*pvd_rms:
            ax4.contour(pvd, colors=['w', ], levels=-pvd_rms * sigma_factor**np.arange(10, 0, -1), linestyles=['dashed', ])
        ax4.set_ylim(ax4.get_ylim()[0], ax4.get_ylim()[1])
        ax4.set_xlim(ax4.get_xlim()[0], ax4.get_xlim()[1])
        #ax4.text(0.5, 0.94, side, color=textcol, fontsize=font_size,transform=ax4.transAxes)
        ax4.text(pos1, pos2, label, fontsize=font_size, transform=ax4.transAxes, 
                color="white", bbox=dict(facecolor="black", edgecolor="none", boxstyle="round,pad=0.3"))
        ax4.invert_yaxis()
        fig2 = plt.gcf()
        fig2.set_figheight(12)
        fig2.set_figwidth(12)
        ax = plt.gca()
        ax.invert_yaxis()
        xmin, xmax = ax.get_xlim()
        print("Current x-axis limits:", xmin, xmax)
        ax.tick_params(direction='in', length=8.7, width=1.3, pad=12)
        ax.tick_params(which='minor', length=5)
        ax.set_xlabel('Offset (arcmin)', labelpad=1.5)

        ax.set_ylabel(r'$\mathrm{Velocity~(km~s^{-1})}$', labelpad=1.5)
        #ax.set_ylabel(r'V [Mhz]', labelpad=1.5)
        output_figure = self.figure_path + self.figure_name
        plt.savefig(output_figure, bbox_inches="tight", dpi=400)


class Plotsofiaproducts:
    def __init__(self, input_args):
        self.input_args = input_args

    def run(self):
        print("Making plots for NGC 1622")
        contours = (self.input_args.filetoplot[:-11] 
                + self.input_args.source_id + '_mom0.fits')
        figure_name = "ngc1622_column_density.pdf"
        plotter = Plotter(figure_name, 
                contours, label='column density')
        base_contour = 5.7E+18
        n = np.arange(0,8)
        input_nhi = base_contour * pow(2, n)
        bmin = plotter.bmin * 3600
        bmaj = plotter.bmaj * 3600
        contlev = nhitomom0(input_nhi, bmin, bmaj)
        print("Plotting column density map", input_nhi)

        rfilter = glob.glob("data/hcg30/optical/ngc1622*-r.fits")[0]

        plotter.plot_momentzero(r_filter=rfilter, 
                contlev=contlev)
        print("Plotting signal-to-noise ratio map")
        plotter.figure_name = "ngc1622_snr.pdf"
        snr_map = (self.input_args.filetoplot[:-11] 
                + self.input_args.source_id + '_snr.fits')
        plotter.filetoplot = snr_map
        plotter.label = "signal-to-noise ratio"
        plotter.plot_sources(colmap=plt.cm.RdYlBu_r, contlev=None, contcol='#2698ff', 
                    beamcolor='black', colorbar_label="Signal-to-noise ratio")
        print("Plotting moment one map")
        plotter.figure_name = "ngc1622_mom1.pdf"
        mom1 = (self.input_args.filetoplot[:-11] 
                + self.input_args.source_id + '_mom1.fits')
        mom1_kms = Fitshandler.dividefits(mom1, mom1[:-5] + '_kms.fits', value=1000) 
        plotter.filetoplot = mom1_kms
        plotter.label = "moment 1"
        plotter.plot_sources(colmap=plt.cm.RdYlBu_r, 
                    contlev=np.arange(4500, 5100, 50),
                    contcol='black', beamcolor='black', 
                    colorbar_label=r'$\rm{Velocity~(km~s^{-1})}$')
        os.system('rm ' + mom1[:-5] + '_kms.fits')
        print("Plotting channel map")
        plotter.figure_name = "ngc1622_chan.pdf"
        cube = (self.input_args.filetoplot[:-11] 
                + self.input_args.source_id + '_cube.fits')
        openf = fits.open(cube)[0]
        img = openf.data
        chan_map = cube[:-5] + '_chanelmap.fits'
        fits.writeto(chan_map, img[97,:,:]*1000, 
                     delheader(openf.header, '3'), overwrite=True)
        plotter.filetoplot = chan_map
        plotter.label = "data cube (channel map)"
        plotter.plot_sources(colmap=plt.cm.RdYlBu_r, 
                    beamcolor='white', 
                    colorbar_label=r'$\rm{Flux~(mJy~beam^{-1})}$', 
                    textlabelcolor='white', text_xpos=0.2)
        os.system('rm ' + chan_map)
        print("Plotting moment 2 map")
        plotter.figure_name = "ngc1622_mom2.pdf"
        mom2 = (self.input_args.filetoplot[:-11] 
                + self.input_args.source_id + '_mom2.fits')
        mom2_kms = Fitshandler.dividefits(mom2, mom2[:-5] + '_kms.fits', value=1000)
        plotter.filetoplot = mom2_kms
        plotter.label = "moment 2"
        plotter.plot_sources(colmap=plt.cm.RdYlBu_r, 
                    beamcolor='black', 
                    colorbar_label=r'$\rm{Velocity~(km~s^{-1})}$')
        os.system('rm ' + mom2_kms)
        print("Ploting integrated spectrum")
        plotter.figure_name = "ngc1622_spec.pdf"
        plotter.label = "Integrated spectrum"
        plotter.filetoplot = (self.input_args.filetoplot[:-11] 
                + self.input_args.source_id + '_spec.txt')
        plotter.plotspec()
        print("Plot position velocity diagram")
        plotter.label = "Integrated spectrum"
        label = ["major axis pv-diagram", "minor axis pv-diagram"]
        pv_label = ["pv", "pv_min"]
        for j in range(2):
            plotter.figure_name = "ngc1622_"+pv_label[j]+".pdf"
            plotter.filetoplot = mom1.replace("mom1", pv_label[j])
            plotter.plotpvd(pv_label=pv_label[j], label=label[j], pos2=0.94) 

def get_args():
    '''This function parses and returns arguments passed in'''
    description = 'Select dataset'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-c', '--filetoplot', dest='filetoplot', 
                        help='Name of file to plot')
    parser.add_argument('-s', '--source_id', dest='source_id', \
                        help='Sofia source id number')
    parser.add_argument('-f', '--catalogue', dest='catalogue', \
                        help='Sofia catalogue')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    #ngc1622_catalog = np.loadtxt(args.catalogue, dtype="str")
    #dec = float(ngc1622_catalog[int(args.source_id) - 1][-8])
    #ra = float(ngc1622_catalog[int(args.source_id) - 1][-9])

    hcg_position = {
            "hcg30": {"ra": "4:36:37", "dec": "-3:11:05"},
    }
    
    band = ['g', 'z', 'i', 'r']
    radec = radec2deg(hcg_position["hcg30"]["ra"], hcg_position["hcg30"]["dec"])
    ra = round(radec["ra"],6)
    dec = round(radec["dec"],6)  
    source_position = {"ngc1622":{"ra": ra, "dec": dec}}
    get_deep_optical(source_position, label="ngc1622", directory="data/hcg30/optical/", size=0.1833, pixscale=0.5)
    app = Plotsofiaproducts(args)
    app.run()
    custom_font.setup_fonts()
    # files_to_remove = [data_i, data_z, data_g]
    # for files in files_to_remove:
    #     if os.path.exists(files):
    #         os.remove(files)
#python workflow/scripts/plot_sofia_products.py --filetoplot outputs/sofia_pbc/hcg30/hcg30_line60_masked.pb_corr_vopt_cubelets/hcg30_line60_masked.pb_corr_vopt_1_mom0.fits --source_id 7  --catalogue outputs/sofia_pbc/hcg30/hcg30_line60_masked.pb_corr_vopt_cat.txt 