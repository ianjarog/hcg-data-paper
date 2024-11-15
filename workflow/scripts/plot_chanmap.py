import os
import yaml
import glob
import aplpy
import numpy
import string
import random
import argparse
import pyspeckit
from PIL import Image
import figure_properties
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as units
import matplotlib.pyplot as plt
from reproject import reproject_interp
from astropy.coordinates import SkyCoord
from utility_functions import get_decals
from analysis_tools.functions import delheader
from analysis_tools.functions import radec2deg
from utility_functions import get_deep_optical
from PIL import Image, ImageEnhance
from utility_functions import invert_image_colors
from utility_functions import extract_medians_from_fits
from astropy.visualization import (SinhStretch, ImageNormalize)
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, ListedColormap
from matplotlib.cm import get_cmap

hcgs_labels = {"hcg16":{"A":(32.35312,-10.13622), "B":(32.3361,-10.13309), "C":(32.41071,-10.14637), "D":(32.42872,-10.18394), "NGC848":(32.5735,-10.32145)},
              "hcg31":{"G":(75.43338,-4.28875), "Q":(75.40974,-4.22245), "A":(75.41146,-4.25946), "B":(75.39756,-4.26401), "C":(75.40751,-4.2577)},
              "hcg30":{"A":(69.07744,-2.83129), "B":(69.12619,-2.86656), "C":(69.097,-2.79985), "D":(69.15276,-2.84302)},
              "hcg90":{"A":(330.50889,-31.87014), "B":(330.53631,-31.99068), "C":(330.51423,-31.97451), "D":(330.52602,-31.99423)},
              "hcg91":{"A":(332.28174,-27.80984), "B":(332.3183,-27.73134), "C":(332.30881,-27.78241), "D":(332.28567,-27.80086), 
                  "LEDA749936":(332.22916666666663, -27.783333333333335)},
              "hcg97":{"A":(356.84591,-2.30096), "B":(356.90753,-2.31716), "C":(356.84884,-2.35143), "D":(356.82867,-2.31329), "E":(356.83253,-2.28102)}}

coord_hcg30c = {"ra":6.910151635171060e+01, "dec":-2.802908880849139e+00} #SoFiA source 5
coord_hcg30a = {"ra":6.907887052453162e+01, "dec":-2.845510123245648e+00} #SoFiA source 14 

zoom_area = open("config/parameters.yaml")
zoom_center = yaml.safe_load(zoom_area) # position of zoom area
fig_prop = open("config/figure_properties.yaml")
fig_prop = yaml.safe_load(fig_prop)

class FitsProcessor:
    def __init__(self, file1, file2, output_file):
        self.file1 = file1
        self.file2 = file2
        self.output_file = output_file

    def read_fits_data(self, fits_file):
        """Read data from a FITS file."""
        with fits.open(fits_file) as hdul:
            data = hdul[0].data
            header = hdul[0].header
        return data, header

    def add_fits_data(self, data1, data2):
        """Add data from two FITS files."""
        return (data1 + data2) / 2

    def write_fits_data(self, data, header):
        """Write data to a new FITS file."""
        hdu = fits.PrimaryHDU(data, header=header)
        hdul = fits.HDUList([hdu])
        hdul.writeto(self.output_file, overwrite=True)

    def process(self):
        """Read, add, and write FITS data."""
        data1, header1 = self.read_fits_data(self.file1)
        data2, _ = self.read_fits_data(self.file2)
        result_data = self.add_fits_data(data1, data2)
        self.write_fits_data(result_data, header1)

def hidelab(fits_fig, hidelabx, hidelaby):
    if (not hidelaby) and (not hidelabx):
        fits_fig.axis_labels.hide_x()
        fits_fig.axis_labels.hide_y()
    elif (not hidelabx):
        fits_fig.axis_labels.hide_x()
    elif (not hidelaby):
        fits_fig.axis_labels.hide_y()
    else: 
        pass

def show_fits(optical_map, contour_map, output,  gal, bmin, bmaj, bpa, img, levels=[1], 
              contcol="black", showcontour='no', text_label='', 
              optical='no', ellipsecolor='black', xspacing=0,
            textcolor='black', png_image='', vmin=0, vmax=0,
              show_yticks=True, show_xticks=True, rectsize=0.5):
        
    fig_size = plt.figure(figsize=(fig_prop['general']['figsize'][0],
                                   fig_prop['general']['figsize'][1]))
    size_deg = rectsize / 60
    fits_fig = aplpy.FITSFigure(optical_map, figure=fig_size, 
                                subplot=[fig_prop['general']['x0'],
                                         fig_prop['general']['y0'],
                                         fig_prop['general']['w'],
                                         fig_prop['general']['h']
                                        ])  
    
    fits_fig.show_grayscale(invert=True, stretch='linear', vmin=vmin, vmax=vmax)
    fits_fig.show_contour(contour_map, levels=levels, colors=contcol, smooth=1, linewidths=2)
    recenter_coord = plot_data[gal]["recenter"]
    if recenter_coord:
        recenter_coord_deg = radec2deg(ra=recenter_coord[0], dec=recenter_coord[1])
        fits_fig.recenter(recenter_coord_deg["ra"], recenter_coord_deg["dec"], width=size_deg, height=size_deg)
    if gal== 'hcg90': 
        png_image = "data/hcg90/optical/png_image_hcg90.png"
        img = Image.open(png_image)
        image = numpy.array(img)[::-1]
        ax3 = fits_fig._figure.axes[0]
        m = numpy.mean(image)
        s = numpy.std(image)
        ax3.imshow(image,vmin=m-s,vmax=m+30*s,cmap='gist_gray')
    figprop = figure_properties.Figureprop(fits_fig)
    figprop.aplpy_figprop()
    #if hide:
    #    print("HIDE IS TRUE")
    fig_size.canvas.draw()
    beampos = figprop.ra_dec_corner(optical_map, edge_offset_percent=5)
    x_word = beampos[0]
    y_word = beampos[1]
    if gal=="hcg30":
        x_word, y_word = radec2deg("4:36:40", "-2:52:51.5")["ra"], radec2deg("4:36:40", "-2:52:51.5")["dec"]
    if gal=="hcg16":
        x_word, y_word = radec2deg("2:10:44", "-10:26:00")["ra"], radec2deg("2:10:44", "-10:26:00")["dec"]
    members_label = list(hcgs_labels[gal].keys())
    offsets = {"hcg16":{"A":(0,1/60.),"B":(0,1/60.), "C":(0,1/60.), "D":(0,1/60.), "NGC848":(-5/60,0.)}, 
               "hcg31":{"Q":(0.5/60,0), "A":(0.5/60,0), "B":(-0.5/60,0), "C":(0,0.5/60), "G":(0.5/60.,0)}, 
               "hcg30":{"A":(1/60,0), "B":(-1/60,0), "C":(+0.6/60,0), "D":(-0.4/60,0)}, 
               "hcg91":{"A":(-1.2/60, 0), "B":(+0.5/60, 0), "C":(+0.7/60,0), "D":(0,0.5/60), "LEDA749936":(0,-0.5/60)},
               "hcg90":{"A":(0,3/60), "B":(0.5/60,0), "C":(-0.5/60,0), "D":(0,-3/60)}, 
               "hcg97":{"A":(0,0), "B":(0,0), "C":(0,0), "D":(0,0), "E":(0,0)}}
    for lab in range(len(members_label)):
        if members_label[lab]=="NGC848":
            textlab = "NGC 848"
        elif members_label[lab]=="LEDA749936":
            textlab = "LEDA 749936"
        else:
            textlab = members_label[lab]
        offset_ra = offsets[gal][members_label[lab]][0]
        offset_dec = offsets[gal][members_label[lab]][1]
        fits_fig.add_label(hcgs_labels[gal][members_label[lab]][0]+offset_ra, hcgs_labels[gal][members_label[lab]][1]+offset_dec, 
                textlab, color="black", weight="bold", size=30)
    text_labels = radec2deg("2:10:44", "-10:26:00")["ra"], radec2deg("2:10:44", "-10:26:00")["dec"]
    # Show beam at the lower left corner of the image
    fits_fig.show_ellipses(x_word, y_word, bmin, bmaj,
                           angle=bpa, edgecolor=ellipsecolor)    
    ax = plt.gca()
    figprop.axprop(ax, secaxtickhide=True, show_yticks=show_yticks, show_xticks=show_xticks) 
    hidelab(fits_fig, hidelaby=show_yticks, hidelabx=show_xticks)
    def custom_tick_formatter(x, pos):
        if pos == 0:  # First tick label
            return ""
        else:
            return f"{x:.1f}"
    
    # Apply custom formatter to the y-axis
    #ax.coords[0].set_major_formatter(plt.FuncFormatter(custom_tick_formatter))
    if xspacing > 0:
        fits_fig.ticks.set_xspacing(xspacing)
    ax.text(0.6, 0.90, text_label, transform=ax.transAxes, 
             fontsize=fig_prop['general']['axtext_fontsize'],
             color="white", bbox=dict(facecolor="black", edgecolor="none", boxstyle="round,pad=0.3"))
    ax.text(0.05, 0.90, gal.upper()[:3] + ' ' + gal[3:], 
             transform=ax.transAxes, 
             fontsize=fig_prop['general']['axtext_fontsize'],
             color="white", bbox=dict(facecolor="black", edgecolor="none", boxstyle="round,pad=0.3"))
    if showcontour == 'yes' and optical != "yes":
        ax.contour(img, levels=levels, colors=contcol)
    plt.tight_layout()
    plt.savefig(output, bbox_inches="tight")# dpi=fig_prop['general']['dpi'])

def get_args():
    '''This function parses and returns arguments passed in'''
    description = 'Select dataset'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-f', '--fitsfile', dest='fitsfile', \
                        help='Fits file to plot')
    parser.add_argument('-n', '--noisecube', dest='noisecube', \
                        help='Fits file to plot')
    args = parser.parse_args()
    return args

custom_font = figure_properties.FontManager(
        ['/opt/fits-tools/src/analysis_tools/'], 'tex gyre heros')

def rect_coord(gal):
    """ gal e.g. hcg16"""
    coord_hcg = SkyCoord(zoom_center[gal]["zoom_center"][0],
                         zoom_center[gal]["zoom_center"][1], frame='icrs')
    center_ra = coord_hcg.ra.deg
    center_dec = coord_hcg.dec.deg
    return {"center_ra": center_ra, "center_dec": center_dec}

def remove_file(f):
    """Removing a file from the disk 
    
    Parameters
    ---------
    f: str
        File to remove from the disk
    """
    if os.path.exists(f):
        os.remove(f)
        print("Removed {0} from disk".format(f))

def normalize_and_save_fits(input_file, output_file, 
                            stretch=SinhStretch(), fill_value=0):
    with fits.open(input_file) as hdul:
        data = numpy.squeeze(hdul[0].data)
        header = hdul[0].header

        # Check and handle NaN or infinite values
        if numpy.any(numpy.isnan(data)) or numpy.any(numpy.isinf(data)):
            data = numpy.nan_to_num(data, nan=fill_value, 
                                    posinf=fill_value, neginf=fill_value)

        norm = ImageNormalize(data, stretch=stretch)
        normalized_data = norm(data)

        # Handle MaskedArray case
        if isinstance(normalized_data, numpy.ma.MaskedArray):
            normalized_data = normalized_data.filled(fill_value)

        hdu = fits.PrimaryHDU(data=normalized_data, header=header)
        fits.writeto(output_file, normalized_data, header, overwrite=True)
    return output_file

def generate_random_name(length=5):
    """Generate temporary file name that will be deleted later on"""
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(length)) + 'delete'

def get_percentile_vmin_vmax(file_name, lower_percentile, upper_percentile):
    """
    Calculate and return the lower and upper percentile values of the data in a FITS file.

    This function opens a FITS file, extracts the data from its first HDU, and computes
    the lower and upper percentile values specified by the user. It handles the data by
    flattening the array and removing NaN values before performing the percentile calculation.

    Parameters
    ----------
    file_name : str
        The path to the FITS file from which data will be read.
    lower_percentile : float
        The lower percentile value for the data range. It should be a number between 0 and 100.
    upper_percentile : float
        The upper percentile value for the data range. It should be a number between 0 and 100,
        and greater than the `lower_percentile`.

    Returns
    -------
    vmin : float
        The value at the lower percentile in the data.
    vmax : float
        The value at the upper percentile in the data.

    Notes
    -----
    This function is primarily used for determining suitable scaling limits (vmin, vmax)
    for visualization of FITS image data.
    """
    with fits.open(file_name) as hdul:
        data = hdul[0].data
        # Flatten the array and remove NaN values for percentile calculation
        data_flat = data.flatten()
        data_flat = data_flat[~numpy.isnan(data_flat)]
        vmin = numpy.percentile(data_flat, lower_percentile)
        vmax = numpy.percentile(data_flat, upper_percentile)

        return vmin, vmax

def reproject_fits(image_data, image_header, base_header, output):
    """Reproject data array from fits to a new header and save it to a new fits file"""
    reprojected_image, footprint = reproject_interp((image_data, WCS(image_header)), base_header)
    reprojected_image[numpy.isnan(reprojected_image)] = numpy.nan
    fits.writeto(output, reprojected_image, base_header, overwrite=True)
    print(repr(base_header))
    return output

plot_data = {
    "hcg16": {
        #"channels": [171, 172, 173, 180, 181, 182], # 183, 200, 201, 202, 203, 204],
        "channels": [171, 172, 173, 180, 181, 182, 183, 200, 201, 202, 203, 204],
        # "channels": numpy.arange(3500,4300)
        "contcol" : ['blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red'],
        "contour_levels" : numpy.array([1.5, 2, 2.5, 3,  4, 8, 16, 32]), 
        "textcol" : "black",
        "xspacing": 0.16, 
        "recenter": ("02h09m52.5s", "-10d13m20s"),
        "pmin": 1,
        "pmax": 99,
        "pixscale" : 2.5,
    },
    "hcg30": {
        "channels": [153, 164],
        "contcol" : ['blue', 'blue', 'red', 'red', 'red', 'red'],
        "contour_levels" : numpy.array([1.5, 2, 3, 4, 5, 6]), 
        "textcol" : "black",
        "xspacing": 0.03, 
        "recenter": ("4:36:27.99", "-2:49:51.8"),
        "pmin": -0.01427534817185605,
        "pmax": 0.04519556238746736,
        "pixscale": 0.262

    },
    "hcg31": {
        "channels": [334, 345],
        "contcol" : ['blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red'],
        "textcol" : "black",
        "contour_levels" : numpy.array([1.5, 2, 2.5, 3, 6, 9, 16, 32]) , 
        "xspacing": 0.04, 
        "recenter": False,
        "pmin": 1,
        "pmax": 99,
        "pixscale": 0.262

    },
    "hcg90": {
        "channels": [476, 487],
        "contour_levels" : numpy.array([1.5, 2, 3, 4, 5, 6]) ,
        "contcol" : ['yellow', '#00BFFF', 'red', 'red', 'red', 'red'],
        "textcol": "white",
        "xspacing": 0.25,
        "recenter": False,
        "pmin": 0.7,
        "pmax": 99.62,
        "pixscale": 0.262
    },
    "hcg91": {
        "channels": [647, 658],
        "contcol" : ['blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red'],
        "contour_levels" : numpy.array([1.5, 2,  2.5, 3, 6, 9, 16, 32]), 
        "textcol" : "black",
        "xspacing": 0.055,
        "recenter": False,
        "pmin": 1,
        "pmax": 99,
        "pixscale": 0.262

    },
    "hcg97": {
        "channels": [143, 316],
        "xspacing": 0,
        "recenter": False,
        "pmin": 1,
        "pmax": 99,
        "pixscale": 0.262
    }
}

if __name__ == '__main__': 
    args = get_args()
    custom_font.setup_fonts()
    f = args.fitsfile
    gal = os.path.basename(f)[:5]
    openf = fits.open(f)[0]
    img_data = numpy.squeeze(openf.data)
    header_data = delheader(openf.header)
    crval3 = header_data["CRVAL3"]
    cdelt3 = header_data["CDELT3"]
    velocities = []
    chan = []
    for vel in range(img_data.shape[0]):
        print("vel", vel)
        velo = crval3
        velocities.append(velo)
        chan.append(vel)
        crval3 = crval3 + cdelt3
    # print(velocities, "VELOCITIES")
    chanmin = 3500000
    chanmax = 4300000
    filtered_chan = [c for v, c in zip(velocities, chan) if chanmin <= v <= chanmax]
    filtered_velocities = [v for v in velocities if chanmin <= v <= chanmax]
    noise = extract_medians_from_fits(args.noisecube)
    velocities = pyspeckit.Cube(f).xarr.value / 1000
    med_noise = float(numpy.median(noise)) * 1000
    noise_val = {gal: med_noise}
    print("MEDIAN NOISE VALUE OF ", gal.upper(), med_noise)
    with open(f[:-5] + '_median_noise.yaml', 'w') as val:
        yaml.dump(noise_val, val, default_flow_style=False)
    if len(plot_data[gal]["channels"]) == 2:
        chanrange = numpy.arange(plot_data[gal]["channels"][0], plot_data[gal]["channels"][1]+1)
    else:
        chanrange = plot_data[gal]["channels"]
    # chanrange = filtered_chan
    show_yticks = [True, False, False]*(len(chanrange)+1) # hide y ticks 
    show_xticks = [False, False, False]*int((len(chanrange)+1)/3) # hide y ticks 
    show_xticks[-1] = True
    show_xticks[-2] = True
    show_xticks[-3] = True
        # prepare optical data
    folder = os.path.join("data", gal, "optical/") 
    g_filter = glob.glob(folder + 'center_legacystamps_*-gfilter.fits')[0]
    i_filter = glob.glob(folder + 'center_legacystamps_*-ifilter.fits')[0]
    r_filter = glob.glob(folder + 'center_legacystamps_*-rfilter.fits')[0]  
    if gal == 'hcg91':
        optical_fits = i_filter
        print("OPTICAL FITS")
    else:
        optical_fits = r_filter
    with fits.open(optical_fits) as hdul:
        data = hdul[0].data  # Access the image data

    # Calculate vmin and vmax using the 5th and 95th percentiles
    vmin = numpy.percentile(data, 1)   # 5th percentile
    vmax = numpy.percentile(data, 99)  # 95th percentile    
    if gal == 'hcg91':
        z_filter = glob.glob(folder + 'center_legacystamps_*-ifilter.fits')[0] # no z filter for hcg 91
    else:    
        z_filter = glob.glob(folder + 'center_legacystamps_*-zfilter.fits')[0]
    if gal == 'hcg30':
        g_filter = glob.glob(folder + 'center_legacystamps_*-gfilter.fits')[0]
        i_filter = glob.glob(folder + 'center_legacystamps_*-ifilter.fits')[0]
        r_filter = glob.glob(folder + 'center_legacystamps_*-rfilter.fits')[0]  
        z_filter = glob.glob(folder + 'center_legacystamps_*-rfilter.fits')[0]  

    ra_deg = openf.header["CRVAL1"]
    dec_deg = openf.header["CRVAL2"]
    hi_pos = SkyCoord(ra_deg, dec_deg, unit="deg")
    zoom_area = open("config/parameters.yaml")
    zoom_rect = yaml.safe_load(zoom_area) # position  and size of zoom area
    opt_view = numpy.array([zoom_rect[gal]["zoom_size"][0]+60,]) * units.arcmin
    image, header_decals = get_decals(hi_pos, opt_view=opt_view, dev_dr=True, pixscale=plot_data[gal]["pixscale"])
    if gal=="hcg90":
        hcg_position = {
            "hcg90": {"ra": 330.542, "dec": -31.987}
            }
        hcg90_opt = get_deep_optical(hcg_position)
        print(hcg90_opt, "HCG90-OPT")
        header_decals = fits.open(hcg90_opt["fits"])[0].header
    rms_cube = med_noise / 1000
    print(rms_cube, " RMS_CUBE")
    t = 0
    while t < len(chanrange):
        j = chanrange[t]
        print ("Plotting channel map number ", j, "for ", gal)
        img = numpy.squeeze(openf.data)[j, :, :]
        vel_j = round(velocities[j])
        hdr = delheader(openf.header)
        bmin = hdr["BMIN"]
        bmaj = hdr["BMAJ"]
        bpa = hdr["BPA"]
        f_zoom = generate_random_name() + str(j) + 'cube_zoom.fits'
        fits.writeto(f_zoom, img, delheader(hdr, "3"), overwrite=True) 
        zoomed_fits = fits.open(f_zoom)[0]
        img_zoom = zoomed_fits.data
        hdr_zoom = zoomed_fits.header
        
        # calculate column density from moment zero value
        fits_zoom = reproject_fits(img_zoom, hdr_zoom, header_decals, 
                                generate_random_name() + '_zoom_chan.fits')
        output_fig = "outputs/publication/figures/" + gal + "_chanmap" + str(j) +"_pbc_zoom.pdf"
        # Plot column density maps for large and zoomed area. 
        show_fits(
                  optical_fits,f_zoom,
                  gal=gal, levels=plot_data[gal]["contour_levels"]*rms_cube,
                  rectsize=zoom_rect[gal]["zoom_size"][0],
                  img=img_zoom, showcontour='yes', png_image=None,
                  bmin=bmin, bmaj=bmaj, bpa=bpa, xspacing=plot_data[gal]["xspacing"],
                  output=output_fig, text_label="V = " + str(vel_j) + "$~\mathrm{km~s^{-1}}$",
                  textcolor=plot_data[gal]["textcol"], optical='yes',vmin=vmin,vmax=vmax,
                  contcol=plot_data[gal]["contcol"], ellipsecolor='black',show_yticks=show_yticks[t], 
                  show_xticks=show_xticks[t]  
                 )
        t = t + 1
#python workflow/scripts/plot_chanmap.py -f outputs/data/hcg30/hcg30_line60_masked.pb_corr_vopt.fits -n outputs/sofia_pbc/hcg30/hcg30_line60_masked.pb_corr_vopt_noise.fits
#python workflow/scripts/plot_chanmap.py -f outputs/data/hcg16/hcg16_line60_masked.pb_corr_vopt.fits -n outputs/sofia_pbc/hcg16/hcg16_line60_masked.pb_corr_vopt_noise.fits
