import os
import yaml
import glob
import aplpy
import numpy
import astropy.units as u
import argparse
from pvextractor.geometry import Path
import figure_properties
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as units
import matplotlib.pyplot as plt
from pvextractor import extract_pv_slice
from reproject import reproject_interp
from astropy.coordinates import SkyCoord
from utility_functions import get_decals
from PIL import Image
from analysis_tools.functions import delheader
from utility_functions import get_deep_optical
from analysis_tools.functions import FitsCutter
from utility_functions import invert_image_colors
from analysis_tools.functions import radec2deg
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, ListedColormap
from matplotlib.cm import get_cmap
from utility_functions import generate_random_name
from astropy.visualization import (SinhStretch, ImageNormalize)
import cv2
from astropy import wcs

zoom_area = open("config/parameters.yaml")
zoom_rect = yaml.safe_load(zoom_area) # position  and size of zoom area

# moment one contour levels for the group center 
mom1_levels = {
    "hcg16": numpy.arange(3600, 4000, 50), 
    "hcg30": numpy.arange(4500, 5000, 50),
    "hcg31": numpy.arange(3900, 4200, 25),
    "hcg91": numpy.arange(6800, 8000, 50),
    "hcg97": numpy.arange(6000, 7400, 50),
    "hcg90": numpy.arange(2000, 3200, 50)
}

hcg_position = {
            "hcg90": {"ra": 330.542, "dec": -31.987}
            }

def reproject_fits_cube(data_cube, base_header, data_cube_out):
    """
    Reprojects each slice of a 3D FITS cube to match a given base header.

    Parameters:
        data_cube (numpy.ndarray): The 3D numpy array from the original FITS file.
        original_header (astropy.io.fits.Header): The header associated with the data_cube.
        base_header (astropy.io.fits.Header): The header to reproject each slice into.

    Returns:
        numpy.ndarray: A reprojected 3D FITS data cube.
    """
    hdu = fits.open(data_cube)[0]
    original_header = hdu.header
    original_wcs = WCS(original_header)
    base_wcs = WCS(base_header)
    shape_out = (base_header['NAXIS2'], base_header['NAXIS1'])  # Output shape from the base header
    
    # Initialize an empty array to hold the reprojected data
    reprojected_data = numpy.zeros((data_cube.shape[0], shape_out[0], shape_out[1]))
    
    # Reproject each slice
    for i in range(data_cube.shape[0]):
        reprojected_slice, _ = reproject_interp((data_cube[i, :, :], original_wcs.celestial),
                                                output_projection=base_wcs.celestial,
                                                shape_out=shape_out)
        reprojected_data[i, :, :] = reprojected_slice
    return reprojected_data

def pvdiagrams(filename, outname, centra, centdec, 
               pa, majlen, minlen, mindist, majlen_opt, minlen_opt, mindist_opt, dirc, pv_width, base_header, factor=1):
    """
    Produce four PV diagrams from filename

    Input:
    filename (str)     : Input file name
    outname (str)      : Prefix to output file names
    centra (float)     : Right Ascension of the centre in degrees
    centdec (float)    : Declination of the centre in degrees
    pa (float)         : Position angle, N over E, of the receding side major axis in degrees
    majlen (float)     : Length of cut along major axis in arcsec
    minlen (float)     : Length of cut along minor axis and along parallel lines in arcsec
    mindist (float)    : Distance of cuts along parallel lines from cut along minor axis

    Will produce four pv diagrams from a cube filename, one along the
    "major" axis, defined by the central coordinates centra (deg),
    centdec (deg), position angle (pa in degrees), and length majlen
    in arcsec, three perpendicular with respect to the major axis, one
    through the centre, two (at each side) at a distance of mindist
    (in arcsec) from the centre, of the length minlen. Output PV
    diagrams will be called outname+'_pvmaj.fits',
    outname+'_pvmin_cent.fits', outname+'_pvmin_left.fits',
    outname+'_pvmin_right.fits'. Position angle should be defined as
    the receding side major axis.
    """
    print('Creating Pvdiagrams')
    from astropy.wcs.utils import skycoord_to_pixel
    from astropy.coordinates import SkyCoord
    import astropy.units as unit
    fitsfile = fits.open(filename)
    hdu = fitsfile[0]
    print("HDU_DATA_SHAPE", hdu.data.shape)
    hdrc = hdu.header.copy()
            
    w = wcs.WCS(hdrc, fitsfile)

    cenco = numpy.array([[centra,centdec]])

    central_pix_x, central_pix_y = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_world2pix(centra,centdec, 0)

    # Convert pixel coordinates back to celestial coordinates (to be sure we're accurate)
    central_ra_dec = w.celestial.wcs_pix2world(central_pix_x, central_pix_y, 0)

    # Now load the 2D base header WCS
    base_wcs = WCS(base_header).celestial

    # Convert the RA and Dec from the original cube to pixel coordinates in the base image
    base_pix_x, base_pix_y = base_wcs.wcs_world2pix(central_ra_dec[0], central_ra_dec[1], 0)
    print("BASE_PIC", base_pix_x, base_pix_y)

    scale_xy = hdu.header['cdelt2']*3600. # convert arcsec to pixel

    # Make a grid
    endpoints_y_b = numpy.array([majlen*factor/2., -majlen/2.,  mindist  ,  mindist ,         0.,        0.,   -mindist,  -mindist])
    endpoints_x_b = numpy.array([       0.,         0., -minlen/2., minlen/2., -minlen/2., minlen/2., -minlen/2.-0.0001, minlen/2.+0.0001])

    endpoints_y_b_opt = numpy.array([majlen_opt*factor/2., -majlen_opt/2.,  mindist_opt  ,  mindist_opt ,         0.,        0.,   -mindist_opt,  -mindist_opt])
    endpoints_x_b_opt = numpy.array([       0.,         0., -minlen_opt/2., minlen_opt/2., -minlen_opt/2., minlen_opt/2., -minlen_opt/2.-0.0001, minlen_opt/2.+0.0001])

    endpoints_y_bt = numpy.array([majlen*factor/2.+2.5, -majlen/2.+2.5,  mindist  ,  mindist ,         0.,        0.,   -mindist,  -mindist])
    endpoints_x_bt = numpy.array([       0.,         0., -minlen/2.-2.5, minlen/2.-2.5, -minlen/2.-2.5, minlen/2.-2.5, -minlen/2.-0.0001-2.5, minlen/2.+0.0001-2.5])

    # Rotate
    endpoints_x = endpoints_x_b*numpy.cos(numpy.pi*pa/180.)-endpoints_y_b*numpy.sin(numpy.pi*pa/180.)+central_pix_x
    endpoints_y = endpoints_x_b*numpy.sin(numpy.pi*pa/180.)+endpoints_y_b*numpy.cos(numpy.pi*pa/180.)+central_pix_y

    # Rotate
    endpoints_x_opt = endpoints_x_b_opt*numpy.cos(numpy.pi*pa/180.)-endpoints_y_b_opt*numpy.sin(numpy.pi*pa/180.)+base_pix_x
    endpoints_y_opt = endpoints_x_b_opt*numpy.sin(numpy.pi*pa/180.)+endpoints_y_b_opt*numpy.cos(numpy.pi*pa/180.)+base_pix_y

    print("xst = ", tuple(endpoints_x))
    print("yst = ", tuple(endpoints_y))
    endpoints_xt = endpoints_x_bt*numpy.cos(numpy.pi*pa/180.)-endpoints_y_bt*numpy.sin(numpy.pi*pa/180.)+central_pix_x
    endpoints_yt = endpoints_x_bt*numpy.sin(numpy.pi*pa/180.)+endpoints_y_bt*numpy.cos(numpy.pi*pa/180.)+central_pix_y

    i = -1
    list_points = []
    list_points_opt = []
    for names in ['_pvmaj.fits','_pvmin_left.fits', '_pvmin_cent.fits', '_pvmin_right.fits']:
        i = i + 1
        endpoints = [(endpoints_x[2*i],endpoints_y[2*i]),(endpoints_x[2*i+1],endpoints_y[2*i+1])]
        endpoints_opt = [(endpoints_x_opt[2*i],endpoints_y_opt[2*i]),(endpoints_x_opt[2*i+1],endpoints_y_opt[2*i+1])]
        endpointst = [(endpoints_xt[2*i],endpoints_yt[2*i]),(endpoints_xt[2*i+1],endpoints_yt[2*i+1])]
        if names == '_pvmaj.fits':
            print("ENDPOINTS")
            print(endpoints)
        xy = Path(endpoints, width=pv_width*u.arcsec)
        xy_opt = Path(endpoints_opt, width=pv_width*u.arcsec)
        pv = extract_pv_slice(hdu, xy)
        list_points.append(xy) 
        list_points_opt.append(xy) 
        header = pv.header
        #    print header
        pixels =  header['NAXIS1']
        pv.header['CRPIX1'] = pixels/2
        pv.header['CDELT1'] = pv.header['CDELT1']
        pv.writeto(dirc + outname + names, overwrite = True)

    return {"xst" : tuple(endpoints_x), "yst": tuple(endpoints_y), "path": list_points[0], "path_reproj": list_points_opt[0]}

def rect_coord(gal):
    """ gal e.g. hcg16"""
    coord_hcg = SkyCoord(zoom_rect[gal]["zoom_center"][0] * u.deg, 
                         zoom_rect[gal]["zoom_center"][1] * u.deg, frame='icrs')
    center_ra = coord_hcg.ra.deg
    center_dec = coord_hcg.dec.deg
    return {"center_ra": center_ra, "center_dec": center_dec}

fig_prop = open("config/figure_properties.yaml")
fig_prop = yaml.safe_load(fig_prop)

def show_fits(fits_map, output,  gal, bmin, bmaj, bpa, img, pvd_path, pathcol="blue", levels=[1], 
              rectsize=0.5, colmap = plt.cm.gray_r,contcol="black", 
              showcontour='no', showzoomarea='no', text_label='', 
              cbar_label='', optical='no', ellipsecolor='black', 
              rgbcube='', rectanglecolor='black', linecolor='black', 
              showcolorbar='yes', textcolor='black', vmin_vmax_plot='no', 
              vmin=1, vmax=1, addlabel='no', ra=[], dec=[], sofia_id='1',
              cat='', mom1_basename='', png_image=''):
        
    size_deg = rectsize / 60
    fig_size = plt.figure(figsize=(fig_prop['general']['figsize'][0],
                                   fig_prop['general']['figsize'][1]))
    if optical != 'yes':
        fits_fig = aplpy.FITSFigure(fits_map, figure=fig_size, 
                     subplot=[fig_prop['general']['x0'],
                     fig_prop['general']['y0'],
                     fig_prop['general']['w'],
                     fig_prop['general']['h']])

        if vmin_vmax_plot == 'no':
            fits_fig.show_colorscale(vmin=vmin, vmax=vmax, aspect='auto', 
                                  cmap=colmap, stretch='linear')
        else:
            plot_mom1_sources(fits_fig, cat, mom1_basename)

    else:
            rgb2d = hcg90_opt["fits"]
            png_image = invert_image_colors(png_image)
            if gal == "hcg90" and 'zoom' in fits_map:
                png_image_opt = hcg90_opt["png"]
                png_image = "data/hcg90/optical/png_image_hcg90.png"

            else: 
                rgb2d = rgbcube[:-5]+'_2d.fits'

            fits_fig = aplpy.FITSFigure(fits_map, figure=fig_size, 
                                    subplot=[fig_prop['general']['x0'],
                                    fig_prop['general']['y0'],
                                    fig_prop['general']['w'],
                                    fig_prop['general']['h']
                                    ])
                ## SHOW PVD path
            ax3 = fits_fig._figure.axes[0]
            fits_fig.show_contour(fits_map, levels=levels, colors="red", smooth=1)
            # fits_fig.show_rgb(png_image)
            pvd_path.show_on_axis(ax3, spacing=1,
                          edgecolor=pathcol, linestyle='-',
                          linewidth=0.75)
    if addlabel == 'yes':
        for j in range(len(ra)):
            fits_fig.add_label(ra[j], dec[j], sofia_id[j], color='red', 
                               size=fig_prop['general']['axtext_fontsize'])
    figprop = figure_properties.Figureprop(fits_fig)
    figprop.aplpy_figprop()
    fig_size.canvas.draw()
    beampos = figprop.ra_dec_corner(fits_map, edge_offset_percent=5)
    x_word = beampos[0]
    y_word = beampos[1]
    # Show beam at the lower left corner of the image
    fits_fig.show_ellipses(x_word, y_word, bmin, bmaj,
                           angle=bpa, edgecolor=ellipsecolor)    
    ax = plt.gca()
    ###
    figprop.axprop(ax, secaxtickhide=True) 
    ax.text(0.6, 0.90, text_label, transform=ax.transAxes, color=textcolor, 
             fontsize=fig_prop['general']['axtext_fontsize'])
    ax.text(0.05, 0.90, gal.upper()[:3] + ' ' + gal[3:], 
             transform=ax.transAxes, color=textcolor, 
             fontsize=fig_prop['general']['axtext_fontsize'])
    if showcontour == 'yes' and optical != "yes":
        ax.contour(img, levels=levels, colors=contcol)
    if showcolorbar == 'yes':
        figprop.aplpy_colorbar(ylabel=cbar_label)
    if showzoomarea == 'yes':
        # Get the limits of the plot
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        
        # Convert the limits to world coordinates
        plot_top_right = fits_fig.pixel2world(xlim[1], ylim[1])
        plot_bottom_right = fits_fig.pixel2world(xlim[1], ylim[0])

        # Coordinates for the dashed lines
        top_right_ra = rect_coord(gal)["center_ra"] + size_deg / 2
        top_right_dec = rect_coord(gal)["center_dec"] + size_deg / 2
        bottom_left_ra = rect_coord(gal)["center_ra"] + size_deg / 2
        bottom_left_dec = rect_coord(gal)["center_dec"] - size_deg / 2

        fits_fig.show_rectangles(rect_coord(gal)["center_ra"], 
            rect_coord(gal)["center_dec"], size_deg, 
            size_deg, edgecolor=rectanglecolor, facecolor='none', alpha=0.8, 
            lw=2, linestyle='dashed')
       
        fits_fig.show_lines([numpy.array([[top_right_ra, plot_top_right[0]], 
                            [top_right_dec, plot_top_right[1]]]),
                            numpy.array([[bottom_left_ra, 
                            plot_bottom_right[0]], 
                            [bottom_left_dec, plot_bottom_right[1]]])],
                            color=linecolor, linestyle='--', linewidth=2, alpha=0.8) 

    plt.tight_layout()
    plt.savefig(output, bbox_inches="tight", dpi=fig_prop['general']['dpi'])

def get_args():
    '''This function parses and returns arguments passed in'''
    description = 'Select dataset'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-f', '--fitsfile', dest='fitsfile', \
                        help='Fits file to plot')
    parser.add_argument('-c', '--cube', dest='fitscube', \
                        help='Fits cube to extract pv diagram')
    args = parser.parse_args()
    return args

custom_font = figure_properties.FontManager(
        ['/opt/fits-tools/src/analysis_tools/'], 'tex gyre heros')

def cut_fits(f, rect_center, rectsize, f_out=None, f_out3d = None, chanstart = None, chanend = None):
    """ Cut a fits file within a box 
    
    Parameters
    ----------
    f: str
        Name of the fits file
    rectsize: float
        Radius of the rectangle in arcmin
    f_out: str:
        Output fits
    """
    cutter_fits = FitsCutter(f)
    coord_hcg = SkyCoord(rect_center["ra"] * u.deg, 
                         rect_center["dec"] * u.deg, frame='icrs')
    center_ra = coord_hcg.ra.deg
    center_dec = coord_hcg.dec.deg
    if f_out:
        cutter_fits.cut_fits2D(center_ra,
                           center_dec, 
                           rectsize, f_out)
        return f_out
    if f_out3d:
        cutter_fits.cut3D(center_ra, center_dec, 
          rectsize, chanstart, chanend, f_out3d)
        return f_out3d

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

def get_outputfig(mom_type, filetoplot, gal):
    """Determine which file are we plotting, is it moment1, moment0, etc? 
    E.g., if it is moment 1, we need to show contours 
    by setting showcontour to yes 
    
    """
    if mom_type in filetoplot:
        print(f"Plotting moment {mom_type[-1]} map of {gal.upper()}")
        file_suffix = f'_mom{mom_type[-1]}_pbc_zoom.pdf'
        output = f"results/publication/figures/{gal}{file_suffix}"
        output_fig = [output, output.replace('_zoom.pdf', '_large.pdf')]
        showcontour = {'mom1': 'yes', 'mom0': 'no'}  
        mom1_label = '$\mathrm{Velocity~(km~s^{-1})}$' 
        cmap_mom1 = plt.cm.coolwarm
        cmap_mom0 = plt.cm.RdYlBu_r
        mom0_label = '$\mathrm{Flux~(Jy~beam^{-1}~km~s^{-1})}$'
        mom0_textlabel = "Moment 0"  
        mom1_textlabel = "Moment 1"  
        cbar_label = {'mom1': mom1_label, 'mom0': mom0_label} 
        text_label = {'mom1': mom1_textlabel, 'mom0': mom0_textlabel}
        cmaps = {'mom1': cmap_mom1, 'mom0': cmap_mom0}
        return {'outputfig': output_fig, 'showcontour': showcontour[mom_type],
                'cbar_label': cbar_label[mom_type], 'cmap_mom': cmaps[mom_type],
                'text_label': text_label[mom_type]}
 
def determine_mom_type(filename):
    if 'mom0' in filename:
        return 'mom0'
    elif 'mom1' in filename:
        return 'mom1'
    else:
        return None  # 

def calculate_coldens(img, bmin, bmaj):
    # convert moment zero in Jy/beam * km/s to column density
    bmin_arc = bmin * 3600
    bmaj_arc = bmaj * 3600
    nhi = 1.104e+21 * img * 1000 / bmin_arc / bmaj_arc
    return nhi    

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
        hdu.writeto(output_file, overwrite=True)
    return output_file

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
    reprojected_image[numpy.isnan(reprojected_image)] = 0
    hdu = fits.PrimaryHDU(reprojected_image, header=base_header)
    # Save the reprojected image as a new FITS file
    hdu.writeto(output, overwrite=True)
    return output

def signal_to_noise_map(noise_cube, mom0, channels):
    """Derive signal to noise map. 

    Parameters
    ----------
    noise_cube: str
        Cube containing the noise values
    mom0: str
        Moment zero maps 
    channels: str
        Maps containing the number of contributing channel 
    """

    cube = fits.open(noise_cube)[0]
    mom0_values = fits.open(mom0)[0]
    chanwdth_opt = abs(cube.header["CDELT3"])
    median_noise = numpy.nanmedian(cube.data, axis=0)
    channel_number = fits.open(channels)[0]
    snmap = mom0_values.data / (median_noise * numpy.sqrt(channel_number.data) * chanwdth_opt)
    snmap_header = delheader(cube.header)  
    snmap_hdr = delheader(snmap_header, "3")
    fits.writeto(mom0[:-5] + '_snmap.fits',
                 snmap, snmap_hdr, overwrite=True)
    return mom0[:-5] + '_snmap.fits'

def generate_blankfits(f):
    original_fits = fits.open(f)[0]
    
    # Create a new FITS image of the same size filled with nan
    data = numpy.nan * numpy.zeros((original_fits.data.shape[0], original_fits.data.shape[1]))
    header = original_fits.header
    new_fits = fits.PrimaryHDU(data=data, header=header)
    
    fblk = generate_random_name() + "_blank.fits"
    # Save this new FITS file
    new_fits.writeto(fblk, overwrite=True)
    return fblk
        
# use larger mom1 step level for these sofia sources  
large_source = {"hcg16": [1, 50],
                 "hcg30": [7, 50],
                 "hcg31": [3, 25],
                 "hcg91": [8, 50]}

plot_data = {
    "hcg16": {
        "pmin": -0.01427534817185605,
        "pmax": 0.04519556238746736

    },
    "hcg30": {
        "pmin": -0.01427534817185605,
        "pmax": 0.04519556238746736
    },
    "hcg31": {
        "pmin": -0.01427534817185605,
        "pmax": 0.04519556238746736
    },
    "hcg90": {
        "pmin": -0.01427534817185605,
        "pmax": 0.04519556238746736
    },
    "hcg91": {
        "pmin": -0.01427534817185605,
        "pmax": 0.04519556238746736
    },
    "hcg97": {
        "pmin": -0.01427534817185605,
        "pmax": 0.04519556238746736
    }
}

def apply_intensity_colormap(png_path, base_colormap='gist_heat', num_colors=8):
    # Load the image and convert to grayscale
    img = Image.open(png_path).convert('L')
    data = numpy.array(img)
    
    # Normalize the image data
    norm = Normalize(vmin=data.min(), vmax=data.max())
    
    # Get the base colormap
    base_cmap = get_cmap(base_colormap, num_colors)
    colors = base_cmap(numpy.linspace(0, 1, num_colors))
    colors[-1] = numpy.array([1.0, 1.0, 1.0, 1.0])  # Set the first color to white
    
    # Create a new custom colormap
    custom_cmap = ListedColormap(colors)
    
    # Apply the custom colormap
    mappable = ScalarMappable(norm=norm, cmap=custom_cmap)
    colored_data = mappable.to_rgba(data, bytes=True)  # Convert to RGBA bytes
    
    # Create a new image from the colored data
    new_img = Image.fromarray(colored_data, 'RGBA')
    new_img_path = png_path[:-4] + '_custom_colored.png'
    new_img.save(new_img_path)
    return new_img_path

if __name__ == '__main__': 
    args = get_args()
    custom_font.setup_fonts()
    f = args.fitsfile
    cube = args.fitscube
    gal = os.path.basename(f)[:5]
    openf = fits.open(f)[0]
    snmap = signal_to_noise_map(f[:-9] + "noise.fits" , 
                                f[:-9] + "mom0.fits", f[:-9] + "chan.fits")
    img = openf.data / 1000
    hdr = openf.header
    bmin = hdr["BMIN"]
    bmaj = hdr["BMAJ"]
    bpa = hdr["BPA"]
    img[img==0] = numpy.nan
    f_nan = f[:-5]+'_nan.fits'
    hdu_nan = fits.PrimaryHDU(data=img, header=hdr)
    hdu_nan.writeto(f_nan, overwrite=True)
    fplot = f_nan
    zoom_size = 14.5
    pvpos = radec2deg('22:02:10.0550438695', '-32:01:09.3611094458')
    cut_pos = radec2deg('22:02:10.0550438695', '-32:01:0.3611094458')
    fplot_large = cut_fits(fplot, cut_pos, zoom_rect[gal]["zoom_size"][1], f_out=fplot[:-5]+'_large.fits')
    fplot_zoom = cut_fits(fplot, cut_pos, zoom_size, f_out=fplot[:-5]+'_zoom.fits')
    snmap_large = cut_fits(snmap, cut_pos, zoom_rect[gal]["zoom_size"][1], f_out=fplot[:-5]+'_snmap_large.fits')
    snmap_zoom = cut_fits(snmap, cut_pos, zoom_size, f_out=fplot[:-5]+'_snmap_zoom.fits')
    snmap_large_img = fits.open(snmap_large)[0].data
    snmap_zoom_img = fits.open(snmap_zoom)[0].data

    mom_type = determine_mom_type(fplot)  
    mom_figs = get_outputfig(mom_type, fplot, gal)
    output_fig = mom_figs['outputfig']
    # whether to show contour or not
    showcontour = mom_figs['showcontour']
    # color bar label
    cbar_label = mom_figs['cbar_label']
    # matplotlib ax.text label
    text_label = mom_figs['text_label']
    fig_cmap = mom_figs['cmap_mom']
    levels = mom1_levels[gal]
    showzoomarea = ['no', 'yes']
    larger_fits = fits.open(fplot_large)[0]
    zoomed_fits = fits.open(fplot_zoom)[0]
    img_zoom = zoomed_fits.data
    hdr_zoom = zoomed_fits.header
    img_large = larger_fits.data
    hdr_large = larger_fits.header
    data_img = [img_zoom, img_large]
    
    sofia_catalogue = f[:-9] + "cat.txt"
    mom1_basename = f[:-9] + "cubelets/" + os.path.basename(f)[:-9] 
    f_sofia = numpy.loadtxt(sofia_catalogue, dtype=str)
    ra_list = []
    dec_list = []
    sofia_id_list = []
    for i in range(len(f_sofia)):
        ra, dec = float(f_sofia[i][38]), float(f_sofia[i][39])
        sofia_id = f_sofia[i][2]
        ra_list.append(ra)
        dec_list.append(dec)
        sofia_id_list.append(sofia_id)

    # prepare optical data
    folder = os.path.join("data", gal, "optical/") 
    g_filter = glob.glob(folder + gal+'legacystamps_*-gfilter.fits')[0]
    i_filter = glob.glob(folder + gal+'legacystamps_*-ifilter.fits')[0] 
    z_filter = glob.glob(folder + gal+ 'legacystamps_*-ifilter.fits')[0]
    g_filter = glob.glob(folder + gal+'legacystamps_*-gfilter.fits')[0]
    data_g = normalize_and_save_fits(g_filter, generate_random_name() + '.fits')
    data_i = normalize_and_save_fits(i_filter, generate_random_name() + '.fits')    
    data_z = normalize_and_save_fits(z_filter, generate_random_name() + '.fits')
    # Cut optical data withn a rectangular box
    data_g_large = cut_fits(data_g, cut_pos, zoom_rect[gal]["zoom_size"][1],
                            generate_random_name() + '_g_large.fits')
    data_g_zoom = cut_fits(data_g, cut_pos, zoom_size, 
                           generate_random_name() + '_g_zoom.fits')         
    data_i_large = cut_fits(data_i, cut_pos, zoom_rect[gal]["zoom_size"][1],
                            generate_random_name() + '_i_large.fits')
    data_i_zoom = cut_fits(data_i, cut_pos, zoom_size, 
                           generate_random_name() + '_i_zoom.fits')         
    data_z_large = cut_fits(data_z, cut_pos, zoom_rect[gal]["zoom_size"][1],
                            generate_random_name() + '_z_large.fits')
    data_z_zoom = cut_fits(data_z, cut_pos, zoom_size, 
                           generate_random_name() + '_z_zoom.fits')

    rgb_cube_large = generate_random_name() + '_rgb_cube_large.fits'
    rgb_cube_zoom = generate_random_name() + '_rgb_cube_zoom.fits'
    rgb_cube = [rgb_cube_zoom, rgb_cube_large]
    aplpy.make_rgb_cube([data_z_large, data_i_large, data_g_large], rgb_cube_large)
    aplpy.make_rgb_cube([data_z_zoom, data_i_zoom, data_g_zoom], rgb_cube_zoom)
    if gal == "hcg90":
        pathdir = "results/publication/figures/"
        hcg90_opt = get_deep_optical(hcg_position)
        base_header_zoom = fits.open(hcg90_opt["fits"])[0].header 
    else: 
        base_header_zoom = fits.open(rgb_cube_zoom[:-5] + "_2d.fits")[0].header
    base_header_large = fits.open(rgb_cube_large[:-5] + "_2d.fits")[0].header
    ## Calculate pv diagrams
    cube_cut = cut_fits(cube, cut_pos, rectsize = zoom_size, f_out3d = cube[:-5]+"_tail_cube.fits", chanstart = 314, chanend = 511)
    pv_angle = 320
    pv_width = 49.7403671758 # in arcsec
    mindistance = 10 # in arcsec
    pixsize = hdr_zoom["CDELT1"] * 3600 # in arcsec
    pixsize_opt = base_header_zoom["CDELT1"] * 3600
    print("PIXSIZE_OPT", pixsize_opt)
    majoraxis = 12.1503041262 * 60 # in arcsec
    majlen = majoraxis / pixsize
    minlen = majoraxis / 3 /  pixsize
    mindist = mindistance / pixsize 

    majlen_opt = majoraxis / pixsize
    minlen_opt = majoraxis / 3 / pixsize
    mindist_opt = mindistance / pixsize

    pvdiag = pvdiagrams(filename=cube_cut, outname="pvd", 
                    centra=pvpos["ra"], centdec=pvpos["dec"], 
                    pa=pv_angle, majlen=majlen, minlen=minlen, 
                    mindist=mindist, majlen_opt=majlen_opt, minlen_opt=minlen_opt, 
                    mindist_opt=mindist_opt, dirc="results/sofia_pbc/hcg90/", 
                    pv_width=pv_width, base_header=base_header_zoom, factor=1)
    if mom_type == 'mom0':
        print("Plotting column density map of {0}".format(gal))

        max_contour = 30E+19
        # calculate column density from moment zero value
        img_large_nhi = calculate_coldens(img_large, bmin, bmaj)
        #img_large_nhi[numpy.isnan(img_large_nhi)] = 0
        img_zoom_nhi = calculate_coldens(img_zoom, bmin, bmaj)
        #img_zoom_nhi[numpy.isnan(img_zoom_nhi)] = 0
        #img_zoom_nhi[snmap_zoom_img < 3] = 0
        #img_large_nhi[snmap_large_img < 3] = 0
        # define lowest contour according to signal-to-noise map: 
        mask = (snmap_large_img >= 2.5) & (snmap_large_img <= 3.5)
        filtered_img_large_nhi = numpy.where(mask, img_large_nhi, numpy.nan)
        base_contour = numpy.nanmedian(filtered_img_large_nhi)
        img_zoom_nhi[numpy.isnan(img_zoom_nhi)] = 0
        img_large_nhi[numpy.isnan(img_large_nhi)] = 0
        data_img_nhi = [img_zoom_nhi, img_large_nhi]
        hdrs = [base_header_zoom, base_header_large]
        opt_view = numpy.array([60,]) * units.arcmin
        ra_deg = hdrs[0]["CRVAL1"]
        dec_deg = hdrs[0]["CRVAL2"]
        hi_pos = SkyCoord(ra_deg, dec_deg, unit="deg")
        png_image = generate_random_name()+"_png.png"
        image, header = get_decals(hi_pos, opt_view=opt_view, dev_dr=True)
        image.save(png_image)
        # fits_large_nhi = reproject_fits(img_large_nhi, 
        #                                 hdr_large, base_header_large, 
        #                                 generate_random_name() + '_large_nhi.fits') 
        # fits_zoom_nhi = reproject_fits(img_zoom_nhi, 
        #                                hdr_zoom, base_header_zoom, 
        #                                generate_random_name() + '_zoom_nhi.fits')
        fits.writeto("image_large_nhi.fits", img_large_nhi, hdr_large, overwrite=True)
        fits.writeto("image_zoom_nhi.fits", img_zoom_nhi, hdr_zoom, overwrite=True)
        fits_large_nhi = 'image_large_nhi.fits'
        fits_zoom_nhi = 'image_zoom_nhi.fits'
        fits_nhi = [fits_zoom_nhi, fits_large_nhi]
        rgb_cubes = [rgb_cube_zoom, rgb_cube_large]
        contour_levels = base_contour * 2 ** numpy.arange(20) 
        contour_levels = contour_levels[contour_levels <= max_contour]
        print(contour_levels, " CONTOUR LEVELS of {0}".format(gal.upper()))
        output_fig_nhi = ["results/publication/figures/" + gal + "_coldens_pbc_zoom_deep_view_slice.pdf",
                          "results/publication/figures/" + gal + "_coldens_pbc_large_deep_view_slice.pdf"]
        # Plot column density maps for large and zoomed area. 
        for k in range(2):
            png_image = generate_random_name()+'delete.png'
            adjust_r = 1.5
            adjust_g = 1
            adjust_b = 0.3  # Reduce blue channel influence
            aplpy.make_rgb_image(rgb_cubes[k], png_image,
                vmin_g = plot_data[gal]["pmin"] * adjust_g,
                vmax_g = plot_data[gal]["pmax"] * adjust_g,
                vmin_r = plot_data[gal]["pmin"] * adjust_r,
                vmax_r = plot_data[gal]["pmax"] * adjust_r,
                vmin_b = plot_data[gal]["pmin"] * adjust_b,
                vmax_b = plot_data[gal]["pmax"] * adjust_b,
            )
            show_fits(
                      fits_nhi[k],
                      gal=gal, levels=contour_levels,
                      img=data_img_nhi[k],
                      pvd_path = pvdiag["path"],
                      showzoomarea=showzoomarea[k],  
                      showcontour='yes',
                      rectsize=zoom_size,
                      bmin=bmin, bmaj=bmaj, bpa=bpa,
                      output=output_fig_nhi[k], text_label="Column density",
                      colmap=fig_cmap, cbar_label=cbar_label, textcolor='white',
                      optical='yes', rgbcube=rgb_cubes[k], rectanglecolor='white', 
                      contcol="#2698ff", ellipsecolor='white', linecolor='white', 
                      showcolorbar='no', png_image=png_image)
#python workflow/scripts/plot_moments_hcg90_view_slice.py --fitsfile results/sofia_pbc/hcg90/hcg90_line60_masked.pb_corr_vopt_mom0.fits -c results/data/hcg90/hcg90_line60_masked.pb_corr_vopt.fits
