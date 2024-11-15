import os
import yaml
import glob
import aplpy
import numpy
import string
import random
import argparse
from PIL import Image
import figure_properties
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as units
import matplotlib.pyplot as plt
from legacystamps import download
from matplotlib.patches import Ellipse
from reproject import reproject_interp
from astropy.coordinates import SkyCoord
from analysis_tools.functions import delheader
from analysis_tools.functions import FitsCutter
from astropy.visualization import (SinhStretch, ImageNormalize)


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
            "hcg90": {"ra": 330.52343, "dec": -31.96680}
            }

def get_deep_optical(hcg_position):
    bands = ['g', 'r']
    
    for gal, pos in hcg_position.items():
        images = {}  # Store downloaded images
        headers = {}  # Store headers for reprojection
    
        for band in bands:
            download(ra=pos["ra"], dec=pos["dec"], mode='fits', bands=band, size=.6, layer='ls-dr10', pixscale=0.75)
    
            # Construct the filename based on the download
            filename_base = f"legacystamps_{pos['ra']}0_{pos['dec']}00_ls-dr10_{band}.fits"
            filename = f"legacystamps_{pos['ra']}0_{pos['dec']}00_ls-dr10.fits"
            os.system("mv " + filename + " " + filename_base)
            # Read the downloaded FITS file
            with fits.open(filename_base) as hdulist:
                images[band] = hdulist[0].data
                headers[band] = hdulist[0].header
    
        # Desired pixel scale in arcseconds per pixel
        target_header = headers['g'].copy()
        arcsec_per_pixel = 0.27 * 8
    
        # Convert arcsec_per_pixel to degrees per pixel since FITS WCS uses degrees
        deg_per_pixel = arcsec_per_pixel / 3600
    
        # Update the target header to reflect the new pixel scale
        target_header['CDELT1'] = -deg_per_pixel  # Negative for RA to indicate the direction
        target_header['CDELT2'] = deg_per_pixel   # Positive for Dec
    
        for band in bands:
            images[band], _ = reproject_interp((images[band], headers[band]), target_header)
    
        # Process the images: r, (g+r)/2, g
        r_image = images['r']
        g_image = (images['g'] + r_image) / 2 
        b_image = images['g']
    
        # Save individual images to use in make_rgbcube
        fits.writeto(f"{gal}_r.fits", r_image, target_header, overwrite=True)
        fits.writeto(f"{gal}_g.fits", g_image, target_header, overwrite=True)
        fits.writeto(f"{gal}_b.fits", b_image, target_header, overwrite=True)
        out_cube = generate_random_name() + "_cube.fits"
        aplpy.make_rgb_cube([f"{gal}_r.fits", f"{gal}_g.fits", f"{gal}_b.fits"], out_cube)
    
        pmax = 99.62
        out_png = generate_random_name() + ".png"
        #Create an RGB image from the cube
        aplpy.make_rgb_image(out_cube, out_png, stretch_r='log', stretch_g='log', stretch_b='log',
                            pmin_r=0.7, pmax_r=pmax, pmin_g=0.7, pmax_g=pmax, pmin_b=0.7, pmax_b=pmax)
    
        ## Clean up individual FITS files if no longer needed
        os.remove(f"{gal}_r.fits")
        os.remove(f"{gal}_g.fits")
        os.remove(f"{gal}_b.fits")
        img = Image.open(out_png)
    
        # Convert the image to an array for manipulation
        data = numpy.array(img)
        # Identify background pixels (assuming they are black [0, 0, 0])
        # You might need to adjust the threshold according to your image's background values
        background_threshold = 55
        mask = (data[:, :, 0] <= background_threshold) & \
            (data[:, :, 1] <= background_threshold) & \
            (data[:, :, 2] <= background_threshold)
    
        # Set identified background pixels to a grayish color, e.g., (70, 70, 70)
        pmaxstr = format(pmax, '.2f')
        j = 246
        data[mask] = [j, j, j]  # Adjust the RGB values to get the desired shade of gray
    
        # Convert the array back to an image
        new_img = Image.fromarray(data)
    
        # Save the modified image
        new_img.save(out_png[:-4] + "_grayish.png")
        return {"fits": out_cube[:-5] + "_2d.fits", "png": out_png[:-4] + "_grayish.png"}

def rect_coord(gal):
    """ gal e.g. hcg16"""
    coord_hcg = SkyCoord(zoom_rect[gal]["zoom_center"][0], 
                         zoom_rect[gal]["zoom_center"][1], frame='icrs')
    center_ra = coord_hcg.ra.deg
    center_dec = coord_hcg.dec.deg
    return {"center_ra": center_ra, "center_dec": center_dec}

fig_prop = open("config/figure_properties.yaml")
fig_prop = yaml.safe_load(fig_prop)

def show_fits(fits_map, output,  gal, bmin, bmaj, bpa, img, levels=[1], 
              rectsize=0.5, colmap = plt.cm.gray_r,contcol="black", 
              showcontour='no', showzoomarea='no', text_label='', 
              cbar_label='', optical='no', ellipsecolor='black', 
              rgbcube='', rectanglecolor='black', linecolor='black', 
              showcolorbar='yes', textcolor='black', vmin_vmax_plot='no', 
              vmin=1, vmax=1, addlabel='no', ra=[], dec=[], sofia_id='1',
              cat='', mom1_basename=''):
        
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
            png_image = generate_random_name()+'delete.png'
            aplpy.make_rgb_image(rgbcube, png_image,
                             pmin_g = 1,
                             pmax_g = 99,
                             pmin_r = 1,
                             pmax_r = 99,
                             pmin_b = 1,
                             pmax_b = 99)
            if gal == "hcg90" and 'zoom' in fits_map:
                # Example usage
                png_image = hcg90_opt["png"]
                rgb2d = hcg90_opt["fits"]
            else: 
                rgb2d = rgbcube[:-5]+'_2d.fits'

            fits_fig = aplpy.FITSFigure(rgb2d, figure=fig_size, 
                                    subplot=[fig_prop['general']['x0'],
                                    fig_prop['general']['y0'],
                                    fig_prop['general']['w'],
                                    fig_prop['general']['h']
                                    ])
            fits_fig.show_contour(fits_map, levels=levels, colors=contcol, smooth=1)
            fits_fig.show_rgb(png_image)
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
    args = parser.parse_args()
    return args

custom_font = figure_properties.FontManager(
        ['/opt/fits-tools/src/analysis_tools/'], 'tex gyre heros')

def cut_fits(f, rectsize, f_out):
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
    cutter_fits.cut_fits2D(rect_coord(gal)["center_ra"],
                           rect_coord(gal)["center_dec"], 
                           rectsize, f_out)
    return f_out

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
        file_suffix = f'_mom{mom_type[-1]}_pbc_zoompres.pdf'
        output = f"results/publication/figures/{gal}{file_suffix}"
        output_fig = [output, output.replace('_zoompres.pdf', '_largepres.pdf')]
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

def generate_random_name(length=5):
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
def plot_mom1_sources(fits_fig, cat, mom1_basename):
    """ Plots mom1 with each sources having their own vmin vmax to 
        highlight their rotational patterns

    Parameters
    ----------
    fits_fig: str
        aplpy.FITSFigure 
    cat: str
        A .txt file containing the SoFiA source catalogue
    mom1: str
        Moment 1 maps of the sources from SoFiA run (in the cubeletes directory) 
    """
    j = 1
    for a in range(len(cat)):
        f = mom1_basename + str(j) + '_mom1.fits'
        gal = os.path.basename(f)[:5]
        # RA, DEC of the sources.
        desired_ra, desired_dec = float(cat[a][38]), float(cat[a][39])
        # load the second image
        skycoord = SkyCoord(ra=desired_ra*units.degree, 
                            dec=desired_dec*units.degree, frame='fk5')
        x, y = skycoord.to_pixel(fits_fig._wcs)
        ax = fits_fig._figure.axes[0]

        new_array = fits.open(f)[0].data / 1000
        vmin, vmax = get_percentile_vmin_vmax(f, 1, 99)
        vmin = vmin / 1000 
        vmax = vmax / 1000
        min_value = numpy.nanmin(new_array)
        max_value = numpy.nanmax(new_array)
        factor = 0.5
        extent=[x-factor*new_array.shape[1], x+factor*new_array.shape[1],
                y-factor*new_array.shape[0], y+factor*new_array.shape[0]]
        # Overlay the new array at the specified position
        ax.imshow(new_array, vmin=vmin, vmax=vmax, extent=extent,
                  origin='lower', transform=ax.transData, cmap=plt.cm.coolwarm)
        ax.set_facecolor('black')
        if gal in large_source.keys():
            if j == large_source[gal][0]:
                step = large_source[gal][1]
            else:
                step = 25
        else:
            step = 25
        levels = numpy.arange(min_value, max_value, step)
        #ax.contour(new_array, levels=levels, colors="black", extent=extent)
        j = j + 1

if __name__ == '__main__': 
    args = get_args()
    custom_font.setup_fonts()
    f = args.fitsfile
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
    
    fplot_large = cut_fits(fplot, zoom_rect[gal]["zoom_size"][1], fplot[:-5]+'_large.fits')
    fplot_zoom = cut_fits(fplot, zoom_rect[gal]["zoom_size"][0], fplot[:-5]+'_zoom.fits')
    snmap_large = cut_fits(snmap, zoom_rect[gal]["zoom_size"][1], fplot[:-5]+'_snmap_large.fits')
    snmap_zoom = cut_fits(snmap, zoom_rect[gal]["zoom_size"][0], fplot[:-5]+'_snmap_zoom.fits')
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
    if gal == 'hcg91':
        z_filter = glob.glob(folder + gal+'legacystamps_*-ifilter.fits')[0] # no z filter for hcg 91
    else:    
        z_filter = glob.glob(folder + gal+ 'legacystamps_*-ifilter.fits')[0]
        g_filter = glob.glob(folder + gal+'legacystamps_*-gfilter.fits')[0]
    data_g = normalize_and_save_fits(g_filter, generate_random_name() + '.fits')
    data_i = normalize_and_save_fits(i_filter, generate_random_name() + '.fits')    
    data_z = normalize_and_save_fits(z_filter, generate_random_name() + '.fits')
    # Cut optical data withn a rectangular box
    data_g_large = cut_fits(data_g, zoom_rect[gal]["zoom_size"][1],
                            generate_random_name() + '_g_large.fits')
    data_g_zoom = cut_fits(data_g, zoom_rect[gal]["zoom_size"][0], 
                           generate_random_name() + '_g_zoom.fits')         
    data_i_large = cut_fits(data_i, zoom_rect[gal]["zoom_size"][1],
                            generate_random_name() + '_i_large.fits')
    data_i_zoom = cut_fits(data_i, zoom_rect[gal]["zoom_size"][0], 
                           generate_random_name() + '_i_zoom.fits')         
    data_z_large = cut_fits(data_z, zoom_rect[gal]["zoom_size"][1],
                            generate_random_name() + '_z_large.fits')
    data_z_zoom = cut_fits(data_z, zoom_rect[gal]["zoom_size"][0], 
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

    j = 0
    for filename in [fplot_zoom, fplot_large]:
        if filename == fplot_large and mom_type == 'mom0':
            # Add source id in moment zero plot. 
            show_fits(
                      filename,
                      vmin=numpy.nanpercentile(img, 10),	
                      vmax=numpy.nanpercentile(img, 90), 
                      gal=gal,
                      img=data_img[j],
                      showzoomarea='no', 
                      showcontour='no',
                      addlabel='yes', ra=ra_list, 
                      dec=dec_list, sofia_id=sofia_id_list, 
                      bmin=bmin, bmaj=bmaj, bpa=bpa,
                      output=output_fig[j][:-4]+'_idpres.pdf', 
                      text_label=text_label, colmap=fig_cmap, 
                      cbar_label=cbar_label
                     )
        if  mom_type == 'mom1':
            # Plot moment 1 map with each source having their own vmin, vmax
            blank_fits = generate_blankfits(filename)
            show_fits(
                      blank_fits,
                      gal=gal,
                      img=data_img[j],
                      vmin_vmax_plot='yes',
                      showcontour='no',
                      showcolorbar='no',
                      optical='no',
                      cat=f_sofia,
                      showzoomarea=showzoomarea[j],
                      mom1_basename=mom1_basename,
                      bmin=bmin, bmaj=bmaj, bpa=bpa,
                      rectsize=zoom_rect[gal]["zoom_size"][0],
                      output=output_fig[j][:-4]+'_sourcespres.pdf', 
                      text_label=text_label, rgbcube=rgb_cube_large, 
                     )
        # plot moment zero map and moment 1 map for both large and zoomed area
        show_fits(
                  filename,
                  vmin=numpy.nanpercentile(img, 10),    
                  vmax=numpy.nanpercentile(img, 90), 
                  gal=gal, levels=levels,
                  img=data_img[j],
                  showzoomarea=showzoomarea[j], 
                  showcontour=showcontour,  
                  rectsize=zoom_rect[gal]["zoom_size"][0],
                  bmin=bmin, bmaj=bmaj, bpa=bpa,
                  output=output_fig[j], text_label=text_label,
                  colmap = fig_cmap, cbar_label=cbar_label
                 )
        j = j + 1

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

        fits_large_nhi = reproject_fits(img_large_nhi, 
                                        hdr_large, base_header_large, 
                                        generate_random_name() + '_large_nhi.fits') 
        fits_zoom_nhi = reproject_fits(img_zoom_nhi, 
                                       hdr_zoom, base_header_zoom, 
                                       generate_random_name() + '_zoom_nhi.fits')
        fits_nhi = [fits_zoom_nhi, fits_large_nhi]
        rgb_cubes = [rgb_cube_zoom, rgb_cube_large]
        contour_levels = base_contour * 2 ** numpy.arange(20) 
        contour_levels = contour_levels[contour_levels <= max_contour]
        print(contour_levels, " CONTOUR LEVELS of {0}".format(gal.upper()));
        output_fig_nhi = ["results/publication/figures/" + gal + "_coldens_pbc_zoompres.pdf",
                          "results/publication/figures/" + gal + "_coldens_pbc_largepres.pdf"]
        # Plot column density maps for large and zoomed area. 
        for k in range(2):
            show_fits(
                      fits_nhi[k],
                      gal=gal, levels=contour_levels,
                      img=data_img_nhi[k],
                      showzoomarea=showzoomarea[k], 
                      showcontour='yes',
                      rectsize=zoom_rect[gal]["zoom_size"][0],
                      bmin=bmin, bmaj=bmaj, bpa=bpa,
                      output=output_fig_nhi[k], text_label="Column density",
                      colmap=fig_cmap, cbar_label=cbar_label, textcolor='white',
                      optical='yes', rgbcube=rgb_cubes[k], rectanglecolor='white', 
                      contcol="#2698ff", ellipsecolor='white', linecolor='white', 
                      showcolorbar='no')
