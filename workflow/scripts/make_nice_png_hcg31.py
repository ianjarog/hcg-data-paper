import os
import yaml
import glob
import aplpy
import numpy
import string
import random
import argparse
import figure_properties
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as units
import matplotlib.pyplot as plt
from PIL import Image, ImageEnhance
from matplotlib.patches import Ellipse
from reproject import reproject_interp
from utility_functions import get_decals
from astropy.coordinates import SkyCoord
from analysis_tools.functions import delheader
from analysis_tools.functions import radec2deg
from utility_functions import invert_image_colors
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, ListedColormap
from matplotlib.cm import get_cmap
from astropy.visualization import (SinhStretch, ImageNormalize)

hcg31_tidal = {"A":[("5:01:38.67", "-4:15:33.6"), "white"],
               "A1":[("5:01:39.74", "-4:15:12.4"), "blue"],
               "B":[("5:01:35.25", "-4:15:51.6"), "white"],
               "C":[("5:01:37.76", "-4:15:28.4"), "white"],
               "E1":[("5:01:37.23", "-4:15:48.9"), "white"],
               "E2":[("5:01:37.48", "-4:15:57.5"), "white"],
               "F1":[("5:01:39.71", "-4:16:22.2"), "blue"],
               "F2":[("5:01:40.15", "-4:16:27.5"), "blue"],
               "F3":[("5:01:40.62", "-4:16:36.5"), "blue"],
               "G":[("5:01:44.01", "-4:17:19.5"), 'white'],
               "H1":[("5:01:38.20", "-4:16:07.2"), "blue"],
               "H2":[("5:01:38.72", "-4:16:13.5"), "blue"],
            #    "Q":("5:01:38.30", "-4:13:20.9"),
            #    "R1":("5:01:34.33", "-4:12:56.7")
}
#file:///Users/ianja/Downloads/0602288.pdf
fig_prop = open("config/figure_properties.yaml")
fig_prop = yaml.safe_load(fig_prop)
zoom_area = open("config/parameters.yaml")
zoom_rect = yaml.safe_load(zoom_area) # position  and size of zoom area
def png_background(out_png, threshold=55, rgbval=245):
        img = Image.open(out_png)

        # Convert the image to an array for manipulation
        data = numpy.array(img)
        # Identify background pixels (assuming they are black [0, 0, 0])
        # You might need to adjust the threshold according to your image's background values
        background_threshold = threshold
        mask = (data[:, :, 0] <= background_threshold) & \
            (data[:, :, 1] <= background_threshold) & \
            (data[:, :, 2] <= background_threshold)

        # Set identified background pixels to a grayish color, e.g., (70, 70, 70)
        j = rgbval
        data[mask] = [j, j, j]  # Adjust the RGB values to get the desired shade of gray

        # Convert the array back to an image
        new_img = Image.fromarray(data)

        # Save the modified image
        new_img.save(out_png[:-4] + "_grayish.png")
        return out_png[:-4] + "_grayish.png"
ellipsepos = {"hcg31": {"ellipsepos": ("5:02:11", "-4:24:00")}}

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

def enhance_image(png_path):
    img = Image.open(png_path)
    enhancer = ImageEnhance.Contrast(img)
    img = enhancer.enhance(3)  # Increase contrast; adjust the factor as needed
    enhancer = ImageEnhance.Brightness(img)
    img = enhancer.enhance(1.2)  # Increase brightness; adjust the factor as needed
    new_path = png_path[:-4] + "_enhanced.png"
    img.save(new_path)
    return new_path

def show_fits(fits_map, output, gal, bmin, bmaj, bpa, png_image,
              ellipsecolor='black', rgbcube='', textcolor='black', 
              pmin_g=1, pmax_g=99, pmin_r=1, pmax_r=99, 
              pmin_b=1, pmax_b=99, contour_levels=[1],
              ):
        
    fig_size = plt.figure(figsize=(fig_prop['general']['figsize'][0],
                                   fig_prop['general']['figsize'][1]))
    aplpy.make_rgb_image(rgbcube, png_image, 
                     stretch_r='linear', 
                     stretch_g='linear', 
                     stretch_b='linear', 
                     pmin_g = pmin_g,
                     pmax_g = pmax_g,
                     pmin_r = pmin_r,
                     pmax_r = pmax_r,
                     pmin_b = pmin_b,
                     pmax_b = pmax_b)
    # png_image = png_background(png_image, threshold=0, rgbval=0)
    if gal != 'hcg31':
        png_image = enhance_image(png_image)
    # if gal=='hcg31':
    #     png_image = png_background(png_image, threshold=40, rgbval=255)
    # else:
        # png_image = "data/hcg31/optical/HCG_31_HST.jpg"
    png_image = invert_image_colors(png_image)
    fits_fig = aplpy.FITSFigure(fits_map, figure=fig_size, 
                            subplot=[fig_prop['general']['x0'],
                            fig_prop['general']['y0'],
                            fig_prop['general']['w'],
                            fig_prop['general']['h']
                            ])
    rectsize = zoom_rect[gal]["zoom_size"][0]
    size_deg = rectsize / 60
    recenter_coord_deg = radec2deg(ra="5:01:38", dec="-4:16:00")
    fits_fig.recenter(recenter_coord_deg["ra"], recenter_coord_deg["dec"], width=0.07, height=0.07)
    print("Contour levels", contour_levels)
    contour_levels = list(contour_levels)[:6]+[3.2e+21]
    print("Contour levels", contour_levels)
    fits_fig.show_contour(fits_map, levels=contour_levels, colors="#00BFFF", linewidths=2)
    for label, pos_col in hcg31_tidal.items():
        ra_deg = radec2deg(pos_col[0][0], pos_col[0][1])["ra"]
        dec_deg = radec2deg(pos_col[0][0], pos_col[0][1])["dec"]
        fits_fig.add_label(ra_deg, dec_deg, label, relative=False, color=pos_col[1], size=25, horizontalalignment='center')

    hdu = fits.open(fits_map)
    image_data = hdu[0].data
    print(numpy.nanmean(image_data), "MEAN_IMAGE_DATA-BEFORE", fits_map)
    # Normalize and scale the image data for display purposes
    # Adjust the scaling as needed for your data
    # if gal not in ['hcg30', 'hcg90']:
    #     vmin, vmax = get_percentile_vmin_vmax(fits_map, 5, 95)
    # else:
    #     vmin, vmax = get_percentile_vmin_vmax(fits_map, 1, 99)
    # image_data = numpy.clip(image_data, vmin, vmax)
    image_data[image_data==0] = numpy.nan
    norm = Normalize(vmin=numpy.nanmin(image_data), vmax=numpy.nanmax(image_data))
    image_data_normalized = norm(image_data)
    # Create an RGBA image with the FITS data and a fixed alpha for transparency
    rgba_image = numpy.zeros((image_data.shape[0], image_data.shape[1], 4))
    # image_data_normalized = image_data
    # image_data_normalized = 1 - image_data_normalized
    # image_data_normalized = apply_colormap(image_data)
    # rgba_image = rgba_image[:, :, :, 0]
    rgba_image[..., 0] = image_data_normalized * 0  # Red channel (set to 0 for cyan)
    rgba_image[..., 1] = image_data_normalized * 0.588  # Green channel
    rgba_image[..., 2] = image_data_normalized * 0.968 # Blue channel
    rgba_image[..., 3] = image_data_normalized * 0.7 # Alpha channel for transparency
    rgba_image[numpy.isnan(image_data), 3] = 0 
    pil_image = Image.fromarray((rgba_image * 255).astype('uint8'), 'RGBA')
    enhancer = ImageEnhance.Contrast(pil_image)
    enhanced_image = enhancer.enhance(2.0)  # Increase contrast; adjust the factor as needed
    enhanced_rgba = numpy.array(enhanced_image) / 255.0
    # rgba_image = apply_colormap(rgba_image)
    # rgba_image = rgba_image[:, :, :, 0]
    print("Shape of rgba_image before imshow:", rgba_image.shape)
    # fits_fig.show_rgb(png_image)
    png_image = apply_intensity_colormap(png_image)
    img = Image.open(png_image)
    image = numpy.array(img)[::-1]
    fits_fig.show_grayscale()
    ax3 = fits_fig._figure.axes[0]
    m = numpy.mean(image)
    s = numpy.std(image)
    ax3.imshow(image,vmin=m-s,vmax=m+30*s,cmap='gray')
    # Step 3: Overlay the modified FITS image
    # Use 'imshow' from matplotlib to overlay the semi-transparent FITS image
    # ax3.imshow(enhanced_rgba, origin='lower')
    figprop = figure_properties.Figureprop(fits_fig)
    fig_size.canvas.draw()
    beampos = figprop.ra_dec_corner(fits_map, edge_offset_percent=5)
    x_word = radec2deg("5:01:31", "-4:17:40")["ra"]
    y_word = radec2deg("5:01:31", "-4:17:40")["dec"]
    # Show beam at the lower left corner of the image
    fits_fig.show_ellipses(x_word, y_word, bmin, bmaj,
                           angle=bpa, edgecolor="black")    
    ax = plt.gca()
    figprop.axprop(ax, secaxtickhide=True) 
    fits_fig.tick_labels.set_font(size=fig_prop['general']['ticklabelsize'])
    fits_fig.axis_labels.set_font(size=fig_prop['general']['xylabel_fontsize'])
    # fits_fig.axis_labels.hide_y()
    # fits_fig.tick_labels.hide_y()
    # ax.text(0.05, 0.90, gal.upper()[:3] + ' ' + gal[3:], 
    #          transform=ax.transAxes, color=textcolor, 
    #          fontsize=fig_prop['general']['axtext_fontsize'])
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

def reproject_fits(image_data, image_header, base_header, output):
    """Reproject data array from fits to a new header and save it to a new fits file"""
    reprojected_image, footprint = reproject_interp((image_data, WCS(image_header)), base_header)
    reprojected_image[numpy.isnan(reprojected_image)] = 0
    hdu = fits.PrimaryHDU(reprojected_image, header=base_header)
    # Save the reprojected image as a new FITS file
    hdu.writeto(output, overwrite=True)
    return output

def calculate_coldens(img, bmin, bmaj):
    # convert moment zero in Jy/beam * km/s to column density
    bmin_arc = bmin * 3600
    bmaj_arc = bmaj * 3600
    nhi = 1.104e+21 * img * 1000 / bmin_arc / bmaj_arc
    return nhi  

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

if __name__ == '__main__': 
    args = get_args()
    custom_font.setup_fonts()
    f = args.fitsfile
    gal = os.path.basename(f)[:5]
    openf = fits.open(f)[0]
    img = openf.data / 1000
    hdr = openf.header
    bmin = hdr["BMIN"]
    bmaj = hdr["BMAJ"]
    bpa = hdr["BPA"]
    img[img==0] = numpy.nan
    img = calculate_coldens(img, bmin, bmaj)
    snmap = signal_to_noise_map(f[:-9] + "noise.fits" , 
                                f[:-9] + "mom0.fits", f[:-9] + "chan.fits")
    snmap_img = fits.open(snmap)[0].data
    mask = (snmap_img >= 2.5) & (snmap_img <= 3.5)
    filtered_img = numpy.where(mask, img, numpy.nan)
    base_contour = numpy.nanmedian(filtered_img)
    img[numpy.isnan(img)] = 0
    f_nan = f[:-5]+'_nan.fits'
    hdu_nan = fits.PrimaryHDU(data=img, header=hdr)
    hdu_nan.writeto(f_nan, overwrite=True)
    fplot = f_nan
    contour_levels = base_contour * 2 ** numpy.arange(20) 
    # prepare optical data
    folder = os.path.join("data", gal, "optical/") 
    g_filter = glob.glob(folder + 'legacystamps_*-gfilter.fits')[0]
    i_filter = glob.glob(folder + 'legacystamps_*-rfilter.fits')[0]
    if gal == 'hcg91':
        z_filter = glob.glob(folder + 'legacystamps_*-ifilter.fits')[0] # no z filter for hcg 91
    else:    
        z_filter = glob.glob(folder + 'legacystamps_*-zfilter.fits')[0]
    data_g = normalize_and_save_fits(g_filter, generate_random_name() + '.fits')
    data_i = normalize_and_save_fits(i_filter, generate_random_name() + '.fits')    
    data_z = normalize_and_save_fits(z_filter, generate_random_name() + '.fits')

    rgb_cube_zoom = generate_random_name() + '_rgb_cube_zoom.fits'
    aplpy.make_rgb_cube([data_z, data_i, data_g], rgb_cube_zoom)
    
    pmin_r, pmax_r = 0.4, 99.5  # Emphasize the brighter features more in the infrared.
    pmin_g, pmax_g = 0.4, 99.5  # Standard settings, balancing the mid-range.
    pmin_b, pmax_b = 0.4, 99.5  # Highlight subtle features in colder areas.  
    # This could be deleted
    zoom_area = open("config/parameters.yaml")
    zoom_rect = yaml.safe_load(zoom_area) # position  and size of zoom area
    opt_view = numpy.array([zoom_rect[gal]["zoom_size"][0]+60,]) * units.arcmin
    ra_dec = radec2deg(zoom_rect[gal]["zoom_center"][0], zoom_rect[gal]["zoom_center"][1])
    hi_pos = SkyCoord(ra_dec["ra"], ra_dec["dec"], unit="deg")
    png_image = generate_random_name()+'delete.png'
    image, header_decals = get_decals(hi_pos, opt_view=opt_view, dev_dr=True, pixscale=0.01)
    image.save(png_image)
    ########
    base_header_zoom = fits.open(rgb_cube_zoom[:-5] + "_2d.fits")[0].header
    f_nan_reproj = reproject_fits(img, hdr, base_header_zoom,
                                  generate_random_name() + '_zoom_hi.fits')
    
    ####
    # Plot column density maps for large and zoomed area. 
    output_fig = "outputs/publication/figures/" + gal + "_optical_mom0.pdf"
    show_fits(
              f_nan_reproj, pmin_r=pmin_r, pmax_r=pmax_r,
              pmin_g=pmin_g, pmax_g=pmax_g,
              pmin_b=pmin_b, pmax_b=pmax_b,
              gal=gal, bmin=bmin, bmaj=bmaj, bpa=bpa, png_image=png_image,
              output=output_fig, textcolor='black',
              rgbcube=rgb_cube_zoom, ellipsecolor='black', 
              contour_levels=contour_levels
             )
#python workflow/scripts/make_nice_png_hcg31.py  -f outputs/sofia_pbc/hcg31/hcg31_line60_masked.pb_corr_vopt_mom0.fits
