import os
import numpy
import aplpy
import random
import string
from PIL import Image
from astropy.io import fits
from legacystamps import download
from reproject import reproject_interp
from astropy import units
import requests
from PIL import Image
from astropy.coordinates import SkyCoord
from analysis_tools.functions import FitsCutter
from io import BytesIO
from urllib.error import HTTPError
from astropy.io import ascii, fits
from astropy import units as u
from astroquery.skyview import SkyView
import numpy as np

def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="fits", color=False):
    """Get URL for images in the table

    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """

    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg", "png", "fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra, dec, size=size, filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0, len(table) // 2, len(table) - 1]]
        for i, param in enumerate(["red", "green", "blue"]):
            url = url + "&{}={}".format(param, table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase + filename)
    return url

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
    coord_hcg = SkyCoord(rect_center["ra"] * units.deg, 
                         rect_center["dec"] * units.deg, frame='icrs')
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

def getcolorim(ra, dec, size=240, output_size=None, filters="grizy", format="jpg"):
    """Get color image at a sky position

    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png")
    Returns the image
    """

    #     if format not in ("jpg","png"):
    #         raise ValueError("format must be jpg or png")
    url = geturl(ra, dec, size=size, filters=filters, output_size=output_size, format=format, color=True)
    r = requests.get(url)
    im = Image.open(BytesIO(r.content))
    return im


def get_skyview(hi_pos, opt_view=6*u.arcmin, survey='DSS2 Blue', cache=True):
    """Retrieve the optical image from a certain pointing.

    :param hi_pos: position in the HI data
    :type hi_pos: astropy SkyCoord object
    :param opt_view: size of the optical image, defaults to 6*u.arcmin
    :type opt_view: astropy.units.arcmin, optional
    :param survey: survey containing optical data, defaults to 'DSS2 Blue'
    :type survey: str, optional
    :param cache: are the downloaded ancillary files cached
    :type cache: bool
    :return: optical image
    :rtype: astropy HDUList
    """
    # DSS2 Blue images have a 1 arc/pix pixel scale, but retrieving ~the pixel scale doesn't work.
    opt_pixels = (opt_view.to(u.arcsec).value * 2).astype(int)

    # Get a survey image from SkyView:
    if hi_pos.frame.name == 'galactic':
        path = SkyView.get_images(position=hi_pos.to_string('hmsdms'), coordinates='galactic', width=opt_view[0],
                                  height=opt_view[-1], survey=[survey], pixels=opt_pixels, cache=cache)
    elif (not hi_pos.equinox) or (hi_pos.frame.name == 'icrs'):
        path = SkyView.get_images(position=hi_pos.to_string('hmsdms'), coordinates='ICRS', width=opt_view[0],
                                  height=opt_view[-1], survey=[survey], pixels=opt_pixels, cache=cache)
    # Note that there seems to be a bug in SkyView that it sometimes won't retrieve non-J2000.0.  Keep an eye on this!
    else:
        path = SkyView.get_images(position=hi_pos.to_string('hmsdms'), coordinates=hi_pos.equinox.value,
                                  width=opt_view[0], height=opt_view[-1], survey=[survey], pixels=opt_pixels,
                                  cache=cache)
    if len(path) != 0:
        print("\tSurvey image retrieved from {}.".format(survey))
        result = path[0]
    else:
        print("\tWARNING: No {} image retrieved.  Bug, or server error?  Try again later?".format(survey))
        result = None

    return result


def get_panstarrs(hi_pos, opt_view=6*u.arcmin):
    """Get PanSTARRS false color image and r-band fits (for the WCS).

    :param hi_pos: position in the HI data
    :type hi_pos: astropy SkyCoord object
    :param opt_view: size of the optical image, defaults to 6*u.arcmin
    :type opt_view: astropy.units.arcmin, optional
    :return: color image and FITS header
    :rtype: Tuple[color_im, FITS_header] TODO check exact types!
    """
    #  Get PanSTARRS false color image and r-band fits (for the WCS).
    pstar_pixsc = 0.25
    if len(opt_view) > 1:
        print("\tWARNING: PanSTARRS only returns square images; taking largest dimension.")
        opt_view = np.max(opt_view)
    path = geturl(hi_pos.ra.deg, hi_pos.dec.deg, size=int(opt_view.to(u.arcsec).value / pstar_pixsc),
                  filters="r", format="fits")

    if len(path) != 0:
        fits_head = fits.getheader(path[0])
        color_im = getcolorim(hi_pos.ra.deg, hi_pos.dec.deg, size=int(opt_view.to(u.arcsec).value / pstar_pixsc),
                              filters="gri")
        print("\tOptical false color image retrieved from PanSTARRS.")
    else:
        print("\tWARNING: No PanSTARRS false color image retrieved.  Server error or no PanSTARRS coverage?")
        fits_head = None
        color_im = None

    return color_im, fits_head


def get_decals(hi_pos, opt_view=6*u.arcmin, dev_dr=False, pixscale = 0.262):
    """Get DECaLS false color image and g-band fits (for the WCS).

    :param hi_pos: position in the HI data
    :type hi_pos: astropy SkyCoord object
    :param opt_view: size of the optical image, defaults to 6*u.arcmin
    :type opt_view: astropy.units.arcmin, optional
    :return: color image and FITS header
    :rtype: Tuple[color_im, FITS_header] TODO check exact types!
    """
    # Get DECaLS false color image and fits (for the WCS). Example URL for this script provided by John Wu.
    # default(?) arcsec/pixel
    dimen = (opt_view.to(u.arcsec).value / pixscale).astype(int)
    if dev_dr:
        url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra={}&dec={}&layer=ls-dr10&' \
              'pixscale={}&width={}&height={}&bands=g'.format(hi_pos.ra.deg, hi_pos.dec.deg, pixscale, dimen[0], dimen[-1])
    else:
        url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra={}&dec={}&layer=ls-dr9&' \
              'pixscale={}&width={}&height={}&bands=g'.format(hi_pos.ra.deg, hi_pos.dec.deg, pixscale, dimen[0], dimen[-1])

    try:
        fits_head = fits.getheader(url)
        r = requests.get(url.replace("fits", "jpg").split('bands')[0])
        color_im = Image.open(BytesIO(r.content))
    except HTTPError:
        print("\tWARNING: HTTP Error, no DECaLS false color image retrieved. Server error or no DECaLS coverage?")
        fits_head = None
        color_im = None

    return color_im, fits_head


def get_wise(hi_pos, opt_view=6*u.arcmin, survey='WISE W1'):
    """Retrieve a WISE image from a certain pointing.

    :param hi_pos: position in the HI data
    :type hi_pos: astropy SkyCoord object
    :param opt_view: size of the optical image, defaults to 6*u.arcmin
    :type opt_view: astropy.units.arcmin, optional
    :param survey: survey containing WISE data, defaults to 'WISE W1'
    :type survey: str, optional
    :return: optical image
    :rtype: astropy HDUList
    """

    # IRSA query only accepts things in J2000/ICRS, so:
    # Note 2MASS is available through astroquery! get_irsa --> get_wise
    hi_pos = hi_pos.transform_to('icrs')
    params = {'POS': '{},{}'.format(hi_pos.ra.deg, hi_pos.dec.deg)}

    api_endpoint = "https://irsa.ipac.caltech.edu/ibe/search/wise/allwise/p3am_cdd"
    r = requests.get(url=api_endpoint, params=params)
    tab = ascii.read(r.content.decode(), header_start=44, data_start=48, format='fixed_width')

    params['coadd_id'] = tab['coadd_id'][0]
    params['coaddgrp'] = tab['coadd_id'][0][:2]
    params['coadd_ra'] = tab['coadd_id'][0][:4]
    params['band'] = int(survey[-1])
    params['size'] = str(opt_view[0].value) + ',' + str(opt_view[-1].value) + str(opt_view[0].unit)

    path = str.format("/{coaddgrp:s}/{coadd_ra:s}/{coadd_id:s}/{coadd_id:s}-w{band:1d}-int-3.fits?"
                      "center={POS:s}&size={size:s}", **params)

    try:
        result = fits.open(api_endpoint.replace('search', 'data') + path)
    except:
        result = None

    return result

def extract_medians_from_fits(fits_file):
    # Open the FITS file
    with fits.open(fits_file) as hdul:
        # Assuming the data cube is in the primary HDU
        data_cube = numpy.squeeze(hdul[0].data)

        # Check if the data cube is 3D
        if data_cube.ndim != 3:
            raise ValueError("The FITS file does not contain a 3D data cube.")

        # Calculate the median of each channel
        medians = numpy.nanmedian(data_cube, axis=(1, 2))

        return medians

def nhitomom0(input_nhi, bmin, bmaj):
    """
    converting column density values to moment zero in Jy/beam*m/s

    Parameters
    ---------
    input_nhi: array 
        column density values to be converted to moment zero
    bmin: float 
        input FWHM of the minor axis of the beam
    bmaj: float 
        input FWHM of the major axis of the beam

    Return
    ------
    momzero: array
        moment zero values
    """
    if isinstance(input_nhi, numpy.ndarray):
        momentzero = input_nhi * (bmin * bmaj) / 1.104e+21
        return momentzero
    else:
        raise TypeError(" The input values are not valid NumPy arrays.")

def invert_image_colors(out_png):
    img = Image.open(out_png)

    # Convert the image to an array for manipulation
    data = numpy.array(img)

    # Invert the colors
    # Subtract the data from 255 to invert it. This makes bright areas dark and dark areas bright
    inverted_data = 255 - data

    # Convert the array back to an image
    new_img = Image.fromarray(inverted_data)
    # Save the modified image
    new_img.save(out_png[:-4] + "_inverted.png")
    return out_png[:-4] + "_inverted.png"

def mom0tonhi(input_mom0, bmin, bmaj):
    """
    converting moment zero in Jy/beam*m/s to column density

    Parameters
    ---------
    input_mom0: array 
        moment zero values
    bmin: float 
        input FWHM of the minor axis of the beam
    bmaj: float 
        input FWHM of the major axis of the beam

    Return
    ------
    momzero: array
        moment zero values
    """
    if isinstance(input_mom0, numpy.ndarray):
        nhi = input_mom0 * 1.104e+21 / (bmin * bmaj)
        return nhi
    else:
        raise TypeError(" The input values are not valid NumPy arrays.")


def generate_random_name(length=5):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(length)) + 'delete'

def get_deep_optical(hcg_position, label="legacystamps", 
        directory=None, size=0.56,pixscale=0.67, make_deep=None):
    bands = ['g', 'i']#, 'z', 'i']

    for gal, pos in hcg_position.items():
        images = {}  # Store downloaded images
        headers = {}  # Store headers for reprojection

        for band in bands:
            download(ra=pos["ra"], dec=pos["dec"], mode='fits', bands=band, size=size, layer='ls-dr10', pixscale=pixscale)
            ra = format(round(pos["ra"], 6), '.6f')
            dec = format(round(pos["dec"], 6), '.6f')
            # Construct the filename based on the download
            filename_base = f"{label}_{ra}_{dec}_ls-dr10-{band}.fits"
            filename = f"legacystamps_{ra}_{dec}_ls-dr10.fits"
            if directory:
                os.system("mv " + filename + " " + directory+filename_base)
                with fits.open(directory + filename_base) as hdulist:
                    images[band] = hdulist[0].data
                    headers[band] = hdulist[0].header
            else:
                os.system("mv " + filename + " " + filename_base)
                with fits.open(filename_base) as hdulist:
                    images[band] = hdulist[0].data
                    headers[band] = hdulist[0].header
            # Read the downloaded FITS file
        if make_deep:
            # Desired pixel scale in arcseconds per pixel
            target_header = headers['g'].copy()
            # arcsec_per_pixel = 0.27 * 8

            # # Convert arcsec_per_pixel to degrees per pixel since FITS WCS uses degrees
            # deg_per_pixel = arcsec_per_pixel / 3600

            # # Update the target header to reflect the new pixel scale
            # target_header['CDELT1'] = -deg_per_pixel  # Negative for RA to indicate the direction
            # target_header['CDELT2'] = deg_per_pixel   # Positive for Dec

            for band in bands:
                images[band], _ = reproject_interp((images[band], headers[band]), target_header)

            # Process the images: r, (g+r)/2, g
            r_image = images[bands[0]]
            g_image = (images['g']) + (r_image/2)
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
            gray_img = img.convert('L')
            hist = gray_img.histogram()

            # Find the most frequent pixel value, which is likely to be the background
            background_value = hist.index(max(hist))

            print("Background pixel value:", background_value)
            # Convert the image to an array for manipulation
            data = numpy.array(img)
            # Identify background pixels (assuming they are black [0, 0, 0])
            # You might need to adjust the threshold according to your image's background values
            background_threshold = background_value*1.2
            mask = (data[:, :, 0] <= background_threshold) & \
                (data[:, :, 1] <= background_threshold) & \
                (data[:, :, 2] <= background_threshold)

            # Set identified background pixels to a grayish color, e.g., (70, 70, 70)
            pmaxstr = format(pmax, '.2f')
            j = 229
            t = 228
            b = 226
            data[mask] = [j, t, b]  # Adjust the RGB values to get the desired shade of gray

            # Convert the array back to an image
            new_img = Image.fromarray(data)

            # Save the modified image
            new_img.save(out_png[:-4] + "_grayish.png")
            return {"fits": out_cube[:-5] + "_2d.fits", "png": out_png[:-4] + "_grayish.png"}

# Function to read outputs from rules
def extract_outputs_from_rules(filename, endline="priority:"):
    outputs = []
    in_output_section = False  # Flag to indicate if we're inside the output section

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()  # Remove leading and trailing whitespace

            # Check if we're entering the output section
            if line.startswith('output:'):
                in_output_section = True
                continue  # Move to the next line

            # If we are in the output section and encounter 'priority:', stop reading
            if in_output_section and line.startswith(endline):
                break  # Exit the loop

            # If we are in the output section, collect output paths
            if in_output_section:
                # Remove any trailing commas and quotes
                outputs.extend([path.strip().strip('"') for path in line.split(',') if path.strip()])

    return outputs

def replace_galaxy_indices(outputs, indices):
    replaced_outputs = []
    for idx in indices:
        for output in outputs:
            # Replace {gal_idx} with the current index
            replaced_output = output.replace("{gal_idx}", idx)
            replaced_outputs.append(replaced_output)
    return replaced_outputs

