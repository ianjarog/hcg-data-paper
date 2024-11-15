import os
import yaml
import glob
import aplpy
import numpy
import argparse
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
from analysis_tools.functions import FitsCutter
from utility_functions import generate_random_name
from math import sin, cos, sqrt, atan2, radians

hcgs_labels = {"hcg16":{"A":(32.35312,-10.13622), "B":(32.3361,-10.13309), "C":(32.41071,-10.14637), "D":(32.42872,-10.18394), "NGC848":(32.5735,-10.32145)},
              "hcg31":{"G":(75.43338,-4.28875), "Q":(75.40974,-4.22245), "A":(75.41146,-4.25946), "B":(75.39756,-4.26401), "C":(75.40751,-4.2577)},
              "hcg30":{"A":(69.07744,-2.83129), "B":(69.12619,-2.86656), "C":(69.097,-2.79985), "D":(69.15276,-2.84302)},
              "hcg90":{"A":(330.50889,-31.87014), "B":(330.53631,-31.99068), "C":(330.51423,-31.97451), "D":(330.52602,-31.99423)},
              "hcg91":{"A":(332.28174,-27.80984), "B":(332.3183,-27.73134), "C":(332.30881,-27.78241), "D":(332.28567,-27.80086),
                  "LEDA749936":(332.22916666666663, -27.783333333333335)},
              "hcg97":{"A":(356.84591,-2.30096), "B":(356.90753,-2.31716), "C":(356.84884,-2.35143), "D":(356.82867,-2.31329), "E":(356.83253,-2.28102)}}

offsets = {"hcg16":{"A":(0,1/60.),"B":(0,1/60.), "C":(0,1/60.), "D":(0,1/60.), "NGC848":(-5/60,0.)},
           "hcg31":{"Q":(0.5/60,0), "A":(0.5/60,0), "B":(-0.5/60,0), "C":(0,0.5/60), "G":(0.5/60.,0)},
           "hcg30":{"A":(1/60,0), "B":(-1/60,0), "C":(+0.6/60,0), "D":(-0.4/60,0)},
           "hcg91":{"A":(-1.2/60, 0), "B":(+0.5/60, 0), "C":(+0.7/60,0), "D":(0,0.5/60), "LEDA749936":(0,-0.5/60)},
           "hcg90":{"A":(0,3/60), "B":(0.5/60,0), "C":(-0.5/60,0), "D":(0,-3/60)},
           "hcg97":{"A":(1.1/60,0), "B":(-0.5/60,-1.2/60), "C":(+0.7/60,0), "D":(-0.8/60,0), "E":(0,0.6/60)}}

coord_hcg30c = {"ra":6.910151635171060e+01, "dec":-2.802908880849139e+00} #SoFiA source 5
coord_hcg30a = {"ra":6.907887052453162e+01, "dec":-2.845510123245648e+00} #SoFiA source 14 

def calculate_distance(ra1, dec1, ra2, dec2):
    """
    Calculates the angular distance between two points (ra1, dec1) (ra2, dec2) 
    given their right ascension and declination.

    Parameters:
    ra1 : float
        Right Ascension of the first point in degrees.
    dec1 : float
        Declination of the first point in degrees.
    ra2 : float
        Right Ascension of the second point in degrees.
    dec2 : float
        Declination of the second point in degrees.
    Returns:
    float
        Angular distance in arcminutes.
    """

    # Convert all coordinates from degrees to radians for calculation
    rad_ra1, rad_dec1, rad_ra2, rad_dec2 = map(radians, [ra1, dec1, ra2, dec2])

    # Calculate the change in right ascension
    delta_ra = rad_ra2 - rad_ra1

    # Calculate the 3D Cartesian coordinates for the vector difference
    x_component = cos(rad_dec2) * sin(delta_ra)
    y_component = cos(rad_dec1) * sin(rad_dec2) - sin(rad_dec1) * cos(rad_dec2) * cos(delta_ra)
    z_component = sin(rad_dec1) * sin(rad_dec2) + cos(rad_dec1) * cos(rad_dec2) * cos(delta_ra)

    # Calculate the angular distance using the atan2 function for a better numerical stability
    angular_distance_radians = atan2(sqrt(x_component ** 2 + y_component ** 2), z_component)

    # Convert the angular distance from radians to degrees and then to arcminutes
    angular_distance_degrees = angular_distance_radians * 180 / numpy.pi
    angular_distance_arcminutes = angular_distance_degrees * 60

    return angular_distance_arcminutes

D30=61
hcg30c_arcmin = calculate_distance(ra1=69.097, dec1=-2.79985, ra2=6.910151635171060e+01, dec2=-2.802908880849139e+00)
hcg30c_kpc = (4.848*D30*hcg30c_arcmin)/1000
hcg30a_arcmin = calculate_distance(ra1=69.07744, dec1=-2.83129, ra2=6.907887052453162e+01, dec2=-2.845510123245648e+00)
hcg30a_kpc = (4.848*D30*hcg30a_arcmin)/1000
print("HCG30C distance, arcmin", hcg30c_arcmin, hcg30c_kpc)
print("HCG30A distance, arcmin", hcg30a_arcmin, hcg30a_kpc)

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

def rect_coord(gal):
    """ gal e.g. hcg16"""
    coord_hcg = SkyCoord(zoom_rect[gal]["zoom_center"][0], 
                         zoom_rect[gal]["zoom_center"][1], frame='icrs')
    center_ra = coord_hcg.ra.deg
    center_dec = coord_hcg.dec.deg
    return {"center_ra": center_ra, "center_dec": center_dec}

fig_prop = open("config/figure_properties.yaml")
fig_prop = yaml.safe_load(fig_prop)

def show_fits(optical_map, fits_map, output,  gal, bmin, bmaj, bpa, img, levels=[1], 
              rectsize=0.5, colmap = plt.cm.gray_r,contcol="black", 
              showcontour='no', showzoomarea='no', text_label='', 
              cbar_label='', optical='no', ellipsecolor='black', 
              rectanglecolor='black', linecolor='black', 
              showcolorbar='yes', textcolor='black', vmin_vmax_plot='no', 
              vmin=1, vmax=1, addlabel='no', ra=[], dec=[], sofia_id='1',
              cat='', mom1_basename='', png_image='', recenter = False,
              recenter_coord_deg='', offset=5, cbarloc='left', x_word=None, y_word=None):
        
    size_deg = rectsize / 60
    fig_size = plt.figure(figsize=(fig_prop['general']['figsize'][0],
                                   fig_prop['general']['figsize'][1]))
    if optical != 'yes':
        fits_fig = aplpy.FITSFigure(optical_map, figure=fig_size, 
                     subplot=[fig_prop['general']['x0'],
                     fig_prop['general']['y0'],
                     fig_prop['general']['w'],
                     fig_prop['general']['h']])
        if recenter:
            fits_fig.recenter(recenter_coord_deg["ra"], recenter_coord_deg["dec"], width=size_deg, height=size_deg)

        if vmin_vmax_plot == 'no':
            fits_fig.show_colorscale(vmin=vmin, vmax=vmax, aspect='auto', 
                                  cmap=colmap, stretch='linear')
        else:
            plot_mom1_sources(fits_fig, cat, mom1_basename)

    else:
        if showzoomarea=="no": 
            folder = os.path.join("data", gal, "optical/") 
            fits_fig = aplpy.FITSFigure(glob.glob(folder + 'center_legacystamps_*-rfilter.fits')[0], figure=fig_size, 
                                    subplot=[fig_prop['general']['x0'],
                                    fig_prop['general']['y0'],
                                    fig_prop['general']['w'],
                                    fig_prop['general']['h']
                                    ])
        else:
            fits_fig = aplpy.FITSFigure(optical_map, figure=fig_size, 
                                subplot=[fig_prop['general']['x0'],
                                fig_prop['general']['y0'],
                                fig_prop['general']['w'],
                                fig_prop['general']['h']
                                ])
        fits_fig.show_grayscale(vmin=vmin, vmax=vmax, invert=True)
        fits_fig.show_contour(fits_map, levels=levels, colors=contcol, smooth=1)
        size_large = {"hcg16":40/60., "hcg91":84/60., "hcg30":95/60, "hcg97":90/60}
        size_deg_large = size_deg + size_large[gal]
        fits_fig.recenter(recenter_coord_deg["ra"], recenter_coord_deg["dec"], width=size_deg_large, height=size_deg_large)
        if recenter:
            members_label = list(hcgs_labels[gal].keys()) 
            fits_fig.recenter(recenter_coord_deg["ra"], recenter_coord_deg["dec"], width=size_deg, height=size_deg)
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

    if addlabel == 'yes':
        for j in range(len(ra)):
            fits_fig.add_label(ra[j], dec[j], sofia_id[j], color='red', 
                               size=fig_prop['general']['axtext_fontsize'])
    figprop = figure_properties.Figureprop(fits_fig)
    figprop.aplpy_figprop()
    fig_size.canvas.draw()
    beampos = figprop.ra_dec_corner(optical_map, edge_offset_percent=offset)
    if not x_word:
        x_word = beampos[0]
        y_word = beampos[1]
        # Show beam at the lower left corner of the image
        fits_fig.show_ellipses(x_word, y_word, bmin, bmaj,
                            angle=bpa, edgecolor=ellipsecolor)    
    else:
        fits_fig.show_ellipses(x_word, y_word, bmin, bmaj,
                    angle=bpa, edgecolor="black")  
    ax = plt.gca()
    # ellipse = Ellipse(xy=(x_word, y_word), width=bmaj, height=bmin,
    #                 angle=bpa, edgecolor='red', fc='None', lw=2)
    # ax.add_patch(ellipse)
    figprop.axprop(ax, secaxtickhide=True) 
    ax.text(0.6, 0.90, text_label, transform=ax.transAxes, color=textcolor, 
             fontsize=fig_prop['general']['axtext_fontsize'])
    ax.text(0.05, 0.90, gal.upper()[:3] + ' ' + gal[3:], 
             transform=ax.transAxes, color=textcolor, 
             fontsize=fig_prop['general']['axtext_fontsize'])
    if showcontour == 'yes' and optical != "yes":
        ax.contour(img, levels=levels, colors=contcol)
    if showcolorbar == 'yes':
        figprop.aplpy_colorbar(ylabel=cbar_label, cbarloc=cbarloc)
    if showzoomarea == 'yes':
        # Get the limits of the plot
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        
        # Convert the limits to world coordinates
        plot_top_right = fits_fig.pixel2world(xlim[1], ylim[1])
        plot_bottom_right = fits_fig.pixel2world(xlim[1], ylim[0])

        # Coordinates for the dashed lines
        if gal == "hcg90":
            top_right_ra = zoom_rect[gal]["zoom_center"][0] + size_deg / 2
            top_right_dec = zoom_rect[gal]["zoom_center"][1] + size_deg / 2
            bottom_left_ra = zoom_rect[gal]["zoom_center"][0] + size_deg / 2
            bottom_left_dec = zoom_rect[gal]["zoom_center"][1] - size_deg / 2

            fits_fig.show_rectangles(zoom_rect[gal]["zoom_center"][0], 
                zoom_rect[gal]["zoom_center"][1], size_deg, 
                size_deg, edgecolor=rectanglecolor, facecolor='none', alpha=0.8, 
                lw=2, linestyle='dashed')
        else:
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
    if gal == 'hcg90':
        cutter_fits.cut_fits2D(zoom_rect[gal]["zoom_center"][0],
                           zoom_rect[gal]["zoom_center"][1], 
                           rectsize, f_out)
    else:        
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
        file_suffix = f'_mom{mom_type[-1]}_pbc_zoom.pdf'
        output = f"outputs/publication/figures/{gal}{file_suffix}"
        output_fig = [output, output.replace('_zoom.pdf', '_large.pdf')]
        showcontour = {'mom1': 'yes', 'mom0': 'no'}  
        mom1_label = '$\mathrm{Velocity~(km~s^{-1})}$' 
        cmap_mom1 = plt.cm.coolwarm
        cmap_mom0 = plt.cm.RdYlBu_r
        mom0_label = '$\mathrm{Flux~(Jy~beam^{-1}~km~s^{-1})}$'
        mom0_textlabel = "moment 0"  
        mom1_textlabel = "moment 1"  
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
def plot_mom1_sources(fits_fig, cat, mom1_basename, ellipse=None):
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
        if gal in large_source.keys():
            if j == large_source[gal][0]:
                step = large_source[gal][1]
            else:
                step = 25
        else:
            step = 25
        levels = numpy.arange(min_value, max_value, step)
        ax.contour(new_array, levels=levels, colors="black", extent=extent)
        j = j + 1
plot_data = {
    "hcg16": {
        "pmin": -0.01427534817185605,
        "pmax": 0.04519556238746736,
        "pixscale": [1, 2.5],
        "ellipsepos": ("2:10:45", "-10:26:00")

    },
    "hcg30": {
        "pmin": -0.01427534817185605,
        "pmax": 0.04519556238746736,
        "pixscale": [0.262, 2.5],
        "ellipsepos": ("4:36:40", "-2:52:50")
    },
    "hcg31": {
        "pmin": -0.01427534817185605,
        "pmax": 0.04519556238746736,
        "pixscale": [1, 2.5],
        "ellipsepos": ("5:02:11", "-4:24:00")
    },
    "hcg90": {
        "pmin": -0.01427534817185605,
        "pmax": 0.04519556238746736,
        "pixscale": [0.262, 2.5],
        "ellipsepos": ("22:03:28", "-32:15:00")
    },
    "hcg91": {
        "pmin": -0.01427534817185605,
        "pmax": 0.04519556238746736,
        "pixscale": [2, 2.5],
        "ellipsepos": ("22:09:30", "-27:51:06")
    },
    "hcg97": {
        "pmin": -0.01427534817185605,
        "pmax": 0.04519556238746736,
        "pixscale": [0.262, 2.5],
        "ellipsepos": ("23:47:51", "-2:26:00"),
        "pixscale": [1, 2.5],
    }
}

def get_subcube(cube, chanstart, chanend):
    hdu_cube = fits.open(cube)[0]
    sub_cube = numpy.squeeze(hdu_cube.data)[chanstart:chanend+1, :, :]
    new_header_cube = hdu_cube.header.copy()
    new_header_cube = delheader(new_header_cube)
    if hdu_cube.header["CDELT3"] < 0 :
        hdrlist = numpy.arange(hdu_cube.header["CRVAL3"], -3.21E+09, hdu_cube.header["CDELT3"])
    else:
        hdrlist = numpy.arange(hdu_cube.header["CRVAL3"], 3.21E+09, hdu_cube.header["CDELT3"])
    new_header_cube['CRVAL3'] = hdrlist[chanstart]
    fits.writeto(cube[:-5]+"_subcube.fits", sub_cube, new_header_cube, overwrite=True) 
    return {"cube": cube[:-5]+"_subcube.fits", "cube_img": sub_cube, "cube_header": new_header_cube}


if __name__ == '__main__': 
    args = get_args()
    custom_font.setup_fonts()
    f = args.fitsfile
    gal = os.path.basename(f)[:5]
    # if gal=='hcg91' # plot only velocities between 6000 and 8200 - channels: 464 to 844
    #     cube = f.replace('sofia_pbc', 'data').replace('_mom0','')
    #     chanstart = 464
    #     chanend = 844
    #     sub_cube = get_subcube(cube, chanstart=chanstart, chanend=chanend)
    #     img_sub_cube = sub_cube["cube_img"]
    folder = os.path.join("data", gal, "optical/") 
    if gal == 'hcg91':
        r_filter = glob.glob(folder + 'legacystamps_*-ifilter.fits')[0]  
    else:
        r_filter = glob.glob(folder + 'legacystamps_*-rfilter.fits')[0]
    with fits.open(r_filter) as hdul:
        data = hdul[0].data  # Access the image data

    # Calculate vmin and vmax using the 5th and 95th percentiles
    vmin = numpy.percentile(data, 1)   # 5th percentile
    vmax = numpy.percentile(data, 99)  # 95th percentile    
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
    snmap_large_img = fits.open(snmap_large)[0].data

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
    img_large = larger_fits.data
    img_zoom = zoomed_fits.data
    imgs = [img_zoom, img_large]
    hdr_large = larger_fits.header
    data_img = [img_large] * 2
    if gal=="hcg90":
        recenter_coord_deg = {"ra": zoom_rect[gal]["zoom_center"][0], "dec": zoom_rect[gal]["zoom_center"][1]}
    else:
        recenter_coord_deg = radec2deg(ra=zoom_rect[gal]["zoom_center"][0], dec=zoom_rect[gal]["zoom_center"][1])
    # Plot column density maps for large and zoomed area. 
    recenter = [True, False]
    offset = [15, 5]  
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
    j = 0
    for filename in [fplot_large, fplot_large]:
        if j==1:
            x_word = None
            y_word = None
        else:
            x_word = radec2deg(plot_data[gal]["ellipsepos"][0], plot_data[gal]["ellipsepos"][1])["ra"]
            y_word = radec2deg(plot_data[gal]["ellipsepos"][0], plot_data[gal]["ellipsepos"][1])["dec"]
        if filename == fplot_large and mom_type == 'mom0':
            # Add source id in moment zero plot. 
            show_fits(
                      r_filter, filename,
                      vmin=numpy.nanpercentile(img, 10),	
                      vmax=numpy.nanpercentile(img, 90), 
                      gal=gal,
                      img=data_img[j],
                      showzoomarea='yes', 
                      showcontour='no',
                      recenter_coord_deg=recenter_coord_deg,
                      recenter=False,
                      addlabel='yes', ra=ra_list, 
                      dec=dec_list, sofia_id=sofia_id_list, 
                      bmin=bmin, bmaj=bmaj, bpa=bpa,
                      x_word=x_word, y_word=y_word,
                      output=output_fig[j][:-4]+'_id.pdf', 
                      text_label=text_label, colmap=fig_cmap, 
                      rectsize=zoom_rect[gal]["zoom_size"][0],
                      cbar_label=cbar_label
                     )
        # x_word, y_word = radec2deg("4:36:40", "-2:52:51.5")["ra"], radec2deg("4:36:40", "-2:52:51.5")["dec"]
        if j==0:
            data_img[0] = img_zoom
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
                      recenter=recenter[j],
                      showzoomarea=showzoomarea[j],
                      mom1_basename=mom1_basename,
                      bmin=bmin, bmaj=bmaj, bpa=bpa,
                      x_word=x_word, y_word=y_word,
                      recenter_coord_deg = recenter_coord_deg,
                      rectsize=zoom_rect[gal]["zoom_size"][0],
                      output=output_fig[j][:-4]+'_sources.pdf', 
                      text_label=text_label 
                     )
        # plot moment zero map and moment 1 map for both large and zoomed area
        show_fits(
                  r_filter, filename,
                  vmin=numpy.nanpercentile(imgs[j], 10),    
                  vmax=numpy.nanpercentile(imgs[j], 90), 
                  gal=gal, levels=levels,
                  img=data_img[j],
                  recenter=recenter[j],
                  offset=offset[j],
                  showzoomarea=showzoomarea[j], 
                  showcontour=showcontour,  
                  recenter_coord_deg = recenter_coord_deg,
                  rectsize=zoom_rect[gal]["zoom_size"][0],
                  bmin=bmin, bmaj=bmaj, bpa=bpa,
                  cbarloc='top',
                  x_word=x_word, y_word=y_word,
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
        #img_zoom_nhi[numpy.isnan(img_zoom_nhi)] = 0
        #img_zoom_nhi[snmap_zoom_img < 3] = 0
        #img_large_nhi[snmap_large_img < 3] = 0
        # define lowest contour according to signal-to-noise map: 
        mask = (snmap_large_img >= 2.5) & (snmap_large_img <= 3.5)
        filtered_img_large_nhi = numpy.where(mask, img_large_nhi, numpy.nan)
        base_contour = numpy.nanmedian(filtered_img_large_nhi)
        img_large_nhi[numpy.isnan(img_large_nhi)] = 0
        data_img_nhi = [img_large_nhi] * 2
        ra_deg = hdr_large["CRVAL1"]
        dec_deg = hdr_large["CRVAL2"]
        hi_pos = SkyCoord(ra_deg, dec_deg, unit="deg")
        contour_levels = base_contour * 2 ** numpy.arange(20) 
        contour_levels = contour_levels[contour_levels <= max_contour]
        print(contour_levels, " CONTOUR LEVELS of {0}".format(gal.upper()))
        output_fig_nhi = ["outputs/publication/figures/" + gal + "_coldens_pbc_zoom.pdf",
                          "outputs/publication/figures/" + gal + "_coldens_pbc_large.pdf"]
        pixscale = plot_data[gal]["pixscale"]
        opt_view = [numpy.array([zoom_rect[gal]["zoom_size"][0]+150,]) * units.arcmin, 
                    numpy.array([zoom_rect[gal]["zoom_size"][1],]) * units.arcmin]
        for k in range(2):
            image, header_decals = get_decals(hi_pos, opt_view=opt_view[k], dev_dr=True, pixscale=pixscale[k])
            fits_nhi = reproject_fits(img_large_nhi, 
                                      hdr_large, header_decals, 
                                      generate_random_name() + '_fits_nhi.fits') 
            if k==1:
                x_word = None
                y_word = None
            else:
                x_word = radec2deg(plot_data[gal]["ellipsepos"][0], plot_data[gal]["ellipsepos"][1])["ra"]
                y_word = radec2deg(plot_data[gal]["ellipsepos"][0], plot_data[gal]["ellipsepos"][1])["dec"]
            show_fits(
                      r_filter, fits_nhi,
                      gal=gal, levels=contour_levels,
                      img=data_img_nhi[k],
                      showzoomarea=showzoomarea[k],  
                      showcontour='yes', vmin=vmin, vmax=vmax,
                      rectsize=zoom_rect[gal]["zoom_size"][0],
                      bmin=bmin, bmaj=bmaj, bpa=bpa,
                      output=output_fig_nhi[k], text_label="column density",
                      colmap=fig_cmap, cbar_label=cbar_label, textcolor='black',
                      optical='yes', rectanglecolor='black', recenter=recenter[k], 
                      x_word=x_word, y_word=y_word,
                      recenter_coord_deg = recenter_coord_deg,
                      contcol="#ff0000", ellipsecolor='black', linecolor='black', 
                      showcolorbar='no', offset=offset[k])
#python workflow/scripts/plot_moments.py --fitsfile outputs/sofia_pbc/hcg16/hcg16_line60_masked.pb_corr_vopt_mom0.fits outputs/data/hcg16/hcg16_line60_masked.pb_corr_vopt.fits
#python workflow/scripts/plot_moments.py --fitsfile outputs/sofia_pbc_subcube/hcg91/hcg91_line60_masked.pb_corr_subcube_vopt_mom0.fits outputs/data/hcg90/hcg90_line60_masked.pb_corr_subcube_vopt.fits
