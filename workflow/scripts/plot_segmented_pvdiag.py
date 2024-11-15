import os
import csv
import yaml
import numpy
import aplpy
import argparse
import figure_properties
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units 
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from pvextractor import extract_pv_slice
from astropy.coordinates import SkyCoord
import pvextractor.geometry.path as pvpath
from analysis_tools.functions import radec2deg
from astropy.wcs.utils import skycoord_to_pixel
from analysis_tools.functions import FitsCutter
from math import sin, cos, sqrt, atan2, radians

# load custom font to change default matplotlib font
custom_font = figure_properties.FontManager(
        ['/opt/fits-tools/src/analysis_tools/'], 'tex gyre heros')

fig_prop = open("config/figure_properties.yaml")
fig_prop = yaml.safe_load(fig_prop)

def get_args():
    '''This function parses and returns arguments passed in'''
    description = 'Select dataset'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-m', '--mom_map', dest='mom_map',
                        help='Name of fits moment zero map')
    parser.add_argument('-c', '--cube', dest='cube',
                        help='Name of MeerKAT fits cube')
    parser.add_argument('-v', '--vlacube', dest='vlacube',
                        help='Name of VLA fits cube')

    args = parser.parse_args()
    return args


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


def plot_coords(ax, x_coords, y_coords, x_coords_hook, y_coords_hook):
    main_path_points = [35, 84, 170, 280, 366, -1]
    hook_path_points = [0, 1, 2]
    ax.plot(
        [x_coords_hook[i] for i in hook_path_points] 
        + [x_coords[i] for i in main_path_points],
        [y_coords_hook[i] for i in hook_path_points]
        + [y_coords[i] for i in main_path_points],
        color='black', lw=2
    )
    for i in main_path_points: 
        ax.plot(x_coords[i], y_coords[i], marker='o', color="black", markersize=8)   
    for i in hook_path_points: 
        ax.plot(x_coords_hook[i], y_coords_hook[i], 
                marker='o', color="black", markersize=8)   

def show_fits(path_coords, path_coords_hook, mom_map, bmin, bmaj, bpa,
              vmin, vmax, output, colmap = plt.cm.gray_r,  pathcol='black'):
    x_coords, y_coords = zip(*path_coords)
    x_coords_hook, y_coords_hook = zip(*path_coords_hook)  

    fig_size = plt.figure(figsize=(fig_prop['general']['figsize'][0],
                                   fig_prop['general']['figsize'][1]))

    fits_fig = aplpy.FITSFigure(mom_map, figure=fig_size,
                 subplot=[fig_prop['general']['x0'],
                 fig_prop['general']['y0'],
                 fig_prop['general']['w'],
                 fig_prop['general']['h']])
    #Add arrows labelling features
    fits_fig.show_arrows([32.5650,32.4848,32.3674,32.4498,32.5119,32.5367, 32.6499],
              [-10.3841,-10.3197,-10.1808,-10.0575,-10.1530,-10.2141, -10.2808],
              [0.0407,0.0248,0.0079,-0.0260,-0.0362,-0.0350, -0.03],
              [0.0212,0.0512,0.0389,-0.0433,-0.0189,0.0011, -0.03], color='r', width=1)
    
    fits_fig.add_label(32.5630, -10.3841, 'NGC848S tail', horizontalalignment='left', size=30)
    fits_fig.add_label(32.4828, -10.3197, 'SE tail', horizontalalignment='left', size=30)
    fits_fig.add_label(32.3654, -10.1808, 'NW tail', horizontalalignment='left', size=30)
    fits_fig.add_label(32.4538, -10.0575, 'NE tail', horizontalalignment='right', size=30)
    fits_fig.add_label(32.5139, -10.1530, 'E clump', horizontalalignment='right', size=30)
    fits_fig.add_label(32.5387, -10.2141, 'S clump', horizontalalignment='right', size=30)
    fits_fig.add_label(32.6499, -10.2808, 'Hook', horizontalalignment='right', size=30)
    fits_fig.show_colorscale(vmin=vmin, vmax=vmax, aspect='auto', cmap=colmap)  
    figprop = figure_properties.Figureprop(fits_fig)
    figprop.aplpy_figprop()
    fig_size.canvas.draw() 
    beampos = figprop.ra_dec_corner(mom_map, edge_offset_percent=5)
    x_word = beampos[0]
    y_word = beampos[1]
    # Show beam at the lower left corner of the image
    fits_fig.show_ellipses(x_word, y_word, bmin, bmaj,
                           angle=bpa, edgecolor="black")
    ax = plt.gca()
    figprop.axprop(ax, secaxtickhide=True)
    figprop.aplpy_colorbar(ylabel="$\mathrm{Flux~(Jy~beam^{-1}~km~s^{-1})}$", cbarloc='top')
    plot_coords(ax, x_coords, y_coords, x_coords_hook, y_coords_hook)
    plt.tight_layout()
    plt.savefig(output, bbox_inches="tight", dpi=fig_prop['general']['dpi']) 

# Read the CSV file
csv_file = 'data/hcg16/segmented_pv_coordinates.csv'
# ra dec position for the hook path
radec1 = radec2deg("2:10:04.5", "-10:10:03")
radec2 = radec2deg("2:10:20.7", "-10:17:05")
radec3 = radec2deg("2:10:36.1", "-10:20:30")

radecs = [
          radec2deg("2:10:04.5", "-10:10:03"),
          radec2deg("2:10:20.7", "-10:17:05"),
          radec2deg("2:10:36.1", "-10:20:30"),
        ]

def get_coordinates(hdu):
    U = list()
    for j in range(len(radecs)):
        Uj = skycoord_to_pixel(SkyCoord(ra=radecs[j]["ra"]*units.deg, 
                                        dec=radecs[j]["dec"]*units.deg,
                                        frame='icrs', equinox='J2000'), 
                                        WCS(hdu.header))
        U.append(Uj)
    
    coordinates_hook = [(int(U[0][0]), int(U[0][1])), 
                        (int(U[1][0]), int(U[1][1])), 
                        (int(U[2][0]), int(U[2][1]))]
    coordinates = []
    coordinates_sky = []
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            ra = float(row['RA'])
            dec = float(row['Dec'])
            U = skycoord_to_pixel(SkyCoord(ra=ra*units.deg, dec=dec*units.deg,
                                           frame='icrs', equinox='J2000'), 
                                           WCS(hdu.header))
            coordinates.append((int(U[0]), int(U[1])))
            coordinates_sky.append((ra, dec))
    return {
            "coordinates": coordinates, 
            "coordinates_hook": coordinates_hook, 
            "coordinates_sky": coordinates_sky
           }

coordinates_sky0 = [(radec1["ra"], radec1["dec"]), 
                    (radec2["ra"], radec2["dec"]), 
                    (radec3["ra"], radec3["dec"])]

def plot_pvdiagram(f, figname, pv, pvdrms_in=None, skycoord=None, 
                   elprop=None, tails=None, skycoord0=None, textlabel=''):
    hdr = pv.header
    cdelt1, cdelt2 = hdr["CDELT1"] * 60, hdr["CDELT2"]
    crval1, crval2 = hdr["CRVAL1"] * 60, hdr["CRVAL2"]
    crpix1, crpix2 = -hdr["CRPIX1"], -hdr["CRPIX2"]

    cdelt2 /= 1000 #convert to km/s
    crval2 /= 1000
    # Prepare velocity and position axes data
    vel = [(i + 1 - crpix2) * cdelt2 + crval2 for i in range(pv.data.shape[0])]
    xpos = [(i + 1 - crpix1) * cdelt1 + crval1 for i in range(pv.data.shape[1])]

    # Update PV header for plotting
    pv.header.update(CDELT1=cdelt1, CRVAL1=crval1, CDELT2=cdelt2, CRVAL2=crval2)
    del pv.header["WCSAXES"]

    # Write modified PV data to file
    pv.writeto(f[:-5] + "_arcsec.fits", overwrite=True)
    fig, ax1 = plt.subplots(figsize=(10, 5), 
                            subplot_kw={'projection': WCS(pv.header, fix=True, translate_units='shd')})

    # Create custom colormap
    colors_noise = plt.cm.gray(numpy.linspace(0, 1, 256))
    colors_galaxy = plt.cm.afmhot(numpy.linspace(1, 0.4, 256))
    all_colors = numpy.vstack((colors_noise, colors_galaxy))
    pvd_map = colors.LinearSegmentedColormap.from_list('pvd_map', all_colors)

    # Determine noise level for contour plotting
    if pvdrms_in is not None:
        pvd_rms = pvdrms_in
    else:
        pvd_rms = 1.4826 * numpy.nanmedian(numpy.abs(pv.data[pv.data < 0]))  # Assuming emission is all noise

    # Plot data with contours
    sigma_factor = 3
    divnorm = colors.TwoSlopeNorm(vmin=-sigma_factor * pvd_rms, 
                                  vcenter=+sigma_factor * pvd_rms, vmax=15 * pvd_rms)
    ax1.imshow(pv.data, cmap=pvd_map, aspect='auto', norm=divnorm)
    if numpy.nanmax(pv.data) > sigma_factor * pvd_rms:
        ax1.contour(pv.data, colors=['b'], levels=sigma_factor ** numpy.arange(1, 10) * pvd_rms)
    if numpy.nanmin(pv.data) < -sigma_factor * pvd_rms:
        ax1.contour(pv.data, colors=['w'], 
                    levels=-pvd_rms * sigma_factor ** numpy.arange(10, 0, -1), linestyles=['dashed'])

    # Customize plot appearance
    ax1.set_ylim(ax1.get_ylim()[0], ax1.get_ylim()[-1])
    ax1.set_xlim(ax1.get_xlim()[0], ax1.get_xlim()[-1])
    ax1.invert_yaxis()
    ax1.tick_params(direction='in', length=8.7, width=1.3, pad=12, 
                    labelsize=15)
    ax1.set_xlabel('Position [arcmin]', labelpad=1.5, 
                    fontsize=15)
    ax1.set_ylabel(r'$\mathrm{V~[km~s^{-1}]}$', labelpad=2, 
                   fontsize=15)
    ax1.text(0.01, 0.90, textlabel,
             transform=ax1.transAxes, color="yellow",
             fontsize=24)
    # Plot additional elements like ellipses and tails, ellipse 20 km/s and one hpbw
    if elprop:
        ellipse = Ellipse(elprop[0], width=elprop[1] / cdelt1, height= 20 / cdelt2, 
                          edgecolor='yellow', fc='none', lw=2)
        ax1.add_patch(ellipse)

    if tails:
        for tail in tails.values():
            ax1.text(tail["coord"][0] / cdelt1, (tail["coord"][1] - crval2) / cdelt2, tail["lab"], 
                     color="white", weight="normal", backgroundcolor="black", fontsize=12)

    #Initial distance set to 0
    distances = []
    
    # List of coordinate pairs for which distances need to be calculated
    coordinate_pairs = [
        (skycoord0[0], skycoord0[1]),
        (skycoord0[1], skycoord0[2]),
        (skycoord0[2], skycoord[35]),
        (skycoord[35], skycoord[84]),
        (skycoord[84], skycoord[170]),
        (skycoord[170], skycoord[280]),
        (skycoord[280], skycoord[366])
    ]
    dist = 0
    for i in range(0, len(coordinate_pairs)):
        ra1, dec1 = coordinate_pairs[i][0]
        ra2, dec2 = coordinate_pairs[i][1]
        # Calculate distance between the current and previous 
        # coordinate pair and adjust for scale
        dist +=  calculate_distance(ra1, dec1, ra2, dec2) / cdelt1
        # Accumulate distance
        distances.append(dist)    
    for dist in distances:
        ax1.vlines(x=dist, ymin=ax1.get_ylim()[0], ymax=ax1.get_ylim()[1]/cdelt2, 
                   linestyle='solid', color='black', lw=1)

    plt.savefig(figname, bbox_inches="tight", dpi=fig_prop['general']['dpi']) 

chan_start_cut = 140
chan_end_cut = 264
coord_hcg16 = SkyCoord('02h09m52.5s', '-10d13m20s', frame='icrs')
rectsize = 30
tails = [
        {
          'tail1' : { 
                 "lab":"Hook",
                 "coord": (8, 4150)
        },
          'tail2' : { 
                 "lab":"NGC 848S",
                 "coord": (13, 4150)
        },
          'tail3' : { 
                 "lab":"NGC 848",
                 "coord": (18, 3780)
        },
          'tail4' : { 
                 "lab":"SE tail",
                 "coord": (24, 3900)
        },
          'tail5' : { 
                 "lab":"HCG 16d",
                 "coord": (29, 3600)
        },
          'tail6' : { 
                 "lab":"HCG 16c",
                 "coord": (33, 4100)
        },
          'tail7' : { 
                 "lab":"NW tail",
                 "coord": (36, 3930)
        }
       },
        {
          'tail1' : { 
                 "lab":"Hook",
                 "coord": (8, 4150)
        },
          'tail2' : { 
                 "lab":"NGC 848S",
                 "coord": (13, 4150)
        },
          'tail3' : { 
                 "lab":"NGC 848",
                 "coord": (18, 3730)
        },
          'tail4' : { 
                 "lab":"SE tail",
                 "coord": (24, 3850)
        },
          'tail5' : { 
                 "lab":"HCG 16d",
                 "coord": (29, 3600)
        },
          'tail6' : { 
                 "lab":"HCG 16c",
                 "coord": (33, 4100)
        },
          'tail7' : { 
                 "lab":"NW tail",
                 "coord": (36, 3930)
        }
       },

]

if __name__ == '__main__':
    args = get_args()
    mom_map = args.mom_map
    hdu_mom = fits.open(mom_map)[0]
    img_mom = hdu_mom.data
    hdr_mom = hdu_mom.header
    # Cut cube in channels 
    cutter_cube_pbc = FitsCutter(args.cube)
    cube_pbc_cut = args.cube[:-5]+"_zoom_chancut.fits"
    cutter_cube_pbc.cut3D(coord_hcg16.ra.deg, coord_hcg16.dec.deg,
                          rectsize, chan_start_cut, chan_end_cut,
                          cube_pbc_cut)
    hdr_cube = fits.open(cube_pbc_cut)[0].header
    bmin = hdr_cube["BMIN"]
    bmaj = hdr_cube["BMAJ"]
    bpa = hdr_cube["BPA"]
    hpbw = numpy.sqrt(bmin * bmaj)
    # VLA data
    bmaj_vla, bmin_vla = 37.1952 / 3600, 30.3229 / 3600
    hpbw_vla = numpy.sqrt(bmin_vla * bmaj_vla)
    vla_cube = args.vlacube
    vla_cubename = os.path.basename(vla_cube)
    vla_outdir = os.path.dirname(args.cube) + "/vla_data/"
    if not os.path.exists(vla_outdir):
        os.makedirs(vla_outdir)
    else:
        pass
    pvdfits_vla = vla_outdir + "/" + vla_cubename[:-5] + "_pvdiag.fits"
    hdu_vla = fits.open(vla_cube)[0]
    hdus = [hdu_mom, hdu_vla]
    cubes = [cube_pbc_cut, vla_cube]
    beamsize = [hpbw, hpbw_vla]
    pv_outs = [cube_pbc_cut[:-5]+"_pvdiag.fits", pvdfits_vla]
    output = "outputs/publication/figures/hcg16_segmented_pv_view_mom0.pdf"
    figure_names = [output[:-14]+".pdf", output[:-14]+"_vla.pdf"]
    beampos = [(6, 118), (12, 38)]
    textlabel = ["MeerKAT", "VLA"]
    for j in range(2):
        coordinates = get_coordinates(hdus[j])["coordinates"]
        coordinates_hook = get_coordinates(hdus[j])["coordinates_hook"]
        coordinates_sky = get_coordinates(hdus[j])["coordinates_sky"]
        custom_font.setup_fonts()
        # Viewing pv diagram path in moment zero map
        if j == 0:
            show_fits(
                      mom_map = mom_map,
                      path_coords=coordinates, bmin=bmin, bpa=bpa, 
                      bmaj=bmaj, path_coords_hook=coordinates_hook, 
                      vmin=numpy.nanpercentile(img_mom, 1), 
                      vmax=numpy.nanpercentile(img_mom, 99), 
                      output=output,
                      colmap = plt.cm.RdYlBu_r
                     )
        # Plotting pv diagrams
        pvdiag_path = pvpath.Path(coordinates_hook + [coordinates[35], 
                                  coordinates[84], coordinates[170], 
                                  coordinates[280], coordinates[366], 
                                  coordinates[-1]])
        pv = extract_pv_slice(cube=cubes[j], 
                              path=pvdiag_path, spacing=1)
        pixels =  pv.header['NAXIS1']
        if j == 1: 
            if "TIMESYS" in pv.header:
                if  pv.header['TIMESYS'] == 'TAI':
                    pv.header['TIMESYS'] = 'tai'  # Correct the scale to lowercase
        pv.writeto(pv_outs[j], overwrite = True)
        plot_pvdiagram(
                       f=pv_outs[j], pv=pv, 
                       figname=figure_names[j],
                       skycoord=coordinates_sky, 
                       skycoord0=coordinates_sky0,
                       tails=tails[j], textlabel=textlabel[j], 
                       pvdrms_in=None, elprop=[beampos[j], beamsize[j] * 60]
                      )
        
#python workflow/scripts/plot_segmented_pvdiag.py -m outputs/sofia_pbc/hcg16/hcg16_line60_masked.pb_corr_vopt_mom0_nan_zoom.fits -c outputs/data/hcg16/hcg16_line60_masked.pb_corr_vopt.fits -v data/hcg16/vla_data/HCG16_CD_rob2_MS.fits
