import os
import csv
import yaml
import numpy
import aplpy
import argparse
import subprocess
import figure_properties
from shutil import which
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as unit
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from utility_functions import mom0tonhi
from pvextractor import extract_pv_slice
from astropy.coordinates import SkyCoord
import pvextractor.geometry.path as pvpath
from analysis_tools.functions import radec2deg
from astropy.wcs.utils import skycoord_to_pixel
from analysis_tools.functions import FitsCutter
from utility_functions import get_deep_optical
from analysis_tools.functions import delheader
from math import sin, cos, sqrt, atan2, radians
from utility_functions import generate_random_name
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

hcg_position = {
            "hcg90": {"ra": 330.52343, "dec": -31.96680}
            }

def show_fits(path_coords, path_coords_sky, mom_map, bmin, bmaj, bpa,
              vmin, vmax, output, colmap = plt.cm.gray_r,  pathcol='black', contour_level=[10]):
    x_coords, y_coords = zip(*path_coords)

    fig_size = plt.figure(figsize=(fig_prop['general']['figsize'][0],
                                   fig_prop['general']['figsize'][1]))
    png_image = hcg90_opt["png"]
    rgb2d = hcg90_opt["fits"]

    fits_fig = aplpy.FITSFigure(rgb2d, figure=fig_size,
                 subplot=[fig_prop['general']['x0'],
                 fig_prop['general']['y0'],
                 fig_prop['general']['w'],
                 fig_prop['general']['h']])
    recenter_coord = ("22:02:05", "-32:01:05")
    recenter_deg = radec2deg(ra=recenter_coord[0], dec=recenter_coord[1])
    print("I AM RECENTERING")
    fits_fig.recenter(recenter_deg["ra"], recenter_deg["dec"], width = 0.28, height = 0.28)
    fits_fig.show_contour(mom_map, aspect='auto', levels=contour_level, colors='#ff1493') 
    print(path_coords_sky)
    ax = plt.gca()
    #fits_fig.show_lines([path_coords_sky], color="red")
    #fits_fig.show_markers(path_coords_sky[0,:], path_coords_sky[1,:], 
    #                      marker='o', s=10, c='red')
    fits_fig.show_rgb(png_image)
    #ax.contour(fits.open(mom_map)[0].data, levels=contour_level, color="red") 
    figprop = figure_properties.Figureprop(fits_fig)
    figprop.aplpy_figprop()
    fig_size.canvas.draw() 
    x_word = 330.83624999999995
    y_word = -32.231388888888894
    # Show beam at the lower left corner of the image
    fits_fig.show_ellipses(x_word, y_word, bmin, bmaj,
                           angle=bpa, edgecolor="red", linewidth=1.8)
    figprop.axprop(ax, secaxtickhide=True)
    plt.tight_layout()
    print("FIGURE SAVED TO ", output)
    plt.savefig(output, bbox_inches="tight", dpi=fig_prop['general']['dpi']) 

# Read the CSV file
csv_file = 'data/hcg90/segmented_pv_coordinates.csv'

def get_coordinates(hdu):
    coordinates_sky = []
    coordinates = []
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            ra = float(row['RA'])
            dec = float(row['Dec'])
            U = skycoord_to_pixel(SkyCoord(ra=ra*unit.deg, dec=dec*unit.deg,
                                           frame='icrs', equinox='J2000'), 
                                           WCS(hdu.header))
            coordinates.append((int(U[0]), int(U[1])))
            coordinates_sky.append((ra, dec))
    return {
            "coordinates": coordinates,
            "coordinates_sky": coordinates_sky,
           }

def plot_pvdiagram(f, figname, pv, pvdrms_in=None, skycoord=None, 
                   elprop=None, skycoord0=None, textlabel=''):
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
    fig, ax1 = plt.subplots(figsize=(10, 10), 
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
    ax1.minorticks_on()
    ax1.coords[0].display_minor_ticks(True)
    ax1.coords[1].display_minor_ticks(True)
    ax1.tick_params(axis='x',direction='in', length=8.7, top=True, right=True)
    ax1.tick_params(axis='y',direction='in', length=8.7, top=True, right=True)
    ax1.tick_params(which='minor', length=5)
    ax1.tick_params(which='minor', length=5)
    ax1.coords[0].set_ticks_position('all')  # RA or Longitude ticks
    ax1.coords[1].set_ticks_position('all')  # DEC or Latitude ticks
    ax1.tick_params(direction='in', width=1.3, pad=12, labelsize=22)
    ax1.set_xlabel('Position [arcmin]', labelpad=1.5, 
                    fontsize=22)
    ax1.set_ylabel(r'$\mathrm{Velocity~[km~s^{-1}]}$', labelpad=2, 
                   fontsize=22)
    ax1.text(0.01, 0.90, textlabel,
             transform=ax1.transAxes, color="yellow",
             fontsize=24)
    # Plot additional elements like ellipses and tails, ellipse 20 km/s and one hpbw
    if elprop:
        ellipse = Ellipse(elprop[0], width=elprop[1] / cdelt1, height= 20 / cdelt2, 
                          edgecolor='red', fc='none', lw=2)
        ax1.add_patch(ellipse)

    #Initial distance set to 0
    distances = []
    
    # List of coordinate pairs for which distances need to be calculated
    coordinate_pairs = [
        (skycoord[0], skycoord[1]),
        (skycoord[1], skycoord[2]),
        (skycoord[2], skycoord[3]),
        (skycoord[3], skycoord[4]),
        (skycoord[4], skycoord[5]),
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
        print(dist, "dist and XLIM", ax1.get_xlim())
        ax1.vlines(x=dist, ymin=ax1.get_ylim()[0], ymax=ax1.get_ylim()[1]/cdelt2, 
                   linestyle='solid', color='red', lw=2)

    plt.savefig(figname, bbox_inches="tight", dpi=fig_prop['general']['dpi']) 

coord_hcg90 = SkyCoord('02h09m52.5s', '-10d13m20s', frame='icrs')
rectsize = 30

def process_fits(base_filename, x1, y1, x2, y2, x3, y3):
    try:
        # Convert FITS file format (assuming xyin means input from xy plane)
        subprocess.run(["fits", f"in={base_filename}.fits", "op=xyin", f"out={base_filename}"], check=True)

        # Perform image subtraction or region selection
        region_spec = f"boxes({x1},{y1},{x2},{y2})({x3}, {y3})"
        subprocess.run(["imsub", f"in={base_filename}", f"out={base_filename}_imsub", f"region='{region_spec}'"], check=True)

        # Convert processed image back to FITS format
        subprocess.run(["fits", f"in={base_filename}_imsub", "op=xyout", f"out={base_filename}_imsub.fits"], check=True)

        # Remove intermediate data files
        subprocess.run(["rm", "-r", f"{base_filename}_imsub"], check=True)
        subprocess.run(["rm", "-r", f"{base_filename}"], check=True)

    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        return False
    return True

def eliminate_time(cat):
    '''Eliminates timestamp from sofia catalog. Updates the file

    Parameters
    ----------
    cat: str
        Path to sofia catalog
    '''
    # Read in the file
    with open(cat, 'r') as infile :
        lines = infile.readlines()
    # Replace the target string
    for i, line in enumerate(lines):
        if line[:7] == '# Time:':
            lines[i] = '# Time:\n'
    # Write the file out again
    with open(cat, 'w') as infile:
        for i in lines:
            infile.write(i)

def execute_sofia(executable, output_catalog, updated_parfile):
    """
    Running sofia.
 
    Parameters
    ----------
    executable: str
        Name of executable, either sofia or sofia2, 
        depends on how SoFia-2 environment was set up
    output_catalog: str
        Name of sofia output catalogue
    updated_parfile: 
        updated sofia parameter file
    """
    print('Executing Sofia-2')
    subprocess.call([executable, updated_parfile])
    print("Catalogue ", output_catalog[:-4]+".xml")
    #subprocess.run(["sofia_image_pipeline", "-c", output_catalog[:-4]+".xml"])
    eliminate_time(output_catalog)

def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None

def update_parfile(parfile, output_path, datacube, weight_cube):
    '''Updates file with paramenters

    Parameters
    ----------
    parfile: str
        File contanining sofia parameters
    output_path: str
        Path of output file
    datacube: str
        Path to datacube
    scfind_threshold: float
        Sofia parameter scfind.threshold
    reliability_threshold: float
        Sofia parameter reliability.threshold
    reliability_minpix: int
        Sofia parameter reliability.minPixels
    reliability_enable: str
        Sofia parameter reliability.enable
    Returns
    -------
    updated_parfile: str
        Path of file with updated parameters
    '''
    datacube_path = datacube
    datacube_name = os.path.basename(datacube)
    updated_parfile = os.path.join(output_path, datacube_name.split('_')[0] + '_sofia_updated_vopt_tail.par')
    print(f'Updated parfile: {updated_parfile}')
    path_norm = os.path.normpath(output_path)
    final_path = os.path.basename(path_norm)
    with open(parfile, 'r') as filein, open(updated_parfile, 'w') as fileout:
        lines = filein.read().replace('output_path', path_norm+"/hcg90-center-sofia")
        lines = lines.replace('datacube', datacube_path+".fits")
        lines = lines.replace('outname', f'{datacube_name}'+".fits")
        lines = lines.replace('weight_cube', weight_cube[:-5]+"_imsub.fits")
        fileout.write(lines)
    return updated_parfile

def squared_cube(cube, output_cube):
    with fits.open(cube) as hdul:
        data_cube = hdul[0].data
        hdr = hdul[0].header
    # Square the data
    squared_data_cube = data_cube * data_cube
    fits.writeto(output_cube, squared_data_cube, hdr, overwrite=True)
    return output_cube

def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Directory created: {result_path}")
    else:
        # The directory already exists
        print(f"Directory already exists: {result_path}")

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
    mom_map = args.mom_map
    hdu_mom = fits.open(mom_map)[0]
    img_mom = hdu_mom.data
    hdr_mom = hdu_mom.header
    # Cut cubes 
    cube_or = args.cube
    weight_or = args.cube[:-15] + "_vopt.fits"
    cube_copy = os.path.basename(cube_or)[:10] + "_central_tail.fits"
    weight_copy = os.path.basename(weight_or)[:-5] + "_central_tail.fits"
    # copy files temporary to a local directory because miriad does not like long characters
    subprocess.run(["cp", cube_or, "./" + cube_copy], capture_output=True, text=True)
    subprocess.run(["cp", weight_or, "./" + weight_copy], capture_output=True, text=True)
    weight_cube = squared_cube(weight_copy, weight_copy[:-5] + "_squared.fits")
    result_path = "results/sofia_pbc/hcg90/hcg90-center/"
    sofia_path = "results/sofia_pbc/hcg90/hcg90-center/hcg90-center-sofia/"
    create_dir(result_path)
    create_dir(sofia_path)
    cut_coord = [("22:02:39.36", "-32:06:24"), ("22:01:39.65", "-31:53:40")]
    radec_deg_lower = radec2deg(cut_coord[0][0], cut_coord[0][1])
    radec_deg_upper = radec2deg(cut_coord[1][0], cut_coord[1][1])
    radec_pix_lower = skycoord_to_pixel(SkyCoord(ra=radec_deg_lower["ra"]*unit.deg, 
                                        dec=radec_deg_lower["dec"]*unit.deg,
                                        frame='icrs', equinox='J2000'),
                                        WCS(hdu_mom))     
    radec_pix_upper = skycoord_to_pixel(SkyCoord(ra=radec_deg_upper["ra"]*unit.deg, 
                                        dec=radec_deg_upper["dec"]*unit.deg,
                                        frame='icrs', equinox='J2000'),
                                        WCS(hdu_mom))    
    xmin,ymin, = int(radec_pix_lower[0]), int(radec_pix_lower[1]) 
    xmax,ymax, = int(radec_pix_upper[0]), int(radec_pix_upper[1]) 
    z1,z2 = 325, 514 
    # run miriad to cut cubes
    process_fits(cube_copy[:-5], xmin, ymin, xmax, ymax, z1, z2)
    process_fits(weight_cube[:-5], xmin, ymin, xmax, ymax, z1, z2)
    subprocess.run(["mv", cube_copy[:-5]+"_imsub.fits", result_path], capture_output=True, text=True)
    subprocess.run(["mv", weight_cube[:-5]+"_imsub.fits", result_path], capture_output=True, text=True)
    datacube_name = cube_copy[:-5]+"_imsub"
    # run SoFiA on the cutted cubes
    parfile = "config/sofia_pbc/hcg90_sofia_cut_center.par"
    output_catalog = os.path.join(sofia_path, f'{datacube_name}_cat.txt')
    if not os.path.isfile(output_catalog):
        print(f'Parfile: {parfile}')
        updated_parfile = update_parfile(parfile, result_path, result_path + datacube_name,
	     result_path + weight_cube)
        if is_tool('sofia'):
            execute_sofia('sofia', output_catalog, f"{updated_parfile}")
        elif is_tool('sofia2'):
            execute_sofia('sofia2', output_catalog, f"{updated_parfile}")
        else:
            print('SoFia-2 not available. Please install it.')
            sys.exit(1)
    else:
        print(f"We have already found the catalogue {output_catalog}. \
                Sofia will not be executed" )
    hdr_cube = fits.open(result_path + cube_copy[:-5]+"_imsub.fits")[0].header
    bmin = hdr_cube["BMIN"]
    bmaj = hdr_cube["BMAJ"]
    mom0_central_tail = sofia_path + cube_copy[:-5]+"_imsub_mom0.fits"
    hdu_mom0_central_tail = fits.open(mom0_central_tail)[0]
    img_mom0_central_tail = hdu_mom0_central_tail.data
    snmap = signal_to_noise_map(sofia_path + cube_copy[:-5]+"_imsub_noise.fits",
                                mom0_central_tail, sofia_path + cube_copy[:-5]+"_imsub_chan.fits")
    snmap_central_tail = fits.open(snmap)[0].data
    mask = (snmap_central_tail >= 2.5) & (snmap_central_tail <= 3.5)
    nhi = mom0tonhi(img_mom0_central_tail, bmin*3600, bmaj*3600)
    nhi_fits = sofia_path + cube_copy[:-5]+"_imsub_nhi.fits"
    fits.writeto(nhi_fits, nhi, hdu_mom0_central_tail.header, overwrite=True)
    filtered_img_large_nhi = numpy.where(mask, nhi, numpy.nan)
    base_contour = numpy.nanmedian(filtered_img_large_nhi)  
    levels = base_contour * 2 ** numpy.arange(10)
    print(levels, "LEVELS")
    bpa = hdr_cube["BPA"]
    hpbw = numpy.sqrt(bmin * bmaj)
    hdus = [hdu_mom0_central_tail]
    cubes = [result_path + cube_copy[:-5]+"_imsub.fits"]
    beamsize = [hpbw]
    pv_outs = [result_path + cube_copy[:-5]+"_imsub_pvdiag.fits"]
    output = "results/publication/figures/hcg90_segmented_pv_view_mom0_tail.pdf"
    figure_names = [output[:-19]+"_tail.pdf"]
    beampos = [(4, 174)]
    textlabel = [""]
    #optical image of hcg90
    pathdir = "results/publication/figures/"
    hcg90_opt = get_deep_optical(hcg_position)
    base_header_zoom = fits.open(hcg90_opt["fits"])[0].header
    for j in range(1):
        coordinates_sky = get_coordinates(hdus[j])["coordinates_sky"]
        coordinates_sky_array = numpy.array(coordinates_sky).T # necessary for show_lines
        coordinates = get_coordinates(hdus[j])["coordinates"]
        custom_font.setup_fonts()
        # Viewing pv diagram path in moment zero map
        if j == 0:
            show_fits(
                      mom_map=nhi_fits,path_coords_sky=coordinates_sky_array,
                      path_coords=coordinates, bmin=bmin, bpa=bpa, 
                      bmaj=bmaj, 
                      vmin=numpy.nanpercentile(img_mom, 10), 
                      vmax=numpy.nanpercentile(img_mom, 90), 
                      output=output,
                      colmap = plt.cm.RdYlBu_r,
                      contour_level = levels
                     )
        # Plotting pv diagrams
        print("PATH ", [coordinates[0], coordinates[1],
                                              coordinates[2], coordinates[3], coordinates[4],
                                                                                coordinates[5]])
        pvdiag_path = pvpath.Path([coordinates[0], coordinates[1], 
                                  coordinates[2], coordinates[3], coordinates[4]])
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
                       skycoord0=None,
                       textlabel=textlabel[j], 
                       pvdrms_in=None, elprop=[beampos[j], beamsize[j] * 60]
                      )
#python workflow/scripts/plot_segmented_pvdiag_hcg90_tail.py -c results/data/hcg90/hcg90_line60_masked.pb_corr_vopt.fits -m results/sofia_pbc/hcg90/hcg90_line60_masked.pb_corr_vopt_mom0.fits
