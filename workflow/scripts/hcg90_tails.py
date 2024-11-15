import os
import yaml
import aplpy
import numpy
import argparse
import pyspeckit
from astropy import wcs
from astropy import units 
import matplotlib as mpl
import figure_properties
from astropy import units
from astropy.io import fits
from matplotlib import colors
import matplotlib.pyplot as plt
from pvextractor.geometry import Path
from scipy.interpolate import interp1d
from pvextractor import extract_pv_slice
from analysis_tools.functions import delheader
from analysis_tools.functions import radec2deg
from utility_functions import cut_fits
from astropy.convolution import convolve, Box1DKernel

zoom_area = open("config/parameters.yaml")
zoom_rect = yaml.safe_load(zoom_area) # position  and size of zoom area

def pvdiagrams(filename, outname, centra, centdec, 
               pa, majlen, minlen, mindist, dirc, pv_width, factor=1):
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
    fitsfile = fits.open(filename)
    hdu = fitsfile[0]
    print("HDU_DATA_SHAPE", hdu.data.shape)
    hdrc = hdu.header.copy()
            
    w = wcs.WCS(hdrc, fitsfile)

    central_pix_x, central_pix_y = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_world2pix(centra,centdec, 0)

    # Make a grid
    endpoints_y_b = numpy.array([majlen*factor/2., -majlen/2.,  mindist  ,  mindist ,         0.,        0.,   -mindist,  -mindist])
    endpoints_x_b = numpy.array([       0.,         0., -minlen/2., minlen/2., -minlen/2., minlen/2., -minlen/2.-0.0001, minlen/2.+0.0001])

    # Rotate
    endpoints_x = endpoints_x_b*numpy.cos(numpy.pi*pa/180.)-endpoints_y_b*numpy.sin(numpy.pi*pa/180.)+central_pix_x
    endpoints_y = endpoints_x_b*numpy.sin(numpy.pi*pa/180.)+endpoints_y_b*numpy.cos(numpy.pi*pa/180.)+central_pix_y

    print("xst = ", tuple(endpoints_x))
    print("yst = ", tuple(endpoints_y))

    i = -1
    list_points = []
    list_points_opt = []
    for names in ['_pvmaj.fits','_pvmin_left.fits', '_pvmin_cent.fits', '_pvmin_right.fits']:
        i = i + 1
        endpoints = [(endpoints_x[2*i],endpoints_y[2*i]),(endpoints_x[2*i+1],endpoints_y[2*i+1])]
        if names == '_pvmaj.fits':
            print("ENDPOINTS")
            print(endpoints)
        xy = Path(endpoints, width=pv_width*units.arcsec)
        pv = extract_pv_slice(hdu, xy)
        list_points.append(xy) 
        list_points_opt.append(xy) 
        header = pv.header
        #    print header
        pixels =  header['NAXIS1']
        pv.header['CRPIX1'] = 0
        pv.header['CRVAL1'] = 0
        pv.writeto(dirc + outname + names, overwrite = True)

    return {"xst" : tuple(endpoints_x), "yst": tuple(endpoints_y), "path": list_points[0], "header": pv.header}

fig_prop = open("config/figure_properties.yaml")
fig_prop = yaml.safe_load(fig_prop)

def plot_spectra(spectrum, smoothed_x, smoothed_y, yerror, xmin, xmax, output):
        fig = plt.figure()
        ax_spec = fig.add_subplot(111)
        ax_spec.set_xlim(xmin, xmax)
        ax_spec.set_ylabel("Integrated flux [Jy]", labelpad=20, fontsize=25)
        ax_spec.set_xlabel(r"$\mathrm{Velocity~[km~s^{-1}]}$", labelpad=20, fontsize=25)
        ax_spec.tick_params(pad=16, labelsize=25)
        ax_spec.plot(spectrum[:,0], spectrum[:,1], color="blue", zorder=1)
        ax_spec.plot([xmin - 10, xmax + 10], [0, 0], '--', color='gray', zorder=4)
        ax_spec.errorbar(spectrum[:,0], spectrum[:,1], yerr=yerror, elinewidth=0.75, capsize=1, mec="blue", mfc="blue", zorder=2)
        ax_spec.plot(smoothed_x, smoothed_y, color="black", lw=2, zorder=3)
        ax_spec.minorticks_on()
        ax_spec.tick_params(direction='in', length=8.7, top=True, right=True)
        ax_spec.tick_params(which='minor', length=5, direction='in', top=True, right=True)
        fig2 = plt.gcf()
        fig2.set_figheight(11)
        fig2.set_figwidth(11)
        plt.tight_layout()
        plt.savefig(output, bbox_inches="tight")

def plotpvd(f, hdr, img, figname, side, pos1=0.08, 
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
    mpl.rc('font', size=15)
    fig = aplpy.FITSFigure(f2, dimensions=[0,1])
    fig.set_theme('publication')
    fig.ticks.set_tick_direction('in')
    fig.axis_labels.set_xpad(2)
    fig.axis_labels.set_ypad(2)
    fig.tick_labels.set_font(size = 25)
    fig.axis_labels.set_font(size = 25)
    colors_noise = plt.cm.gray(numpy.linspace(0, 1, 256))
    colors_galaxy = plt.cm.afmhot(numpy.linspace(1, 0.4, 256))
    all_colors = numpy.vstack((colors_noise, colors_galaxy))
    pvd_map = colors.LinearSegmentedColormap.from_list('pvd_map', all_colors)
    pvd = img
    if pvdrms_in:
        pvd_rms = pvdrms_in
    else:
        pvd_rms = 1.4826 * numpy.nanmedian(numpy.abs(pvd[pvd < 0])) # assuming emission is all noise 
    sigma_factor = 3
    divnorm = colors.TwoSlopeNorm(vmin=-sigma_factor*pvd_rms, vcenter=+sigma_factor*pvd_rms, vmax=15*pvd_rms)
    plt.gca().invert_yaxis()
    print(pvd_rms, "PVDRMS")
    ax4 = plt.gca()
    ax4.imshow(pvd, cmap=pvd_map, aspect='auto', norm=divnorm)
    if numpy.nanmax(pvd) > sigma_factor*pvd_rms:
        ax4.contour(pvd, colors=['b', ], levels=sigma_factor**numpy.arange(1, 10)*pvd_rms)
    # Plot negative contours
    if numpy.nanmin(pvd) < -sigma_factor*pvd_rms:
        ax4.contour(pvd, colors=['w', ], levels=-pvd_rms * sigma_factor**numpy.arange(10, 0, -1), linestyles=['dashed', ])
    ax4.set_ylim(ax4.get_ylim()[0], ax4.get_ylim()[1])
    ax4.set_xlim(ax4.get_xlim()[0], ax4.get_xlim()[1])
    ax4.text(0.8, 0.9, side, color=textcol, fontsize=20,transform=ax4.transAxes)
    ax4.text(pos1, pos2, label, transform=ax4.transAxes, color=textcol)
    ax4.invert_yaxis()
    fig2 = plt.gcf()
    fig2.set_figheight(12)
    fig2.set_figwidth(12)
    ax = plt.gca()
    ax.invert_yaxis()
    ax.tick_params(direction='in', length=8.7, width=1.3, pad=12)
    ax.tick_params(which='minor', length=5)
    ax.set_xlabel('Offset (arcmin)', labelpad=1.5)

    ax.set_ylabel(r'$\mathrm{Velocity~(km~s^{-1})}$', labelpad=1.5)
    #ax.set_ylabel(r'V [Mhz]', labelpad=1.5)

    plt.savefig(figname, bbox_inches="tight", dpi=400)

def show_fits(fits_map, output,  gal, bmin, bmaj, bpa, img, 
              beampos, pathcol="blue", levels=[1], 
              colmap = plt.cm.gray_r,contcol="black", 
              showcontour=None, text_label='', 
              cbar_label='', ellipsecolor='black', 
              showcolorbar='yes', textcolor='black', 
              vmin=1, vmax=1, showcolorscale='no', pvd_path=None,
              base_contour_snmap=None, show_yticks=True, ylabel='Dec (J2000)'):
        
    fig_size = plt.figure(figsize=(fig_prop['general']['figsize'][0],
                                   fig_prop['general']['figsize'][1]))
    fits_fig = aplpy.FITSFigure(fits_map, figure=fig_size, 
                    subplot=[fig_prop['general']['x0'],
                    fig_prop['general']['y0'],
                    fig_prop['general']['w'],
                    fig_prop['general']['h']])
    ax3 = fits_fig._figure.axes[0]
    if showcolorscale == 'yes' and 'snmap' in fits_map :
        im = ax3.imshow(numpy.abs(img), norm=norm, cmap=colmap, origin='lower')
    if showcolorscale == 'yes' and 'mom1' in fits_map :
        im = ax3.imshow(img, vmin=vmin, vmax=vmax, 
                        cmap=colmap, origin='lower')
    if base_contour_snmap:
        img[numpy.isnan(img)] = 0
        ax3.contour(numpy.abs(img), linewidths=2, levels=[base_contour_snmap ], colors=['k', ])
        # fits_fig.show_colorscale(vmin=vmin, vmax=vmax, aspect='auto', 
        #                     cmap=colmap, stretch='linear', norm=norm)

    if showcontour:
        fits_fig.show_contour(fits_map, levels=levels, colors=contcol)
    ## SHOW PVD path
    if pvd_path:
        pvd_path.show_on_axis(ax3, spacing=1,
                    edgecolor=pathcol, linestyle='-',
                    linewidth=0.75)
    figprop = figure_properties.Figureprop(fits_fig)
    figprop.aplpy_figprop()
    fig_size.canvas.draw()
    fits_fig.show_ellipses(beampos["ra"], beampos["dec"], bmin, bmaj,
                           angle=bpa, edgecolor=ellipsecolor)    
    ax = plt.gca()
    ###
    figprop.axprop(ax, secaxtickhide=True, show_yticks=show_yticks, ylabel=ylabel) 
    ax.text(0.6, 0.90, text_label, transform=ax.transAxes, color=textcolor, 
             fontsize=fig_prop['general']['axtext_fontsize'])
    # ax.text(0.05, 0.90, gal.upper()[:3] + ' ' + gal[3:], 
    #          transform=ax.transAxes, color=textcolor, 
    #          fontsize=fig_prop['general']['axtext_fontsize'])
    fig = plt.gcf()
    if showcolorbar == 'yes':
        cb_ax = fig.add_axes([0.1, 0.9, 0.8, 0.025])
        cbar = fig.colorbar(im, cax=cb_ax) 
        cbar = fig.colorbar(im, cax=cb_ax, orientation='horizontal')
        cbar.ax.tick_params(direction='in', pad=10, labelsize=25)
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')
        cbar.set_label(cbar_label, fontsize=25, labelpad=16)

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

def calculate_coldens(img, bmin, bmaj):
    # convert moment zero in Jy/beam * km/s to column density
    bmin_arc = bmin * 3600
    bmaj_arc = bmaj * 3600
    nhi = 1.104e+21 * img * 1000 / bmin_arc / bmaj_arc
    return nhi    

def calc_pvdrms(f, mode=None):
    pvd = fits.open(f)[0].data
    if mode == 'negative':
        pvd_rms = 1.4826 * numpy.nanmedian(numpy.abs(pvd[pvd < 0]))
    elif mode == 'all':
        pvd_rms = numpy.nanmedian(abs(pvd - numpy.nanmedian(pvd))) / 0.6745  
    else:
        raise ValueError('Invalide input mode must be negative or all')
    return pvd_rms

if __name__ == '__main__': 
    args = get_args()
    custom_font.setup_fonts()
    f = args.fitsfile
    cube = args.fitscube
    mask_cube = cube.replace('data', 'sofia_pbc').replace('vopt','vopt_mask')
    gal = os.path.basename(f)[:5]
    openf = fits.open(f)[0]
    hdr = openf.header
    cut_pos = radec2deg('22:02:10.0550438695', '-32:00:20') 
    beampos = radec2deg('22:02:40', '-32:06:40') # position the beam within the plo
    zoom_size = 14.5 # cut cube within zoom_size in arcmin
    cube_cut = cut_fits(cube, cut_pos, rectsize = zoom_size, 
                        f_out3d = cube[:-5]+"_tail_cube.fits", chanstart = 314, chanend = 511)
    noise_cube_cut = cut_fits(f[:-9] + "noise.fits", cut_pos, rectsize = zoom_size, 
                        f_out3d = cube[:-5]+"_tail_cube_noise.fits", chanstart = 314, chanend = 511)
    mask_cube_cut = cut_fits(mask_cube, cut_pos, rectsize = zoom_size, 
                        f_out3d = cube[:-5]+"_tail_cube_mask.fits", chanstart = 314, chanend = 511)
    hdu_cube_cut = fits.open(cube_cut)[0]
    hdu_mask_cube_cut = fits.open(mask_cube_cut)[0]
    hdu_noise_cube_cut = fits.open(noise_cube_cut)[0]
    hdr_cube_cut = hdu_cube_cut.header
    hdr2d = delheader(hdr_cube_cut) # delete any 4rth dimensions in header
    hdr2d = delheader(hdr2d, "3") # delete any 3rd dimensions in header
    img_cube_cut = numpy.squeeze(hdu_cube_cut.data)
    img_mask_cube_cut = numpy.squeeze(hdu_mask_cube_cut.data) 
    img_noise_cube_cut = numpy.squeeze(hdu_noise_cube_cut.data)
    beam_fac = 2*numpy.pi/((numpy.sqrt(8*numpy.log(2)))**2) 
    bmin = hdr["BMIN"]
    bmaj = hdr["BMAJ"]
    pixsize = hdr["CDELT1"] * 3600 # in arcsec
    pix1 = abs(hdr["CDELT1"]) 
    pix2 = abs(hdr["CDELT1"])
    Npix_beam = (beam_fac*bmaj*bmin)/(pix1*pix2)
    mask2d = numpy.nansum(img_mask_cube_cut, axis=0)
    spectrum = numpy.nansum(img_cube_cut[:, mask2d != 0], axis=1)
    vels = pyspeckit.Cube(cube_cut).xarr.value / 1000
    spectrum_jy =  spectrum / Npix_beam
    specs = numpy.c_[vels, spectrum_jy]
    velocity_resolution = numpy.mean(numpy.diff(vels))  # Assuming x is evenly spaced
    new_velres = 20
    kernel_size = int(numpy.round(new_velres / velocity_resolution))
    if kernel_size % 2 == 0:
        kernel_size += 1
    box_kernel = Box1DKernel(abs(kernel_size))

    y_smoothed = convolve(spectrum_jy, box_kernel)
    new_x = numpy.arange(min(vels), max(vels), new_velres)
    interpolate_flux = interp1d(vels, y_smoothed, kind='linear', fill_value="extrapolate")
    new_y = interpolate_flux(new_x)
    Npix = numpy.sum(mask2d != 0)
    median_noise = numpy.nanmedian(img_noise_cube_cut, axis=0)
    print(median_noise, Npix, Npix_beam, "PLOTS")
    yerror = vels * 0 + numpy.nanmedian(median_noise) * numpy.sqrt(Npix / Npix_beam)
    plot_spectra(specs, new_x, new_y, yerror, min(vels), max(vels), "outputs/publication/figures/" + gal + "_spectrum_tails.pdf")
    img_cube_cut[img_mask_cube_cut==0] = numpy.nan
    contributing_channels = numpy.sum(~numpy.isnan(img_cube_cut), axis=0)
    # calculate moment zero. 
    img_mom0 = numpy.nansum(img_cube_cut, axis=0) * abs(hdr_cube_cut['CDELT3'])  
    # calculate moment 1. 
    nz, ny, nx = img_cube_cut.shape
    vmap = numpy.ones([nz, ny, nx], dtype=numpy.float64)
    for vel in range(nz):
        vmap[vel, :, :] = vels[vel]
    img_vel = numpy.where(numpy.isnan(img_cube_cut), numpy.nan, img_cube_cut*vmap)
    mom1 = f[:-9] + "tails_mom1.fits"
    img_mom1 = numpy.nansum(img_vel , axis=0) / numpy.nansum(img_cube_cut, axis=0)
    img_mom1[img_mom1==0]=numpy.nan
    fits.writeto(mom1, img_mom1, hdr2d, overwrite=True)
    snmap_img = img_mom0 / (median_noise * numpy.sqrt(contributing_channels) * abs(hdr_cube_cut['CDELT3']))
    snmap = f[:-9]+"tails_snmap.fits"
    fits.writeto(snmap, snmap_img, hdr2d, overwrite=True)
    bpa = hdr["BPA"]
    pvpos = radec2deg('22:02:10.0550438695', '-32:01:09.3611094458') # position velocity diagram center

    ## Calculate pv diagrams

    pv_angle = 320
    pv_width = 49.7403671758 # in arcsec
    mindistance = 10 # in arcsec
    majoraxis = 12.1503041262 * 60 # in arcsec
    majlen = majoraxis / pixsize
    minlen = majoraxis / 3 /  pixsize
    mindist = mindistance / pixsize 

    pvdiag = pvdiagrams(filename=cube_cut, outname="pvd", 
                    centra=pvpos["ra"], centdec=pvpos["dec"], 
                    pa=pv_angle, majlen=majlen, minlen=minlen, 
                    mindist=mindist, dirc="outputs/sofia_pbc/hcg90/", 
                    pv_width=pv_width, factor=1)
    
    fig_pvd = "outputs/publication/figures/" + gal + "_pvd_tails.pdf"

    f_pv = "outputs/sofia_pbc/hcg90/pvd_pvmaj.fits"
    f_pv2 = "outputs/sofia_pbc/hcg90/pvd_pvmaj2.fits"
    # pvdiag2 = create_pv_diagram(cube_cut, f_pv2, pvpos["ra"], pvpos["dec"], 
    #                             length=majoraxis, angle=pv_angle, width=pv_width*units.arcsec)
    plotpvd(f=f_pv, hdr=pvdiag["header"], 
        img=fits.open(f_pv)[0].data, figname=fig_pvd, 
        side="", pvdrms_in=calc_pvdrms(f_pv, "negative"))

    mask = (snmap_img >= 2.5) & (snmap_img <= 3.5)
    output_fig_snmap = "outputs/publication/figures/" + gal + "_snmap_tails.pdf"
    output_fig_mom1 = "outputs/publication/figures/" + gal + "_mom1_tails.pdf"
    cmap_snr = colors.ListedColormap(['w', 'royalblue', 'limegreen', 'yellow', 'orange', 'red'])
    boundaries = [0, 1, 2, 3, 4, 5, 6]
    norm = colors.BoundaryNorm(boundaries, cmap_snr.N, clip=True)
    filtered_snmap_img = numpy.where(mask, snmap_img, numpy.nan)
    base_contour_snmap = numpy.nanmedian(filtered_snmap_img)
    print("Plotting S/N map of {0}".format(gal))
    show_fits(
            snmap, gal=gal, img=snmap_img,
            beampos=beampos,
            bmin=bmin, bmaj=bmaj, bpa=bpa,
            output=output_fig_snmap, text_label="",
            colmap=cmap_snr, cbar_label="Pixel SNR", textcolor='white',
            ellipsecolor='black', 
            showcolorbar='yes', showcolorscale='yes', 
            base_contour_snmap=base_contour_snmap, show_yticks=False, ylabel=' '
            )
    print("Plotting Moment 1 map of {0}".format(gal)) #numpy.arange(2300,2600,20)
    show_fits(
            mom1, gal=gal, img=img_mom1,
            beampos=beampos,
            bmin=bmin, bmaj=bmaj, bpa=bpa,
            output=output_fig_mom1, text_label="",
            colmap=plt.cm.coolwarm, cbar_label="$\mathrm{Velocity~[km~s^{-1}]}$", 
            vmin=numpy.nanpercentile(img_mom1, 10), vmax=numpy.nanpercentile(img_mom1, 90),
            textcolor='white', ellipsecolor='black', 
            showcolorbar='yes', showcolorscale='yes', show_yticks=False,
            showcontour='yes', levels=numpy.arange(2300,2600,20)[:12], ylabel=' '
            )
    print("CONTOUR LEVELS OF MOM1", numpy.arange(2300,2600,20)[:12])
    print("Plotting column density map of {0}".format(gal))
    max_contour = 30E+19
    # calculate column density from moment zero value
    img_zoom_nhi = calculate_coldens(img_mom0 / 1000, bmin, bmaj)
    filtered_img_large_nhi = numpy.where(mask, img_zoom_nhi, numpy.nan)
    base_contour = numpy.nanmedian(filtered_img_large_nhi)
    img_zoom_nhi[numpy.isnan(img_zoom_nhi)] = 0
    fits_zoom_nhi = f[:-9] + "tails_coldens.fits"
    fits.writeto(fits_zoom_nhi, img_zoom_nhi, hdr2d, overwrite=True)
    print(base_contour, "BASE_CONTOUR")
    contour_levels = base_contour * 2 ** numpy.arange(20) 
    contour_levels = contour_levels[contour_levels <= max_contour]
    print(contour_levels, " CONTOUR LEVELS of {0}".format(gal.upper()))
    output_fig_nhi = "outputs/publication/figures/" + gal + "_coldens_tails.pdf"
    # Plot column density maps for large and zoomed area. 
    show_fits(
            fits_zoom_nhi,
            gal=gal, levels=contour_levels[:4],
            img=img_zoom_nhi,
            showcontour='yes', beampos=beampos,
            bmin=bmin, bmaj=bmaj, bpa=bpa,
            output=output_fig_nhi, text_label="",
            textcolor='black',
            contcol="red", ellipsecolor='black', 
            showcolorbar='no', pvd_path = pvdiag["path"]
            )
#python workflow/scripts/hcg90_tails.pdf --fitsfile outputs/sofia_pbc/hcg90/hcg90_line60_masked.pb_corr_vopt_mom0.fits -c outputs/data/hcg90/hcg90_line60_masked.pb_corr_vopt.fits
