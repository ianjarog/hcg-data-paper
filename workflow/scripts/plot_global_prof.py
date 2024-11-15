import os
import json
import yaml
import numpy
import string
import argparse
import pyspeckit
import figure_properties
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def get_args():
    '''This function parses and returns arguments passed in'''
    description = 'Select dataset to plot'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-c', '--cube', dest='cube', \
                        help='Input data cube')
    parser.add_argument('-m', '--mask', dest='mask', \
                        help='Input 3D mask')
    parser.add_argument('-s', '--vlaspec', dest='vlaspec', \
                        help='VLA spectrum')
    args = parser.parse_args()
    return args

# load figure properties
fig_prop = open("config/figure_properties.yaml")
fig_prop = yaml.safe_load(fig_prop)

def axprop(axlist):
    for axes in axlist:
        axes.minorticks_on()
        axes.tick_params(which='both', direction='in', 
                length=fig_prop['axes']['major_length'], 
                width=fig_prop['axes']['tick_width'], 
                pad=fig_prop['axes']['tick_pad'], 
                labelsize=fig_prop['general']['ticklabelsize'])
        axes.tick_params(which='minor', 
                length=fig_prop['axes']['minor_length'])


f = open("data/hcgs/hcgs_gal.json")
hcgs_gal = json.load(f)

f_prop = open("config/parameters.yaml")
hcgs_prop = yaml.safe_load(f_prop)


fig = plt.figure(figsize=(10.5,10.5))

y_lim1 = {"hcg16": -25, 
          "hcg30": -10, 
          "hcg31": -25,
          "hcg91": -19,
          "hcg90": -49, 
          "hcg97": -9}

y_lim2 = {"hcg16": 200, 
          "hcg30": 100, 
          "hcg31": 400,
          "hcg91": 100,
          "hcg97": 50,
          "hcg90": 150}

vmin_vmax = {"hcg16": [3400, 4400], 
             "hcg30": [3500, 5500], 
             "hcg31": [3000, 5000],
             "hcg91": [6000, 8400],
             "hcg90": [1500, 3500],
             "hcg97": [5600, 8000]}

ax_bottom_ylim = {
          "hcg16": [-50, 50],
          "hcg30": [-60, 60],
          "hcg31": [-200, 200],
          "hcg91": [-50, 50],
          "hcg90": [-200, 200],
          "hcg97": [-25, 25]
        }
# def plot_global_prof(cube, cube_mask, vla_spec):
#     """Ploting global profile 
    
#     Parameters
#     ----------
#     cube: str
#        Data cube 
#     cube_mask: str
#        3D mask to mask the data cube
#     vla_spec: VLA integrated spectrum
#     """
#     cube_name = os.path.basename(cube) 
#     vla_spec = numpy.loadtxt(vla_spec)
#     gal = cube_name[:5]
#     cube_mask_hdu = fits.open(cube_mask)[0]
#     cube_mask_img = numpy.squeeze(cube_mask_hdu.data)
#     data_cube_hdu = fits.open(cube)[0]
#     data_cube_img = numpy.squeeze(data_cube_hdu.data)
#     cube_mask_img_float = cube_mask_img.astype(float)
#     cube_mask_img_float[cube_mask_img_float==0] = numpy.nan
#     data_cube_img_masked = data_cube_img.copy()
#     data_cube_img_masked[numpy.isnan(cube_mask_img_float)] = numpy.nan
#     masked_prof = numpy.nansum(numpy.nansum(data_cube_img_masked, axis=2), axis=1)
#     hdr = data_cube_hdu.header
#     bmaj = hdr["BMAJ"] * 3600
#     bmin = hdr["BMIN"] * 3600
#     pix1 = abs(hdr["CDELT1"] * 3600) 
#     distance = hcgs_prop[gal]['distance']
#     pix2 = abs(hdr["CDELT2"] * 3600) 
#     mom0 = cube_mask.replace("_mask.fits", "_mom0.fits")
#     img_mom0 = fits.open(mom0[:-5]+"_nan.fits")[0].data
#     mass_mom0_meerkat = numpy.nansum(img_mom0)
#     beam_fac = 2*numpy.pi / ((numpy.sqrt(8*numpy.log(2)))**2) 
#     mom0_jy = mass_mom0_meerkat / (beam_fac*bmaj*bmin / (pix1*pix2))
#     mom0_mass = 2.36E+5*distance**2*mom0_jy
#     masked_prof_jy = masked_prof / (beam_fac*bmaj*bmin / (pix1*pix2))
#     data_cube = pyspeckit.Cube(cube)
#     velocities = data_cube.xarr.value / 1000
#     vmin = min(velocities) #vmin_vmax[gal][0]
#     vmax = max(velocities) #vmin_vmax[gal][1]
#     vla_min = numpy.min(vla_spec[:, 0])
#     vla_max = numpy.max(vla_spec[:, 0])
#     velocity_range = (velocities >= vmin) & (velocities <= vmax)
#     velocity_range_vla = (vla_spec[:,0] >= vla_min) & (vla_spec[:,0] <= vla_max)
#     flux_range = masked_prof_jy[velocity_range]
#     dv = abs(velocities[0]-velocities[1])
#     dv_vla = abs(vla_spec[:,0][0] - vla_spec[:,0][1])
#     tot_flux = numpy.nansum(flux_range) * dv
#     vla_flux_range = vla_spec[:,1][velocity_range_vla]
#     tot_flux_vla = numpy.nansum(vla_flux_range) * dv_vla
#     mass = 235600 * (distance ** 2) * tot_flux  
#     mass_vla = 235600 * (distance ** 2) * tot_flux_vla / 1000
#     error_jy = 0.3349e-3 * dv / (beam_fac*bmaj*bmin/(pix1*pix2)) 
#     mass_error = error_jy * 235600 * (distance ** 2)
#     print(tot_flux, tot_flux_vla, error_jy, dv, dv_vla, mass, mom0_mass, mass_vla, mass_error)

#     fig = plt.figure(figsize=(fig_prop['general']['figsize'][0],
#                               fig_prop['general']['figsize'][1]))
#     ax = fig.add_axes([fig_prop['general']['x0'],
#                        fig_prop['general']['y0'],
#                        fig_prop['general']['w'],
#                        fig_prop['general']['h']])

#     ax.plot(velocities, masked_prof_jy * 1000, '-', color='#468499', 
#             lw=fig_prop['line_plot']['line_width'], label="MeerKAT")
#     ax.plot(vla_spec[:,0], vla_spec[:,1], '-', color='#ae0e52', 
#             lw=fig_prop['line_plot']['line_width'], label="VLA")
#     ax.set_ylim(y_lim1[gal], y_lim2[gal])
#     ax.hlines(0, xmin=velocities.min()-250, 
#         xmax=velocities.max()+250, color="#803367", lw=2, ls='--')
#     ax.legend(frameon=True, loc='upper left', bbox_to_anchor=(0, 0.98), 
#         fontsize=fig_prop['general']['axtext_fontsize'])
#     #ax = plt.gca()
#     #fig = plt.gcf()
#     pos_main = ax.get_position()
#     ax_bottom = fig.add_axes([pos_main.x0, pos_main.y0 - pos_main.height * 0.25, 
#         pos_main.width, pos_main.height * 0.25])
#     ax.set_position([pos_main.x0, pos_main.y0, pos_main.width, pos_main.height * 0.75])
#     # Ensure the bottom panel shares the x-axis with the main plot
#     plt.setp(ax.get_xticklabels(), visible=False)
#     ax.set_xlim(vmin, vmax)
#     ax_bottom.set_xlim(ax.get_xlim())
#     for k in range(len(hcgs_gal[gal])):
#         delta_ym = hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["delta_ym"]
#         delta_x = hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["delta_x"]
#         ax.vlines(hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["vel"], 
#             ymin=y_lim1[gal], ymax=y_lim2[gal], color="black", ls='dotted')
#         ax_bottom.vlines(hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["vel"], 
#             ymin=-300, ymax=y_lim2[gal], color="black", ls='dotted')
#         alphabet = list(string.ascii_lowercase)
#         if list(hcgs_gal[gal].keys())[k][-1] in alphabet:
#             ax.text(hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["vel"]+delta_x, 
#                 y_lim2[gal]+delta_ym , list(hcgs_gal[gal].keys())[k][-1], color="black",
#                 fontsize=fig_prop['general']['axtext_fontsize'])
#         else:
#             gal_label = (list(hcgs_gal[gal].keys())[k][:-3] 
#                         + " " + list(hcgs_gal[gal].keys())[k][-3:])
#             ax.text(hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["vel"] + delta_x, 
#                 y_lim2[gal]+delta_ym, gal_label, color="black",
#                 fontsize=fig_prop['general']['axtext_fontsize'])
#     plt.text(0.77, 0.86, gal[:-2].upper() + ' ' + gal[-2:], 
#              transform=ax.transAxes, color="black", 
#              fontsize=fig_prop['general']['axtext_fontsize'])
#     ax_bottom.hlines(0, xmin=velocities.min()-250, 
#              xmax=velocities.max()+250, color="k", lw=2, ls='--')
#     secax = ax.secondary_xaxis('top')
#     secaxr = ax.secondary_yaxis('right')
#     bottom_secax = ax_bottom.secondary_xaxis('top')
#     bottom_secaxr = ax_bottom.secondary_yaxis('right')
#     plt.setp(bottom_secaxr.get_yticklabels(), visible=False)
#     plt.setp(bottom_secax.get_xticklabels(), visible=False)
#     velocities_interp = velocities[velocity_range_vla]
#     spline_interp = interp1d(vla_spec[:, 0], vla_spec[:, 1], kind='linear')
#     vla_interp = spline_interp(velocities_interp)
#     residuals = (vla_interp - masked_prof_jy[velocity_range_vla] * 1000)
#     ax_bottom.plot(velocities_interp, residuals, '-', 
#             color='k', lw=fig_prop['line_plot']['line_width'])
#     secax.tick_params(labeltop=False)
#     secaxr.tick_params(labelright=False)
#     ax_bottom.set_xlabel(r"$\mathrm{Velocity~(km~s^{-1})}$", 
#             labelpad=fig_prop['general']['mat_xlabelpad'], 
#             fontsize=fig_prop['general']['xylabel_fontsize'])
#     ax.set_ylabel(r"$\mathrm{Flux~(mJy)}$", 
#             labelpad=fig_prop['general']['mat_ylabelpad'],
#             fontsize=fig_prop['general']['xylabel_fontsize'])
#     ax_bottom.text(0.03, 0.73, 'Residuals', transform=ax_bottom.transAxes, 
#           color="black", fontsize=fig_prop['general']['axtext_fontsize'])
#     ax_bottom.set_ylim(ax_bottom_ylim[gal][0], ax_bottom_ylim[gal][1])
#     axlist = [ax, secax, secaxr, ax_bottom, bottom_secax, bottom_secaxr]
#     axprop(axlist)
#     fig2 = plt.gcf()
#     fig2.set_figheight(fig_prop['general']['figsize'][0])
#     fig2.set_figwidth(fig_prop['general']['figsize'][1])
#     figure_path = "outputs/publication/figures/"
#     plt.savefig(figure_path + gal + '_global_profile.pdf', 
#             bbox_inches="tight", dpi=fig_prop['general']['dpi'])

def zero2nan(file_in, file_out):
    hdu = fits.open(file_in)[0]
    hdr = hdu.header
    hdrc = hdr.copy()
    img = numpy.squeeze(hdu.data)
    img[img==0] = numpy.nan
    fits.writeto(file_out, img, hdrc, overwrite=True)
    return file_out

def plot_global_prof(cube, cube_mask, vla_spec):
    """Plotting global profile 

    Parameters
    ----------
    cube: str
       Data cube 
    cube_mask: str
       3D mask to mask the data cube
    vla_spec: str
       VLA integrated spectrum file path
    """
    cube_name = os.path.basename(cube) 
    vla_spec = numpy.loadtxt(vla_spec)
    gal = cube_name[:5]
    cube_mask_hdu = fits.open(cube_mask)[0]
    cube_mask_img = numpy.squeeze(cube_mask_hdu.data)
    data_cube_hdu = fits.open(cube)[0]
    data_cube_img = numpy.squeeze(data_cube_hdu.data)
    cube_mask_img_float = cube_mask_img.astype(float)
    cube_mask_img_float[cube_mask_img_float == 0] = numpy.nan
    data_cube_img_masked = data_cube_img.copy()
    data_cube_img_masked[numpy.isnan(cube_mask_img_float)] = numpy.nan
    masked_prof = numpy.nansum(numpy.nansum(data_cube_img_masked, axis=2), axis=1)
    hdr = data_cube_hdu.header
    bmaj = hdr["BMAJ"] * 3600
    bmin = hdr["BMIN"] * 3600
    pix1 = abs(hdr["CDELT1"] * 3600) 
    distance = hcgs_prop[gal]['distance']
    pix2 = abs(hdr["CDELT2"] * 3600) 
    mom0 = cube_mask.replace("_mask.fits", "_mom0.fits")
    mom0_nan = zero2nan(mom0, mom0[:-5]+"_nan.fits") 
    img_mom0 = fits.open(mom0_nan)[0].data
    mass_mom0_meerkat = numpy.nansum(img_mom0)
    beam_fac = 2 * numpy.pi / ((numpy.sqrt(8 * numpy.log(2))) ** 2) 
    mom0_jy = mass_mom0_meerkat / (beam_fac * bmaj * bmin / (pix1 * pix2))
    mom0_mass = 2.36E+5 * distance ** 2 * mom0_jy
    masked_prof_jy = masked_prof / (beam_fac * bmaj * bmin / (pix1 * pix2))
    data_cube = pyspeckit.Cube(cube)
    velocities = data_cube.xarr.value / 1000
    vmin = min(velocities)
    vmax = max(velocities)
    vla_min = numpy.min(vla_spec[:, 0])
    vla_max = numpy.max(vla_spec[:, 0])
    velocity_range = (velocities >= vmin) & (velocities <= vmax)
    velocity_range_vla = (vla_spec[:,0] >= vmin) & (vla_spec[:,0] <= vmax)
    flux_range = masked_prof_jy[velocity_range]
    dv = abs(velocities[1] - velocities[0])
    dv_vla = abs(vla_spec[1, 0] - vla_spec[0, 0])
    tot_flux = numpy.nansum(flux_range) * dv
    vla_flux_range = vla_spec[velocity_range_vla, 1]
    tot_flux_vla = numpy.nansum(vla_flux_range) * dv_vla
    mass = 235600 * (distance ** 2) * tot_flux  
    mass_vla = 235600 * (distance ** 2) * tot_flux_vla / 1000
    error_jy = 0.3349e-3 * dv / (beam_fac * bmaj * bmin / (pix1 * pix2)) 
    mass_error = error_jy * 235600 * (distance ** 2)
    print(tot_flux, tot_flux_vla, error_jy, dv, dv_vla, mass, mom0_mass, mass_vla, mass_error)

    fig = plt.figure(figsize=(fig_prop['general']['figsize'][0],
                              fig_prop['general']['figsize'][1]))
    ax = fig.add_axes([fig_prop['general']['x0'],
                       fig_prop['general']['y0'],
                       fig_prop['general']['w'],
                       fig_prop['general']['h']])

    ax.plot(velocities, masked_prof_jy * 1000, '-', color='#468499', 
            lw=fig_prop['line_plot']['line_width'], label="MeerKAT")
    ax.plot(vla_spec[:,0], vla_spec[:,1], '-', color='#ae0e52', 
            lw=fig_prop['line_plot']['line_width'], label="VLA")
    if gal in ["hcg31"]:
        ax.set_xlim(3500, 4500)
    if gal in ["hcg91"]:
        ax.set_xlim(4000, 10000)
    if gal in ["hcg90"]:
        ax.set_xlim(1000, 5000)
    ax.set_ylim(y_lim1[gal], y_lim2[gal])
    ax.hlines(0, xmin=velocities.min() - 250, 
              xmax=velocities.max() + 250, color="#803367", lw=2, ls='--')
    ax.legend(frameon=True, loc='upper left', bbox_to_anchor=(0, 0.98), 
              fontsize=fig_prop['general']['axtext_fontsize'])

    pos_main = ax.get_position()
    ax_bottom = fig.add_axes([pos_main.x0, pos_main.y0 - pos_main.height * 0.25, 
                              pos_main.width, pos_main.height * 0.25])
    ax.set_position([pos_main.x0, pos_main.y0, pos_main.width, pos_main.height * 0.75])
    plt.setp(ax.get_xticklabels(), visible=False)
    # ax.set_xlim(vmin, vmax)
    ax_bottom.set_xlim(ax.get_xlim())

    for k in range(len(hcgs_gal[gal])):
        delta_ym = hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["delta_ym"]
        delta_x = hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["delta_x"]
        ax.vlines(hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["vel"], 
                  ymin=y_lim1[gal], ymax=y_lim2[gal], color="black", ls='dotted')
        ax_bottom.vlines(hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["vel"], 
                         ymin=-300, ymax=y_lim2[gal], color="black", ls='dotted')
        alphabet = list(string.ascii_lowercase)
        if list(hcgs_gal[gal].keys())[k][-1] in alphabet:
            ax.text(hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["vel"] + delta_x, 
                    y_lim2[gal] + delta_ym, list(hcgs_gal[gal].keys())[k][-1], color="black",
                    fontsize=fig_prop['general']['axtext_fontsize'])
        else:
            gal_label = (list(hcgs_gal[gal].keys())[k][:-3] 
                         + " " + list(hcgs_gal[gal].keys())[k][-3:])
            ax.text(hcgs_gal[gal][list(hcgs_gal[gal].keys())[k]]["vel"] + delta_x, 
                    y_lim2[gal] + delta_ym, gal_label, color="black",
                    fontsize=fig_prop['general']['axtext_fontsize'])
    plt.text(0.77, 0.86, gal[:-2].upper() + ' ' + gal[-2:], 
             transform=ax.transAxes, color="black", 
             fontsize=fig_prop['general']['axtext_fontsize'])
    ax_bottom.hlines(0, xmin=velocities.min() - 250, 
                     xmax=velocities.max() + 250, color="k", lw=2, ls='--')
    secax = ax.secondary_xaxis('top')
    secaxr = ax.secondary_yaxis('right')
    bottom_secax = ax_bottom.secondary_xaxis('top')
    bottom_secaxr = ax_bottom.secondary_yaxis('right')
    plt.setp(bottom_secaxr.get_yticklabels(), visible=False)
    plt.setp(bottom_secax.get_xticklabels(), visible=False)

    # Interpolating VLA data to match MeerKAT velocities within range
    velocities_interp = velocities[(velocities >= vla_min) & (velocities <= vla_max)]
    masked_prof_jy_interp = masked_prof_jy[(velocities >= vla_min) & (velocities <= vla_max)]
    spline_interp = interp1d(vla_spec[:, 0], vla_spec[:, 1], kind='linear', bounds_error=False, fill_value="extrapolate")
    vla_interp = spline_interp(velocities_interp)
    residuals = (vla_interp - masked_prof_jy_interp * 1000)

    ax_bottom.plot(velocities_interp, residuals, '-', color='k', lw=fig_prop['line_plot']['line_width'])
    secax.tick_params(labeltop=False)
    secaxr.tick_params(labelright=False)
    ax_bottom.set_xlabel(r"$\mathrm{Velocity~(km~s^{-1})}$", 
                         labelpad=fig_prop['general']['mat_xlabelpad'], 
                         fontsize=fig_prop['general']['xylabel_fontsize'])
    ax.set_ylabel(r"$\mathrm{Flux~(mJy)}$", 
                  labelpad=fig_prop['general']['mat_ylabelpad'],
                  fontsize=fig_prop['general']['xylabel_fontsize'])
    ax_bottom.text(0.03, 0.73, 'Residuals', transform=ax_bottom.transAxes, 
                   color="black", fontsize=fig_prop['general']['axtext_fontsize'])
    ax_bottom.set_ylim(ax_bottom_ylim[gal][0], ax_bottom_ylim[gal][1])
    axlist = [ax, secax, secaxr, ax_bottom, bottom_secax, bottom_secaxr]
    axprop(axlist)
    fig2 = plt.gcf()
    fig2.set_figheight(fig_prop['general']['figsize'][0])
    fig2.set_figwidth(fig_prop['general']['figsize'][1])
    figure_path = "outputs/publication/figures/"
    plt.savefig(figure_path + gal + '_global_profile.pdf', 
                bbox_inches="tight", dpi=fig_prop['general']['dpi'])

custom_font = figure_properties.FontManager(['/opt/fits-tools/src/analysis_tools/'], 'tex gyre heros')

if __name__ == '__main__':
    args = get_args()
    custom_font.setup_fonts()
    plot_global_prof(args.cube, args.mask, args.vlaspec)

# python workflow/scripts/plot_global_prof.py --cube outputs/data/hcg97/hcg97_line60_masked.pb_corr_vopt.fits --mask outputs/sofia_pbc/hcg97/hcg97_line60_masked.pb_corr_vopt_mask.fits --vlaspec data/hcg97/vla_data/hcg97_spec.txt 
