import os
import legacystamps
from analysis_tools.functions import radec2deg

hcg_position = {
        #"hcg16": {"ra": "2:09:57.8648796950", "dec": "-10:13:24.2543007855"},
        #"hcg16": {"ra": "02:09:52.5", "dec": "-10:13:20"},
        #"hcg30": {"ra": "4:36:40.1", "dec": "-2:49:21"},
        #"hcg31": {"ra": "05:01:38.3", "dec": "-04:15:25.00"},
        #"hcg90": {"ra": "22:02:06.6", "dec": "-31:58:14"},
        #"hcg91": {"ra": "22:08:54.8", "dec": "-27:47:35"},
        #'22h08m54.8s', '-27d47m35s'
        "hcg97": {"ra": "23:47:23.4", "dec": "-2:19:09.41"},
}

band = ['g', 'z', 'i']
band = ['r', 'g', 'z', 'i']

for gal in hcg_position.keys():
    if gal == "hcg31":
        size = 0.22
        pixscale = 0.35
    elif gal == "hcg91":
        size = 14.5/60
        pixscale = 0.25
    elif gal == "hcg90":
        size = 0.6
        pixscale = 0.9
    elif gal == "hcg30":
        size = 0.25
        pixscale = 0.3
    else:
        size = 0.5
        pixscale = 0.6
    radec = radec2deg(hcg_position[gal]["ra"], hcg_position[gal]["dec"])
    ra = format(round(radec["ra"], 6), '.6f')
    dec = format(round(radec["dec"], 6), '.6f')
    print(ra, dec, " RADEC")
    for filter in band:
        legacystamps.download(ra=float(ra), dec=float(dec), mode='fits', bands=filter, 
                size=size, layer='ls-dr10', pixscale=pixscale)
        filename = ("legacystamps_" + ra + '_' 
                + dec + "_ls-dr10.fits")
        os.system('mv ' + filename + " data/" + gal + "/optical/" + "center_" + filename[:-5] + "-" + filter + "filter.fits")
