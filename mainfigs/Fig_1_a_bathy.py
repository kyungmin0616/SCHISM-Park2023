from pylib import *
import cmocean
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter

##########################################################################################
sname='~/Dropbox (GaTech)/Old_Mac/Ga_tech/Journal/2nd_GulfStream/Figure/Fig_1_a'
gd=read_schism_hgrid('./ECGOM_GA/grd/hgrid.gr3')
bp=read_schism_bpfile('./ECGOM_GA/station.in')
# bp2=read_shapefile_data('./ECGOM_GA/SHP/SAB_COAST.shp')
bp3=read_shapefile_data('./obs/NOAA/Matthew/best_track/al142016_pts.shp')
# px,py=bp2.xy.T
# px = px[~isnan(px)]; py = py[~isnan(py)]
# px,py=proj_pts(px,py,'epsg:26918','epsg:4326')
px2,py2=array(bp3.xy.T,dtype=float)
px2 = px2[~isnan(px2)]; py2 = py2[~isnan(py2)]
# noaast=array([5,6,7,8,10,11,12,13,15,16,17,19,20,21,22,23,24,26,27,28,30,31,32,1,2,78,77,76,75,74,73,71,67,68,69,66,64,61,62,59,60,58,57,55,56,54,53,51,50,48,47,45,44])
# noaast=array([5,6,7,8,10,11,12,13,15,16,17,19,20,21,22,23,24,26,27,28,30,31,32])
noaast= array([7,8,10,11,12,13,15,16]), # 16 ,17 need to be included in KMP

##########################################################################################
SMALL_SIZE = 8
MEDIUM_SIZE = 8
BIGGER_SIZE = 8
rc('font', family='Helvetica')
rc('font', size=SMALL_SIZE)  # controls default text sizes
rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
#################### google map
extent = [-83.5, -72, 24.5, 35.5] # SAB

tiles = cimgt.GoogleTiles(style = 'satellite')
tiles_res = 8
states110 = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='110m',
            facecolor='none',
            edgecolor='w')

figure(1, figsize=[3.5, 3.5])
clf()
ax = plt.axes(projection=tiles.crs)
ax.set_extent(extent, crs=ccrs.PlateCarree());
ax.add_image(tiles, tiles_res,interpolation='spline36')
ax.add_feature(states110, zorder=1, linewidth=1.5)
ax.set_xticks(linspace(extent[0],extent[1],num=5), crs=ccrs.PlateCarree())
ax.set_yticks(linspace(extent[2],extent[3],num=6), crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.1f',zero_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.1f')
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
# xlabel('Longitude'); ylabel('Latitude')
################################### NOAA ###################################
# Hurricane
for nn,sssc in enumerate(bp3.attvalue[10]):
    if sssc=='Hurricane1': color='yellow'
    elif sssc=='Hurricane2': color='orange'
    elif sssc == 'Tropical Storm': color = 'yellowgreen'
    elif sssc == 'Extratropical Cyclone': color = 'grey'
    elif sssc=='Hurricane3': color='orangered'
    elif sssc=='Hurricane4': color='fuchsia'
    elif sssc=='Hurricane5': color='purple'
    else: color='k'
    if nn==49: continue
    else: plot([px2[nn],px2[nn+1]],[py2[nn],py2[nn+1]],marker='o',lw=2,color=color,transform = ccrs.PlateCarree())

# Tide gauges
tx=bp.x[noaast[0]-1]; ty=bp.y[noaast[0]-1]
for i in arange(len(tx)):
    plot(tx[i], ty[i], marker='s', linestyle='None', color='white', transform=ccrs.PlateCarree())

# for i in noaast:
#     plot(bp.x[i - 1], bp.y[i - 1], marker='s', linestyle='None', color='black', transform=ccrs.PlateCarree())

tight_layout()
savefig('{}'.format(sname) + '.png',dpi=900,bbox_inches='tight')
close()
# #
# ######### time with location
# for i,time in enumerate(bp3.attvalue[1]):
#     plot(px2[i],py2[i],'r.',transform=ccrs.PlateCarree())
#     at_x, at_y = ax.projection.transform_point(px2[i], py2[i],src_crs=ccrs.PlateCarree())
#     annotate(time, xy=(at_x, at_y))


