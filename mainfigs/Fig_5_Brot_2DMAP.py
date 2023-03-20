from pylib import *
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature

S=loadz('./ECGOM_GA/npz/PaperExp/Comb-exp-elev-33_45.npz')
gd=read_schism_hgrid('./ECGOM_GA/grd/hgrid.gr3')
sname='/Users/kpark350/Dropbox (GaTech)/Old_Mac/Ga_tech/Journal/2nd_GulfStream/Figure/Fig_5'
cmap='jet'
lon_formatter = LongitudeFormatter(number_format='.0f', zero_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f')

scale = '110m'
states110 = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale=scale,
            facecolor='none',
            edgecolor='k')
lon_formatter = LongitudeFormatter(number_format='.2f', zero_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.2f')
###################masking########################################################################
# # inundation depth
fpt=gd.dp<0
S.ctrl[:,fpt]=S.ctrl[:,fpt]+gd.dp[fpt]
S.atf[:,fpt]=S.atf[:,fpt]+gd.dp[fpt]
S.bcadj[:,fpt]=S.bcadj[:,fpt]+gd.dp[fpt]
S.btadj[:,fpt]=S.btadj[:,fpt]+gd.dp[fpt]

fpt=S.ctrl<0.01; S.ctrl[fpt]=NaN
# fpt=S.atf<0.01; S.atf[fpt]=NaN
fpt=gd.dp<0; S.atf[:,fpt]=NaN
fpt=S.bcadj<0.01; S.bcadj[fpt]=NaN
fpt=S.btadj<0.01; S.btadj[fpt]=NaN

fpt=S.ctrl<0; S.ctrl[fpt]=NaN
# fpt=S.atf<0; S.atf[fpt]=NaN
fpt=gd.dp<0; S.atf[:,fpt]=NaN
fpt=S.bcadj<0; S.bcadj[fpt]=NaN
fpt=S.btadj<0; S.btadj[fpt]=NaN


# Land masking
# fpt=gd.dp<0
# S.ctrl[:,fpt]=NaN
# S.atf[:,fpt]=NaN
# S.bcadj[:,fpt]=NaN
# S.btadj[:,fpt]=NaN

#Gulf coast
fpt = (gd.y < 30) * (gd.x < -82);
S.ctrl[:,fpt]=NaN
S.atf[:,fpt]=NaN
S.bcadj[:,fpt]=NaN
S.btadj[:,fpt]=NaN
fpt = (gd.y < 27) * (gd.x < -81);
S.ctrl[:,fpt]=NaN
S.atf[:,fpt]=NaN
S.bcadj[:,fpt]=NaN
S.btadj[:,fpt]=NaN

# Only continental shelf
fpt=(gd.dp<=600); S.ctrl[:,~fpt]=NaN; S.atf[:,~fpt]=NaN; S.bcadj[:,~fpt]=NaN; S.btadj[:,~fpt]=NaN
# Cut Bahama
fpt=(gd.x>-79.52)*(gd.y<28); S.ctrl[:,fpt]=NaN; S.btadj[:,fpt]=NaN; S.bcadj[:,fpt]=NaN; S.atf[:,fpt]=NaN
###################################################################################################

xl=[-81.2406,-80.8215]; yl=[31.8,32.1930]; # NGA
# xl=[-81.66,-81.3956]; yl=[30.6333,30.886]; # SGA
# xl=[-83, -75];yl=[25, 35.5]; #SAB
#xl=[-81.53,-80.81]; yl=[31.75,32.28]; # Sav2 for HWMs
#xl=[-81.72,-81.42]; yl=[30.64,30.89]; # King's bay
c1=[0,45]


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
##########################################################################################


fs=datenum('2016-10-11 00:00:00'); fe=datenum('2016-10-11 00:00:00')
fpt=(S.time>=fs)*(S.time<=fe)
fpt=np.where(fpt)[0]

tiles = cimgt.GoogleTiles(style = 'satellite')
tiles_res = 10

c1=[0,40]

######################################################
###################### Plot ##########################
######################################################

c2=c1
extent = [xl[0], xl[1], yl[0], yl[1]]
figure(1, figsize=[7.3, 7.3/3])
clf()
ax=subplot(1,4,1,projection=tiles.crs)
ax.set_extent(extent,crs=ccrs.PlateCarree())
ax.add_image(tiles, tiles_res,interpolation='spline36')
ax.set_xticks(linspace(extent[0], extent[1], num=3), crs=ccrs.PlateCarree())
ax.set_yticks(linspace(extent[2], extent[3], num=4), crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
gd.plot(fmt=1,value=S.ctrl[fpt[0],:]*100,clim=c2, cmap=cmap,cb=False,zorder=10,transform = ccrs.PlateCarree())
plot(-81.1383333, 32.16555556, 'w.', ms=10, zorder=11, transform = ccrs.PlateCarree())


ax2=subplot(1,4,2,projection=tiles.crs)
ax2.set_extent(extent,crs=ccrs.PlateCarree())
ax2.add_image(tiles, tiles_res,interpolation='spline36')
ax2.set_xticks(linspace(extent[0], extent[1], num=3), crs=ccrs.PlateCarree())
ax2.set_yticks(linspace(extent[2], extent[3], num=4), crs=ccrs.PlateCarree())
ax2.xaxis.set_major_formatter(lon_formatter)
ax2.yaxis.set_major_formatter(lat_formatter)
gd.plot(fmt=1,value=S.atf[fpt[0],:]*100,clim=c2, cmap=cmap,cb=False,zorder=10,transform = ccrs.PlateCarree())
plot(-81.1383333, 32.16555556, 'w.', ms=10, zorder=11, transform = ccrs.PlateCarree())
setp(gca(), xticklabels=[], yticklabels=[])

ax3=subplot(1,4,3,projection=tiles.crs)
ax3.set_extent(extent,crs=ccrs.PlateCarree())
ax3.add_image(tiles, tiles_res,interpolation='spline36')
ax3.set_xticks(linspace(extent[0], extent[1], num=3), crs=ccrs.PlateCarree())
ax3.set_yticks(linspace(extent[2], extent[3], num=4), crs=ccrs.PlateCarree())
ax3.xaxis.set_major_formatter(lon_formatter)
ax3.yaxis.set_major_formatter(lat_formatter)
gd.plot(fmt=1,value=S.bcadj[fpt[0],:]*100,clim=c2, cmap=cmap,cb=False,zorder=10,transform = ccrs.PlateCarree())
plot(-81.1383333, 32.16555556, 'w.', ms=10, zorder=11, transform = ccrs.PlateCarree())
setp(gca(), xticklabels=[], yticklabels=[])

ax4=subplot(1,4,4,projection=tiles.crs)
ax4.set_extent(extent,crs=ccrs.PlateCarree())
ax4.add_image(tiles, tiles_res,interpolation='spline36')
ax4.set_xticks(linspace(extent[0], extent[1], num=3), crs=ccrs.PlateCarree())
ax4.set_yticks(linspace(extent[2], extent[3], num=4), crs=ccrs.PlateCarree())
ax4.xaxis.set_major_formatter(lon_formatter)
ax4.yaxis.set_major_formatter(lat_formatter)
gd.plot(fmt=1,value=S.btadj[fpt[0], :]*100,clim=c2, cmap=cmap,cb=False,zorder=10,transform = ccrs.PlateCarree())
plot(-81.1383333, 32.16555556, 'w.', ms=10, zorder=11, transform = ccrs.PlateCarree())
setp(gca(), xticklabels=[], yticklabels=[])

tight_layout()
savefig('{}'.format(sname) + '.png',dpi=900,bbox_inches='tight')
close()
