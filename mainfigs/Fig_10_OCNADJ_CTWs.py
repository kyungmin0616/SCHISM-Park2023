from pylib import *
import cmocean
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
S=loadz('./ECGOM_GA/npz/PaperExp/Comb-exp-elev-33_45.npz')
gd=read_schism_hgrid('./ECGOM_GA/grd/hgrid.gr3')
sname='/Users/kpark350/Downloads/2nd_Paper/Figure/CTS_OCNADJ-ver3'
cmap='jet'

###################masking########################################################################
# inundation depth
fpt=gd.dp<0
S.ctrl[:,fpt]=S.ctrl[:,fpt]+gd.dp[fpt]
S.atf[:,fpt]=S.atf[:,fpt]+gd.dp[fpt]
S.bcadj[:,fpt]=S.bcadj[:,fpt]+gd.dp[fpt]
S.btadj[:,fpt]=S.btadj[:,fpt]+gd.dp[fpt]

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
# Inland masking
fpt=gd.dp<0
S.ctrl[:,fpt]=S.ctrl[:,fpt]+gd.dp[fpt]; S.atf[:,fpt]=S.atf[:,fpt]+gd.dp[fpt]; S.bcadj[:,fpt]=S.bcadj[:,fpt]+gd.dp[fpt]; S.btadj[:,fpt]=S.btadj[:,fpt]+gd.dp[fpt];
###################################################################################################



xl=[-83, -75];yl=[25, 35.5]; #SAB
#xl=[-81.53,-80.81]; yl=[31.75,32.28]; # Sav2 for HWMs
#xl=[-81.72,-81.42]; yl=[30.64,30.89]; # King's bay

## CTW animation
c2=[0,45]
clf()
gz = gd.plot(fmt=1, value=S.ctrl[0,:], clim=c2, cmap=cmap, cb=False)
cbarlabels = linspace(c2[0], c2[1], 10)
colorbar(ticks=cbarlabels, orientation='horizontal')


figure(1, figsize=[9, 9])
clf()
tiles = cimgt.GoogleTiles(style = 'satellite')
tiles_res = 11
extent = [xl[0], xl[1], yl[0], yl[1]] # ECGOM
ax = plt.axes(projection=tiles.crs)
ax.set_extent(extent, crs=ccrs.PlateCarree())
ax.add_image(tiles, tiles_res,interpolation='spline36')
ax.set_xticks(arange(extent[0],extent[1],abs(extent[0]-extent[1])/5), crs=ccrs.PlateCarree())
ax.set_yticks(arange(extent[2],extent[3],abs(extent[2]-extent[3])/5), crs=ccrs.PlateCarree())
ax.set_xticks(linspace(extent[0],extent[1],num=5), crs=ccrs.PlateCarree())
ax.set_yticks(linspace(extent[2],extent[3],num=6), crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.1f',zero_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.1f')
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ioff()
for nn, datelist in enumerate(S.time):
    nn=nn+10
    stime=S.time[nn]
    fpt = (gd.dp < 0) * (S.bcadj[nn,:] < 0); S.bcadj[nn, fpt] = NaN
    fpt = (gd.dp < 0) * (S.btadj[nn,:] < 0); S.btadj[nn, fpt] = NaN
    tmp=(S.bcadj[nn,:]+S.btadj[nn,:])*100
    gz=gd.plot(fmt=1,value=tmp,clim=c2, cmap=cmap,cb=False,transform = ccrs.PlateCarree())
    # cbar = colorbar(ticks=arange(0,60+10,10), orientation='vertical')
    fpt = (gd.dp>70); tmp[fpt]=NaN
    fpt = (gd.dp<5); tmp[fpt]=NaN
    fpt = (gd.dp>350); tmp[fpt]=NaN
    fpt = (gd.y<30)*(gd.x<-82); tmp[fpt]=NaN
    fpt = (gd.y<27)*(gd.x<-81); tmp[fpt]=NaN
    CS=gd.plot(fmt=2,value=tmp,levels=[34,36,38,40,42,44,46,48,50],colors=['k','k','k'], linestyles=['solid','solid','solid'],linewidths=0.8,transform = ccrs.PlateCarree())
    # clabel(CS, inline=1,inline_spacing=-8, fontsize=10)
    title(num2date(stime))
    # xlim(xl);ylim(yl)
    # xticks(arange(xl[0], xl[1], 2));
    # yticks(arange(yl[0], yl[1], 2))
    # gcf().tight_layout()
    savefig('{}/{}'.format(sname,num2date(stime).strftime('%Y%m%d%H%M%S'))+'.png')
    for coll in gz.collections: coll.remove()
    for coll in CS.collections: coll.remove()

    # close()


SMALL_SIZE = 21
MEDIUM_SIZE = 13
BIGGER_SIZE = 13
rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
rc('font', size=MEDIUM_SIZE)  # controls default text sizes
rc('legend', fontsize=MEDIUM_SIZE)  # legend fontsize
rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title