from pylib import *
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt

S=loadz('./ECGOM_GA/npz/Allforing-w-tide-elev-1_55.npz')
sname='/Users/park075/Dropbox (GaTech)/Old_Mac/Ga_tech/Journal/2nd_GulfStream/Figure/Fig_3_a'
gd=read_schism_hgrid('./ECGOM_GA/grd/hgrid.gr3')

####################################################################################
st=datenum('2016-9-18'); et=datenum('2016-10-2')
fpt=(S.time>=st)*(S.time<=et)
S.elev_max_normal=[]
for i in arange(len(gd.x)):
    S.elev_max_normal.append(max(S.elev[fpt,i]))
S.elev_max_normal=array(S.elev_max_normal)

st=datenum('2016-10-11'); et=datenum('2016-10-22')
fpt=(S.time>=st)*(S.time<=et)
S.elev_max_post=[]
for i in arange(len(gd.x)):
    S.elev_max_post.append(max(S.elev[fpt,i]))
S.elev_max_post=array(S.elev_max_post)
####################################################################################
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
####################################################################################
fpt=gd.dp<0
S.elev_max_normal[fpt]=S.elev_max_normal[fpt]+gd.dp[fpt]
S.elev_max_post[fpt]=S.elev_max_post[fpt]+gd.dp[fpt]
fpt2=(gd.dp<0)*(S.elev_max_normal<=0)
S.elev_max_normal[fpt2]=0
fpt2=(gd.dp<0)*(S.elev_max_post<=0)
S.elev_max_post[fpt2]=0
fpt=gd.dp>=0
S.elev_max_normal[fpt]=NaN
S.elev_max_post[fpt]=NaN
####################################################################################
diff=S.elev_max_post-S.elev_max_normal
fpt2=(gd.dp<0)*(diff<=0)
diff[fpt2]=NaN
####################################################################################

c1=[-0.4,0.4]

# xl=[-81.66,-81.3956]; yl=[30.6333,30.886]; # SGA
xl=[-81.2406,-80.8215]; yl=[31.8,32.1930]; # NGA

tiles = cimgt.GoogleTiles(style = 'satellite') #style (optional) – The style for the Google Maps tiles. One of ‘street’, ‘satellite’, ‘terrain’, and ‘only_streets’. Defaults to ‘street’.
extent = [xl[0], xl[1], yl[0], yl[1]]
tiles_res = 12
figure(1, figsize=(3.5,3.5))
clf()
ax = plt.axes(projection=tiles.crs)
ax.set_extent(extent, crs=ccrs.PlateCarree());
ax.add_image(tiles, tiles_res,interpolation='spline36')
ax.set_xticks(arange(extent[0],extent[1],abs(extent[0]-extent[1])/5), crs=ccrs.PlateCarree())
ax.set_yticks(arange(extent[2],extent[3],abs(extent[2]-extent[3])/5), crs=ccrs.PlateCarree())
ax.set_xticks(linspace(extent[0],extent[1],num=5), crs=ccrs.PlateCarree())
ax.set_yticks(linspace(extent[2],extent[3],num=6), crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(number_format='.2f',zero_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.2f')
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
gd.plot(fmt=1,value=diff,clim=c1,levels=100,cmap='bwr',cb=False,transform = ccrs.PlateCarree())
# cbar=colorbar(ticks=cbarlabels,orientation='horizontal',fraction=0.05, pad=0.2)
# cbar.ax.tick_params(labelsize=18)
# cbar.ax.set_title('m',fontsize=18)
tight_layout()
savefig('{}'.format(sname) + '.png',dpi=900,bbox_inches='tight')
close()


gd.plot(fmt=1,value=diff,clim=c1,levels=100,cmap='bwr',cb=False)
