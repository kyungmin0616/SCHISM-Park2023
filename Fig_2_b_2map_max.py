from pylib import *
import cmocean
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature

# S=loadz('./ECGOM_GA/npz/PaperExp/CTRL-elev-vel-1_45.npz')
S=loadz('./ECGOM_GA_GS/npz/RUN01e-elev-hvel-1_51.npz'); S.time = S.time + datenum(2016,9,8)

sname=os.path.expanduser('~/Dropbox (GaTech)/GaTech/Journal/2nd_GulfStream/Figure_GS/Fig_2_b_newgrid')
# gd=read_schism_hgrid('./ECGOM_GA/grd/hgrid.gr3')
gd=read_schism_hgrid('./ECGOM_GA_GS/grd/01.gr3')

# cmap='nipy_spectral'
cmap='jet'


st=datenum('2016-9-18'); et=datenum('2016-10-2')
fpt=(S.time>=st)*(S.time<=et)
S.elev_ano=[]
tmean=S.elev[fpt,:].mean(axis=0)
for i in arange(0,len(S.time)):
   S.elev_ano.append(S.elev[i,:]-tmean)
S.elev_ano=array(S.elev_ano)

st=datenum('2016-10-11'); et=datenum('2016-10-22')
fpt=(S.time>=st)*(S.time<=et)
S.elev_max=amax(S.elev_ano[fpt,:],axis=0)
S.elev_max=array(S.elev_max)

fpt=gd.dp<0
S.elev_max[fpt]=S.elev_max[fpt]+gd.dp[fpt]
fpt2=(gd.dp<0)*(S.elev_max<0)
S.elev_max[fpt2]=NaN
# fpt=(gd.dp<=500); S.elev_max[~fpt]=NaN

# sname='/Users/kpark350/Downloads/2nd_Paper/Figure/fig_1_d_SGA2'
xl=[-83, -72];yl=[25, 35.5]; #SAB
# xl=[-81.38,-80.79]; yl=[31.79,32.24]; # Sav
# xl=[-81.7,-80.4]; yl=[30.8,32.28]; # GA
# xl=[-81.66,-81.3956]; yl=[30.6333,30.886]; # SGA
# xl=[-81.2406,-80.8215]; yl=[31.8,32.1930]; # NGA

# c1=[-0.4,0.4]
c1=[-0.5,0.5]

# c1=[0.2,0.4] # SGA
# c1=[0.2,0.3] # NGA

levels=linspace(c1[0],c1[1],101)
cbarlabels = linspace(levels.min(),levels.max(),9)

tiles = cimgt.GoogleTiles(style = 'satellite')
tiles_res = 9

extent = [xl[0], xl[1], yl[0], yl[1]] # ECGOM

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

figure(1, figsize=(4,4))
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
# xlabel('Longitude'); ylabel('Latitude')
# ioff()
tricontour(gd.x,gd.y,gd.dp,levels=[200, 600, 800, 1000],colors=['k','k','k','k'],linestyles=['solid','solid','solid','solid'],linewidths=1.4,transform = ccrs.PlateCarree())
gd.plot(fmt=1,value=S.elev_max,clim=c1,levels=100,cmap=cmap,cb=False,transform = ccrs.PlateCarree())
# cbar=colorbar(ticks=cbarlabels,orientation='vertical',fraction=0.05, pad=0.015)
# cbar.ax.tick_params(labelsize=18)
# cbar.ax.set_title('m',fontsize=18)
tight_layout()
savefig('{}'.format(sname) + '.png',dpi=900,bbox_inches='tight')
close()

clf()
gd.plot(fmt=1,value=S.elev_max,clim=c1,levels=100,cmap=cmap,cb=False)
tricontour(gd.x,gd.y,gd.dp,levels=[200, 600, 800, 1000],colors=['k','k','k','k'],linestyles=['solid','solid','solid','solid'],linewidths=1.4)

Xflat, Yflat, Zflat = gd.x.flatten(), gd.y.flatten(), S.elev_max.flatten()
def fmt(x, y):
    # get closest point with known data
    dist = np.linalg.norm(np.vstack([Xflat - x, Yflat - y]), axis=0)
    idx = np.argmin(dist)
    z = Zflat[idx]
    return 'x={x:.5f}  y={y:.5f}  z={z:.5f}'.format(x=x, y=y, z=z)

gca().format_coord = fmt
show()

# clf()
# gd2.plot(fmt=1 ,cmap='jet')
# tricontour(gd.x,gd.y,gd.dp,levels=[200, 600, 800, 1000],colors=['k','k','k','k'],linestyles=['solid','solid','solid','solid'],linewidths=1.4)