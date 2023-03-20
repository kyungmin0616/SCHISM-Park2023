from pylib import *
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
from pyproj import Geod

geod = Geod(ellps='WGS84')

##########################################################################################
# S=loadz('./ECGOM_GA/npz/PaperExp/CTRL-elev-vel-1_45.npz')
# gd=read_schism_hgrid('./ECGOM_GA/grd/hgrid.gr3')
bp=read_shapefile_data('./ECGOM_GA/SHP/Paper_SAB_COAST.shp')
# sname='/Users/kpark350/Dropbox (GaTech)/Old_Mac/Ga_tech/Journal/2nd_GulfStream/Figure/Fig_2_a'

S=loadz('./ECGOM_GA_GS/npz/RUN01e-elev-hvel-1_51.npz'); S.time = S.time + datenum(2016,9,8)
sname=os.path.expanduser('~/Dropbox (GaTech)/GaTech/Journal/2nd_GulfStream/Figure/Fig_2_a_newgrid')
# gd=read_schism_hgrid('./ECGOM_GA/grd/hgrid.gr3')
gd=read_schism_hgrid('./ECGOM_GA_GS/grd/01.gr3')

st=datenum(2016,9,22)
se=datenum(2016,10,22)
fpt=(S.time>=st)*(S.time<=se); S.time=S.time[fpt]; S.elev=S.elev[fpt]

cmap='nipy_spectral'
cmap='jet'

levels=linspace(-0.5,0.5,101)
cbarlabels = linspace(levels.min(),levels.max(),5)
##########################################################################################

st=datenum('2016-9-18'); et=datenum('2016-10-2')
fpt=(S.time>=st)*(S.time<=et)
S.elev_ano=[]
tmean=S.elev[fpt].mean(axis=0)
for i in arange(0,len(S.time)):
   S.elev_ano.append(S.elev[i,:]-tmean)
S.elev_ano=array(S.elev_ano)
st=datenum(2016,9,22)
se=datenum(2016,10,22)

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

px,py=bp.xy.T
px = px[~isnan(px.astype(float64))].astype(float64); py = py[~isnan(py.astype(float64))].astype(float64)
px,py=proj_pts(px,py,'epsg:26918','epsg:4326')

sindp=near_pts(c_[px,py],c_[gd.x,gd.y])
sindp = [i for n, i in enumerate(sindp) if i not in sindp[:n]] # delete dupulicate remain the oder of elements

hx=gd.x[sindp]
hy=gd.y[sindp]
helev1=S.elev[:,sindp]
helev1=helev1.transpose()
helev1_ano=S.elev_ano[:,sindp]
helev1_ano=helev1_ano.transpose()
htime=S.time

hdist=[0]

for i in arange(len(hx)):
    if i == len(hx)-1: continue
    azimuth1, azimuth2, distance = geod.inv(hx[i], hy[i], hx[i+1], hy[i+1])
    hdist.append(hdist[i]+distance)
xv,yv=meshgrid(htime,hdist)

lc_NOAA=array([0, 219427.93280598585, 494839.81714877044, 661747.3108285896, 789227.0178567843, 885159.510347716, 1068278.0446788708, 1180779.5138946348])

#### plot Hovmoller diagram
figure(1, figsize=[7.3, 7.3/3])
clf()
contourf(xv, yv/1000000, helev1_ano, levels=levels,cmap=cmap,extend='both')
plot(xlim(),[498.591495804421e-3, 498.591495804421e-3],'k',linestyle='dotted',lw=2);
plot(xlim(),[658.0023447720065e-3, 658.0023447720065e-3],'k',linestyle='dotted',lw=2)
plot(xlim(),[965.3986439407874e-3, 965.3986439407874e-3],'k',linestyle='dotted',lw=2)

for i in lc_NOAA:
    tdist = abs(hdist - i)
    nn=argmin(tdist)
    nn2 =where(helev1_ano[nn,:] == helev1_ano[nn,:].max())[0][0]
    plot(htime[nn2],i/1000000,'k.',ms=10)
    # print(num2date(htime[nn2]))
    # print(i/1000000)

plot([datenum('2016-10-4'), datenum('2016-10-4')],ylim(),'k',linestyle='dotted',lw=2)
# plot([datenum('2016-10-7'), datenum('2016-10-7')],ylim(),'k',linestyle='dotted',lw=2)
# plot([datenum('2016-10-8'), datenum('2016-10-8')],ylim(),'k',linestyle='dotted',lw=2)
# plot([datenum('2016-10-9'), datenum('2016-10-9')],ylim(),'r',linestyle='dotted',lw=2)

plot([datenum('2016-10-10'), datenum('2016-10-10')],ylim(),'k',linestyle='dotted',lw=2)

# plot([datenum('2016-10-12'), datenum('2016-10-12')],ylim(),'k',linestyle='dotted',lw=2)
# plot([datenum('2016-10-14'), datenum('2016-10-14')],ylim(),'k',linestyle='dotted',lw=2)
# plot([datenum('2016-10-13'), datenum('2016-10-13')],ylim(),'k',linestyle='dotted',lw=2)

# yvz=yv/1000000
# Xflat, Yflat, Zflat = xv.flatten(), yvz.flatten(), helev1_ano.flatten()
# def fmt(x, y):
#     # get closest point with known data
#     dist = np.linalg.norm(np.vstack([Xflat - x, Yflat - y]), axis=0)
#     idx = np.argmin(dist)
#     z = Zflat[idx]
#     return 'x={x:.5f}  y={y:.5f}  z={z:.5f}'.format(x=x, y=y, z=z)
#
# gca().format_coord = fmt
# show()

# cbar=colorbar(ticks=cbarlabels,orientation='horizontal',fraction=0.05, pad=0.1)
xts, xls = get_xtick(fmt=2, xts=[*arange(st,se+5,5)], str='%m-%d')
setp(gca(), xticks=xts, xticklabels=xls, xlim=[st, se],yticks=arange(0,1.4,0.2),ylim=[-0.01-0.04,1.2+0.04])
ylabel('Distance (km)$\cdot 10^3$')
xlabel('Date (2016-)')

tight_layout()
savefig('{}'.format(sname) + '.png',dpi=900,bbox_inches='tight')
close()
#
# ######### find boundary
# conus_proj = ccrs.LambertConformal(central_longitude=-96,central_latitude=39.0)
# fig = figure(figsize=(10,8))
# ax = fig.add_subplot(1,1,1,projection=conus_proj)
# ax.set_extent([-83,-72,22,40])
# ax.add_feature(cfeature.STATES, edgecolor='black', zorder=10)
# for i in arange(0,len(hx)):
#     plot(hx[i],hy[i],'r.',transform=ccrs.PlateCarree())
#     at_x, at_y = ax.projection.transform_point(hx[i], hy[i],src_crs=ccrs.PlateCarree())
#     annotate(hdist[i], xy=(at_x, at_y))