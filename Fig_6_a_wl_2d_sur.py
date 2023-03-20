from pylib import *
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import cmocean
import matplotlib.gridspec as gridspec



###################################################################################
cmap='nipy_spectral'
stt=datenum('2016-10-11'); edt=datenum('2016-10-22')
gd=read_schism_hgrid('./ECGOM_GA/grd/hgrid_utm.gr3')
gd2=read_schism_hgrid('./ECGOM_GA/grd/hgrid.gr3')
###################################################################################
# V=loadz('./ECGOM_GA/npz/PaperExp/CTRL-elev-vel-33_45.npz')
# V=loadz('./ECGOM_GA/npz/PaperExp/CTRL-WO-AIR-elev-vel-33_45.npz')
V=loadz('./ECGOM_GA/npz/PaperExp/BCADJ-elev-vel-33_45.npz')

# V2=loadz('./ECGOM_GA/obs/npz/AVISO.npz')
# V2=loadz('./ECGOM_GA/obs/npz/HYCOM.npz')

###################################################################################
omega=7.29e-5
f=2*omega*sin(deg2rad(gd2.y))
###################################################################################
st=datenum(2016,10,11)
se=datenum(2016,10,22)
fpt=(V.time>=st)*(V.time<=se);
V.elev=V.elev[fpt,:]; V.time=V.time[fpt];
V.u=V.u[fpt,:]; V.v=V.v[fpt,:];

gv=[]
for i in arange(len(V.time)):
    [dhdx, dhdy, dhdxy] = gd.compute_gradient(value=V.elev[i, :])
    gv.append((9.8/f)*(dhdx))
gv=array(gv)
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

dx=3000
tnames=['NC','GA','FL']
# txys = array([[[-77.779926, 34.206790], [-74.716270,34.206790]], # NC
#       [[-79.776090, 32.701826],[-74.716270, 32.701826]], # SC
#       [[-81.122218, 31.403865],[-74.716270, 31.403865]], # GA
#       [[-80.641043, 28.706134],[-74.716270, 28.706134]]]) # FL
# # ver2: FL location revised (up more)
# txys = array([[[-77.779926, 34.206790], [-74.716270,34.206790]], # NC
#       [[-79.776090, 32.701826],[-74.716270, 32.701826]], # SC
#       [[-81.122218, 31.403865],[-74.716270, 31.403865]], # GA
#       [[-81.003, 29.267],[-74.716270, 29.267]]]) # FL

# # ver3: FL location revised (up more) and SC location revised (up more)
# txys = array([[[-77.779926, 34.206790], [-74.716270,34.206790]], # NC
#       [[-79.116597, 33.090667],[-74.716270, 33.090667]], # SC
#       [[-81.122218, 31.403865],[-74.716270, 31.403865]], # GA
#       [[-81.003, 29.267],[-74.716270, 29.267]]]) # FL

# # ver3: FL location revised (up more), without SC, GA location revised (up more)
txys = array([[[-77.779926, 34.206790], [-74.716270,34.206790]], # NC
      [[-81.031962, 31.653739],[-74.716270, 31.653739]], # GA
      [[-81.003, 29.267],[-74.716270, 29.267]]]) # FL

# ver4: FL location revised (up more), without SC, GA location revised (32N to be comparable with HYCOM)
# txys = array([[[-77.779926, 34.206790], [-74.716270,34.206790]], # NC
#       [[-80.752394, 32.0],[-74.716270, 32.0]], # GA
#       [[-81.003, 29.267],[-74.716270, 29.267]]]) # FL

T=zdata(); T.x=[]; T.y=[]; T.z=[]; T.dist=[]; T.dx=[]; T.xv=[]; T.yv=[]; T.elev=[]; T.tnames=[]; T.gv=[]; T.v=[]

for m,tname in enumerate(tnames):
    sxy=txys[m]; sx,sy=sxy.T; sx,sy=proj_pts(sx,sy,'epsg:4326','epsg:26918'); sind=gd.compute_acor(array([[sx[0],sy[0]],[sx[1],sy[1]]]))[0]
    if sum(sind == -1) != 0: sys.exit('pts outside of domain: {},{} {}'.format(m, sxy, sind))
    tdist = abs(diff(sx) + 1j * diff(sy))[0];
    npt = int(around(tdist / dx)) + 1;
    dx = tdist / npt
    sxi = linspace(*sx, npt);
    syi = linspace(*sy, npt);
    dist = linspace(0, tdist, npt)
    sindp = near_pts(c_[sxi, syi], c_[gd.x, gd.y])
    szi= gd.dp[sindp]
    # pie, pip, pacor = gd.compute_acor(c_[sxi, syi])
    # szi = (gd.dp[pip] * pacor).sum(axis=1);
    elevi=[];gvi=[]; vi=[]
    xv, yv = np.meshgrid(dist,V.time)
    for i,time in enumerate(V.time):
        # ei = (V.elev[i][pip] * pacor).sum(axis=1)
        gvii = gv[i,sindp]
        vii= V.v[i,sindp]
        ei = V.elev[i,sindp]
        elevi.append(ei); gvi.append(gvii); vi.append(vii)
    xv = array(xv); yv = array(yv); elevi = array(elevi); gvi = array(gvi); vi = array(vi)
    T.x.append(sxi); T.y.append(syi); T.z.append(szi); T.dist.append(dist); T.dx.append(dx)
    T.elev.append(elevi); T.xv.append(xv); T.yv.append(yv); T.tnames.append(tname); T.gv.append(gvi); T.v.append(vi)

lfcutoff=17/24
T.elev_sm=[]
for i in arange(len(T.tnames)):
    eiss=[]
    for j in arange(shape(T.elev[i])[1]):
        eis=lpfilt(T.elev[i][:,j], 1/24,lfcutoff)
        eiss.append(eis)
    eiss=transpose(eiss); eiss=array(eiss)
    T.elev_sm.append(eiss)

################################################### 2D for low frequency
stt = datenum('2016-10-12 12:00:00');
edt = datenum('2016-10-21')
tidx = (stt < V.time) * (V.time < edt)
bdgv=0.1
bddp=600
c1=linspace(11,23,101)
levels = linspace(c1[0], c1[-1], 101)
cbarlabels = arange(levels.min(), levels.max()+2, 2)
levels2=arange(16,30,1)
fig = plt.figure(1,figsize=[3.2,4.1])
clf()
for i, tname in enumerate(T.tnames):
    # fpt = (T.z[i] < dep[1]) * (T.gv[i].mean(axis=0) < 0.3)  # Elevation
    fpt = (T.z[i] < bddp) * (T.gv[i].mean(axis=0) < bdgv)  # ver 2 Elevation

    ax = fig.add_subplot(3, 1, i+1)
    CS=ax.contour(T.xv[i][:,fpt][tidx]/1000, T.yv[i][:,fpt][tidx], T.elev_sm[i][:,fpt][tidx]*100, levels=levels2, colors='k', linestyles='-', linewidths=1)
    clabel(CS,fmt='%2.2f', colors='k', fontsize=8)
    map2=ax.contourf(T.xv[i][:,fpt][tidx]/1000, T.yv[i][:,fpt][tidx], T.elev_sm[i][:,fpt][tidx]*100, levels=c1, cmap=cmap, extend='both')
    xts, xls = get_xtick(fmt=2, xts=arange(stt, edt + 2, 2), str='%d')
    setp(gca(), yticks=xts, yticklabels=xls, ylim=[stt, edt],xlim=[0,145])
    # cbar = colorbar(map2,ticks=cbarlabels, orientation='horizontal', fraction=0.05, pad=0.2)
    # cbar.ax.tick_params(labelsize=10)
    if i<=1:setp(gca(), xticklabels=[])
    ax.set_ylabel('Date (2016-10-)')
    if i >= 2: ax.set_xlabel('Distance (km)')

fig.tight_layout()
savefig('/Users/kpark350/Dropbox (GaTech)/Old_Mac/Ga_tech/Journal/2nd_GulfStream/Figure/Fig_6_a.png',dpi=900,bbox_inches='tight')
close()

############################################ Check transects
# xl=[-1.99e5, 5.61e5];yl=[2.753e6, 3.947e6]; #SAB
xl=[-83, -75];yl=[25, 35.5]; #SAB
bdgv=0.1
colors='rrrr'
# masking land part for velocity field
fpt=gd2.dp<5 *(29.772<gd2.y)*(gd2.y<33.1287); V.v[:,fpt]=NaN
c1=[0,1]
levels = linspace(c1[0], c1[-1], 101)
cbarlabels = arange(levels.min(), levels.max()+0.2, 0.2)

figure(2, figsize=[3.5, 5])
clf()
scale = '110m'
states110 = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale=scale,
            facecolor='none',
            edgecolor='k')
extent = [xl[0], xl[1], yl[0], yl[1]]
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(extent);
# ax.coastlines()
ax.add_feature(states110, zorder=1, linewidth=1.5)
# ax.add_feature(cart.feature.LAND, zorder=100, edgecolor='k')
ax.set_xticks(linspace(extent[0], extent[1], 5))
ax.set_yticks(linspace(extent[2], extent[3], 4))
lon_formatter = LongitudeFormatter(number_format='.1f',zero_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.1f')
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

tricontour(gd2.x, gd2.y, gd2.dp, levels=[200,600, 800, 1000], colors=['w', 'w', 'w', 'w'],  #colors=['w', '', 'k', 'k']
                linestyles=['solid', 'solid', 'solid', 'solid'], linewidths=1.4)
# gd.plot(fmt=1, value=gv.mean(axis=0), clim=[0,1.], cmap=cmocean.cm.thermal)
gd2.plot(fmt=1, value=V.v.mean(axis=0), clim=c1, cb=False, cmap=cmocean.cm.thermal)
# CS = gd2.plot(fmt=2, value=V.elev.mean(axis=0), levels=[0.2,0.3,0.4], colors=['w', 'w', 'w'], linestyles=['solid', 'solid', 'solid'],
#              linewidths=0.8)
# clabel(CS, inline=1, inline_spacing=-8, fontsize=8)
# colorbar(ticks=cbarlabels, orientation='horizontal')

for i, tname in enumerate(T.tnames):

    fpt = (T.z[i] < bddp) * (T.gv[i].mean(axis=0) < bdgv) # Elevation
    sx,sy=proj_pts(T.x[i][fpt], T.y[i][fpt], 'epsg:26918','epsg:4326')
    plot(sx, sy, linestyle='dashed', color=colors[i])

    tmp=T.v[i].mean(axis=0);
    idx=where(tmp == max(tmp))[0][0]
    # if i==1: idx=80
    fpt=(T.gv[i].mean(axis=0) > bdgv)*(T.x[i]<T.x[i][idx])  # Geo velocity
    sx,sy=proj_pts(T.x[i][fpt], T.y[i][fpt], 'epsg:26918','epsg:4326')
    plot(sx,sy,color=colors[i])

    plot(sx.min(),sy.min(),"^",color=colors[i], markersize=7) # Maximum Geo velocity point
    plot(sx.max(),sy.max(),"o",color=colors[i], markersize=7) # First Geo velocity point
xlim(xl);ylim(yl)

fig.tight_layout()
savefig('/Users/kpark350/Dropbox (GaTech)/Old_Mac/Ga_tech/Journal/2nd_GulfStream/Figure/Fig_6_b.png',dpi=900,bbox_inches='tight')
close()

################################################### Elevation vs Reproduced elev using Geostrophic balance

bdgv=0.1
fig = plt.figure(1,figsize=[7.3,4])
clf()
for i, tname in enumerate(T.tnames):

    ### index
    fpt = (T.z[i] < bddp) * (T.gv[i].mean(axis=0) < bdgv) # Elevation
    tmp=T.gv[i].mean(axis=0);
    idx=where(tmp == max(tmp))[0][0]
    fpt2=(T.gv[i].mean(axis=0) > bdgv)*(T.x[i]<T.x[i][idx])  # Geo velocity
    sx, sy = proj_pts(T.x[i], T.y[i], 'epsg:26918','epsg:4326');
    ####################################################################################
    pel = T.elev[i][:, fpt2][:, 0];
    pgv = T.v[i][:, fpt2].mean(axis=1)
    gf=2*omega*sin(deg2rad(sy[fpt2][-1]))
    gdist=(T.dist[i][fpt2][-1]-T.dist[i][fpt2][0])
    gelev = T.elev[i][:, fpt2][:, -1] -(T.v[i][:, fpt2].mean(axis=1)*gf*gdist/9.8)
    ####################################################################################
    if i==0: lfcutoff=17/24
    if i == 0: lfcutoff = 17 / 24
    if i == 0: lfcutoff = 17 / 24
    subplot(2,3,i+1)
    pel = lpfilt(pel,1/24,lfcutoff) # low-pass
    gelev = lpfilt(gelev,1/24,lfcutoff) # low-pass
    stt = datenum('2016-10-14'); edt = datenum('2016-10-21')
    fpt=(stt<V.time)*(V.time<edt)
    # plot(V.time[fpt],(pel[fpt]-pel[fpt][0])*100,'k',lw=2)
    # plot(V.time[fpt],(gelev[fpt]-gelev[fpt][0])*100, 'r',lw=2)
    plot(V.time[fpt],(pel[fpt])*100,'k',lw=2)
    plot(V.time[fpt],(gelev[fpt])*100, 'r',lw=2)
    st = get_stat(pel[fpt], gelev[fpt])
    # legend(['Model','Reconstruction'])
    xts, xls = get_xtick(fmt=2, xts=arange(stt, edt + 1, 1), str='%d')
    setp(gca(), xticks=xts,xticklabels=xls, xlim=[stt, edt-1])
    setp(gca(), xticklabels=[])
    if i==0: ylabel('Water level (cm)')
    # title('{}, R: {:.2f}'.format(T.tnames[i],st.R))
    # if i==0: ylabel('Water level (m)');
    if i==0: ylim([16,20])
    if i==1: ylim([17,24])
    if i==2: ylim([19,23])

    # xlabel('Date (2016-10-)');
    # ylabel('Water level (cm)');

    gca().xaxis.grid('on')
    gca().yaxis.grid('on')

    subplot(2, 3, i + 4)
    # ax2.plot(V.time,pgv,'g')
    # pgv = smooth(pgv, 36) # smooth
    pgv = lpfilt(pgv,1/24,lfcutoff) # low-pass
    # ax2.plot(V.time, pgv[fpt]-pgv[fpt].mean(axis=0), 'g',lw=3)
    # plot(V.time[fpt], (pgv[fpt]-pgv[fpt][0])*100, 'g',lw=3)
    plot(V.time[fpt], (pgv[fpt])*100, 'g',lw=3)

    # if i==0: ylabel('Velocity (m/s)');
    if i==0: setp(gca(), yticks=arange(52,77+5,5), ylim=[52,77])
    if i==1: ylim([35,45])
    if i==2: ylim([57.5,67.5])
    # CTRL-WO-AIR
    # ylim([-0.12, 0.4])
    # legend(['Latitudinal velocity'])
    st = get_stat(pel[fpt], pgv[fpt])
    # title('{}, R: {:.2f}'.format(T.tnames[i],st.R))
    xts, xls = get_xtick(fmt=2, xts=arange(stt, edt + 1, 1), str='%d')
    setp(gca(), xticks=xts,xticklabels=xls, xlim=[stt, edt-1])
    xlabel('Date (2016-10-)');
    if i==0: ylabel('Velocity (cm/s)')
    # ylabel('Velocity (cm/s)');

    gca().xaxis.grid('on')
    gca().yaxis.grid('on')
    # colorbar(surf)
    # plt.show()

fig.tight_layout()
savefig('/Users/kpark350/Dropbox (GaTech)/Old_Mac/Ga_tech/Journal/2nd_GulfStream/Figure/Fig_7.png',dpi=900,bbox_inches='tight')
close()

