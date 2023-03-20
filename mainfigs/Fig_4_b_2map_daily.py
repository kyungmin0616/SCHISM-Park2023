from pylib import *
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter

S=loadz('./ECGOM_GA/npz/PaperExp/Comb-exp-elev-33_45.npz')
gd=read_schism_hgrid('./ECGOM_GA/grd/hgrid.gr3')
sname='/Users/kpark350/Dropbox (GaTech)/Old_Mac/Ga_tech/Journal/2nd_GulfStream/Figure/Fig_4_b'
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

###################masking########################################################################
# # inundation depth
# fpt=gd.dp<0
# S.ctrl[:,fpt]=S.ctrl[:,fpt]+gd.dp[fpt]
# S.atf[:,fpt]=S.atf[:,fpt]+gd.dp[fpt]
# S.bcadj[:,fpt]=S.bcadj[:,fpt]+gd.dp[fpt]
# S.btadj[:,fpt]=S.btadj[:,fpt]+gd.dp[fpt]

# Land masking
fpt=gd.dp<0
S.ctrl[:,fpt]=NaN
S.atf[:,fpt]=NaN
S.bcadj[:,fpt]=NaN
S.btadj[:,fpt]=NaN

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

xl=[-83, -75];yl=[25, 35.5]; #SAB
#xl=[-81.53,-80.81]; yl=[31.75,32.28]; # Sav2 for HWMs
#xl=[-81.72,-81.42]; yl=[30.64,30.89]; # King's bay
c1=[0,60]

##########################################################################################
SMALL_SIZE = 6.5
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

## Magnitude - selected date
datelists=array([datenum('2016-10-13 00:00:00'), datenum('2016-10-16 00:00:00'), datenum('2016-10-19 00:00:00')])
c2=c1
fs=6.5
extent = [xl[0], xl[1], yl[0], yl[1]]

figure(1, figsize=[3.5, 5.5])
clf()
for figi, datelist in enumerate(datelists):
    nn=where(S.time==datelists[figi])[0][0]
    stime=S.time[nn]

    exec("ax{}=subplot(4,3,{},projection=ccrs.PlateCarree())".format(figi,figi+1))
    exec("ax{}.set_extent(extent)".format(figi))
    exec("ax{}.coastlines()".format(figi))
    exec("ax{}.add_feature(states110, zorder=0, linewidth=1)".format(figi))
    exec("ax{}.set_xticks(linspace(extent[0], extent[1], num=3), crs=ccrs.PlateCarree())".format(figi))
    exec("ax{}.set_yticks(linspace(extent[2], extent[3], num=4), crs=ccrs.PlateCarree())".format(figi))
    exec("ax{}.xaxis.set_major_formatter(lon_formatter)".format(figi))
    exec("ax{}.yaxis.set_major_formatter(lat_formatter)".format(figi))

    gd.plot(fmt=1,value=S.ctrl[nn,:]*100,clim=c2, cmap=cmap,cb=False,zorder=10)
    fpt = gd.dp < 10; tmp = S.ctrl[nn, :].copy(); tmp[fpt] = NaN
    CS=gd.plot(fmt=2,value=tmp*100,levels=[30,40,50],colors=['k','k','k'], linestyles=['solid','solid','solid'],linewidths=0.8,zorder=11)
    clabel(CS, inline=1,inline_spacing=-8, fontsize=fs)

    # title('Barotropic Adjustmnet')
    # xlim(xl);ylim(yl)
    # xticks(arange(xl[0], xl[1], 2));
    # yticks(arange(yl[0], yl[1], 2))
    if not figi==0: setp(gca(), xticklabels=[], yticklabels=[])
    # if figi==0: setp(gca(), xticklabels=[])

    exec("ax{}=subplot(4,3,{},projection=ccrs.PlateCarree())".format(figi,figi+4))
    exec("ax{}.set_extent(extent)".format(figi))
    exec("ax{}.coastlines()".format(figi))
    exec("ax{}.add_feature(states110, zorder=0, linewidth=1)".format(figi))
    exec("ax{}.set_xticks(linspace(extent[0], extent[1], num=3), crs=ccrs.PlateCarree())".format(figi))
    exec("ax{}.set_yticks(linspace(extent[2], extent[3], num=4), crs=ccrs.PlateCarree())".format(figi))
    exec("ax{}.xaxis.set_major_formatter(lon_formatter)".format(figi))
    exec("ax{}.yaxis.set_major_formatter(lat_formatter)".format(figi))
    gd.plot(fmt=1,value=S.atf[nn,:]*100,clim=c2, cmap=cmap,cb=False,zorder=10)
    fpt = gd.dp < 10; tmp = S.atf[nn, :].copy(); tmp[fpt] = NaN
    CS=gd.plot(fmt=2,value=tmp*100,levels=[20],colors=['k'], linestyles=['solid'],linewidths=0.8,zorder=11)
    clabel(CS, inline=1,inline_spacing=-8, fontsize=fs)
    # title('Wind')
    # xlim(xl);ylim(yl)
    # xticks(arange(xl[0], xl[1], 2));
    # yticks(arange(yl[0], yl[1], 2))
    if not figi==0: setp(gca(), xticklabels=[], yticklabels=[])
    if figi==0: setp(gca(), xticklabels=[],yticklabels=[])

    # setp(gca(), xticklabels=[])
    # colorbar()

    exec("ax{}=subplot(4,3,{},projection=ccrs.PlateCarree())".format(figi,figi+7))
    exec("ax{}.set_extent(extent)".format(figi))
    exec("ax{}.coastlines()".format(figi))
    exec("ax{}.add_feature(states110, zorder=0, linewidth=1)".format(figi))
    exec("ax{}.set_xticks(linspace(extent[0], extent[1], num=3), crs=ccrs.PlateCarree())".format(figi))
    exec("ax{}.set_yticks(linspace(extent[2], extent[3], num=4), crs=ccrs.PlateCarree())".format(figi))
    exec("ax{}.xaxis.set_major_formatter(lon_formatter)".format(figi))
    exec("ax{}.yaxis.set_major_formatter(lat_formatter)".format(figi))
    gd.plot(fmt=1,value=S.bcadj[nn,:]*100,clim=c2, cmap=cmap,cb=False,zorder=10)
    fpt=gd.dp<10; tmp=S.bcadj[nn,:].copy(); tmp[fpt]=NaN
    CS=gd.plot(fmt=2,value=tmp*100,levels=[20],colors=['k', 'k', 'k'], linestyles=['solid', 'solid', 'solid'],linewidths=0.8,zorder=11)
    clabel(CS, inline=1,inline_spacing=-8, fontsize=fs)
    # title('Baroclinic adjustment')
    # xlim(xl);ylim(yl)
    # xticks(arange(xl[0], xl[1], 2));
    # yticks(arange(yl[0], yl[1], 2))
    if not figi==0: setp(gca(), xticklabels=[], yticklabels=[])
    if figi==0: setp(gca(), xticklabels=[],yticklabels=[])

    # setp(gca(), xticklabels=[])
    # colorbar()

    exec("ax{}=subplot(4,3,{},projection=ccrs.PlateCarree())".format(figi,figi+10))
    exec("ax{}.set_extent(extent)".format(figi))
    exec("ax{}.coastlines()".format(figi))
    exec("ax{}.add_feature(states110, zorder=0, linewidth=1)".format(figi))
    exec("ax{}.set_xticks(linspace(extent[0], extent[1], num=3), crs=ccrs.PlateCarree())".format(figi))
    exec("ax{}.set_yticks(linspace(extent[2], extent[3], num=4), crs=ccrs.PlateCarree())".format(figi))
    exec("ax{}.xaxis.set_major_formatter(lon_formatter)".format(figi))
    exec("ax{}.yaxis.set_major_formatter(lat_formatter)".format(figi))
    tmp = S.btadj[nn, :]
    fpt = (tmp > 0.5); tmp[fpt]=NaN
    gd.plot(fmt=1,value=tmp*100,clim=c2, cmap=cmap,cb=False,zorder=10)
    # title('All forcing')
    # xlim(xl);ylim(yl)
    # xticks(arange(xl[0],xl[1],2)); yticks(arange(yl[0],yl[1],2))
    if not figi==0: setp(gca(), xticklabels=[], yticklabels=[])
    if figi==0: setp(gca(), xticklabels=[],yticklabels=[])

    # setp(gca(), xticklabels=[])
    # colorbar()
    # title('{}'.format(num2date(stime).strftime('%Y-%m-%d %H:%M:%S')), fontsize=11)
    # cbar = colorbar(ticks=arange(0,60+10,10), orientation='vertical', fraction=0.05, pad=0.1)

tight_layout()
savefig('{}'.format(sname) + '.png',dpi=900,bbox_inches='tight')
close()
