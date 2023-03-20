from pylib import *
import cmocean
from pyproj import Geod
geod = Geod(ellps='WGS84')

##########################################################################################
sname='/Users/park075/Dropbox (GaTech)/Old_Mac/Ga_tech/Journal/2nd_GulfStream/Figure/Fig_4_a'
S=loadz('./ECGOM_GA/npz/PaperExp/Comb-exp-elev-33_45.npz')
bp=read_shapefile_data('./ECGOM_GA/SHP/Paper_SAB_COAST.shp')
gd=read_schism_hgrid('./ECGOM_GA/grd/hgrid.gr3')
st=datenum(2016,10,11)
se=datenum(2016,10,22)
ym=[0,1]
##########################################################################################
px,py=array(bp.xy.T,dtype=float)
px = px[~isnan(px)]; py = py[~isnan(py)]
px,py=proj_pts(px,py,'epsg:26918','epsg:4326')
sindp=near_pts(c_[px,py],c_[gd.x,gd.y])
sindp = [i for n, i in enumerate(sindp) if i not in sindp[:n]] # delete dupulicate remain the oder of elements

fpt=(S.time>=st)*(S.time<=se); time=S.time[fpt]
ctrl=S.ctrl[fpt,:]; atf=S.atf[fpt,:]; btadj=S.btadj[fpt,:]; bcadj=S.bcadj[fpt,:]
ctrl=ctrl[:,sindp]; atf=atf[:,sindp]; btadj=btadj[:,sindp]; bcadj=bcadj[:,sindp]
ctrl_mean=ctrl.mean(axis=1); atf_mean=atf.mean(axis=1); btadj_mean=btadj.mean(axis=1); bcadj_mean=bcadj.mean(axis=1)
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
# plotting
figure(1, figsize=[3.5, 5.5])
clf()
plot(time,ctrl_mean*100,lw=2,color='k')
plot(time,atf_mean*100,lw=2,color='b')
plot(time,bcadj_mean*100,lw=2,color='r')
plot(time,btadj_mean*100,lw=2,color='g')
plot(time,(btadj_mean+bcadj_mean)*100,lw=2,color='m')

# legend(['CTRL','ATF','BC-ADJ','BT-ADJ'],loc='upper right')
xts, xls = get_xtick(fmt=2, xts=[*arange(st,se+1,1)], str='%d')
setp(gca(), xticks=xts, xticklabels=xls, xlim=[st, se], ylim=[-5,50])
gca().xaxis.grid('on');gca().yaxis.grid('on')
xlabel('Date (2016-10-)'); ylabel('Water level (cm)')

tight_layout()
# savefig('{}'.format(sname) + '.png',dpi=900,bbox_inches='tight')
# close()

#
# ########## Percentage
# # plotting
# SMALL_SIZE = 18
# MEDIUM_SIZE = 20
# BIGGER_SIZE = 20
# rc('font', size=SMALL_SIZE)  # controls default text sizes
# rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
# rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
# rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
# rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
# rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
# rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
figure(2, figsize=[19, 6])
clf()
plot(time,ctrl_mean/ctrl_mean*100,lw=2,color='k')
plot(time,atf_mean/ctrl_mean*100,lw=2,color='b')
plot(time,bcadj_mean/ctrl_mean*100,lw=2,color='r')
plot(time,btadj_mean/ctrl_mean*100,lw=2,color='g')
plot(time,(btadj_mean+bcadj_mean)/ctrl_mean*100,lw=2,color='m')

legend(['CTRL','ATF','BC-ADJ','BT-ADJ'],loc='upper right')
xts, xls = get_xtick(fmt=2, xts=[*arange(st,se+1,1)], str='%m-%d')
setp(gca(), xticks=xts, xticklabels=xls, xlim=[st, se])
xlabel('Date (2016)'); ylabel('Contribution (%)')
gca().xaxis.grid('on');gca().yaxis.grid('on')