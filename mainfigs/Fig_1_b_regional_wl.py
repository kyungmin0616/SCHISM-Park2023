#!/usr/bin/env python3
from pylib import *
from matplotlib.backends.backend_pdf import PdfPages
close("all")

#------------------------------------------------------------------------------
#input
# ------------------------------------------------------------------------------
StartT_model=datenum(2016,9,8) #plot start time, model start time
runs=['RUN05b-wo-tide-elev-xyz.npz']
tags=['SCHISM']
bpfile='./ECGOM_GA/station.in'
stt=datenum('2016-9-22'); edt=datenum('2016-10-22')
fst=datenum('2016-9-18'); fet=datenum('2016-10-2')
regions=['SAB']
#figure save names; not save if it is None
#sname=None
sname='/Users/kpark350/Dropbox (GaTech)/Old_Mac/Ga_tech/Journal/2nd_GulfStream/Figure/Fig_1_b'

#linestypes for different runs
colors='kgbcm'; lstyle=['-','-','-','-','-']; markers=['None','None','None','None','None']
# colors='gbc'; lstyle=['-','-','-']; markers=['None','.','None']

#msl_to_navd shift
# msl_shift_station=array([8760721]); msl_shift=[-0.5]
msl_shift_station=None; msl_shift=None

#shift all data for each model run
all_shift=[0,0,0,0,0]

#Axis limits
ym=None; xm=[1,50]
ex = linspace(xm[0],xm[1],10)

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

##################################  read obs data
C = loadz('./obs/NOAA/Matthew/WL/noaa_elev_navd.npz')
C1 = loadz('./obs/NOAA/Matthew/WL/noaa_elev_msl.npz')
C2=loadz('./obs/NOAA/Matthew/TIDE/noaa_elev_navd.npz')
C3=loadz('./obs/NOAA/Matthew/TIDE/noaa_elev_msl.npz')


################################## stations to plots
noaa_stations_groups = {

'FT': array([11]),

'SAB': array([7,8,10,11,12,13,15,16]), # 16 ,17 need to be included in KMP

'MAB': array([19,21,23,26,27,28,30]),

'GOME': array([33,36,37,38]),

'GOMX': array([44, 49, 57, 59, 64, 67, 74, 77, 78, 2, 3]), #44~64 west; 67~3 east

'BPT': array([40,151,152,153,155,161] ),

'ECOAST': array([10,11,12,13,15,19,21,23,26,27,28,30,33,36,37,38] ),

'ECOAST-GOMX': array([7,8,10,11,12,13,15,19,21,23,26,27,28,30,33,36,37,38,44, 49, 57, 59, 64, 67, 74, 77, 78, 2, 3] ),

'All': array([2, 7,8,10,11,12,13,15,19,21,23,26,27,28,30,33,36,37,38,3, 78, 74, 64, 59, 49, 44, 40, 151,152,153,155,161,57,77,67])

}

for region in regions:
    if region == 'None': stations=None
    else: stations=noaa_stations_groups[region]

    #------------------------------------------------------------------------------
    #read station data
    #------------------------------------------------------------------------------
    #read station.bp
    fid=open('./ECGOM_GA/stanames.txt'); staname=dict(zip(*array([i.strip().split(',') for i in fid.readlines()]).T)); fid.close()
    bp=read_schism_bpfile(bpfile);  bp.staname=array([staname[i] for i in bp.station])

    #for subset of stations
    if stations is not None:
        staid=array(stations).astype('int')-1
        bp.nsta=len(stations); bp.x=bp.x[staid]; bp.y=bp.y[staid]; bp.z=bp.z[staid]; bp.station=bp.station[staid]; bp.staname=bp.staname[staid]
    else:
        staid=arange(bp.nsta).astype('int')
        stations=arange(bp.nsta)+1

    ##################################  read model results
    Model = []
    for m, run in enumerate(runs):
        Model.append(loadz('./ECGOM_GA/npz/{}'.format(run)))

    #-----------------------plot---------------------------------------------------
    figure(2, figsize=[3.5, 3.7])
    # clf()
    offset = 0
    iflag=0
    MAE=[]; Corr=[]; ME=[]; RMSE=[];
    clf()
    for i in arange(bp.nsta):
        iflag=iflag+1
        station=int(bp.station[i])
        #find target time and station of obs
        lobs='None'
        fp=(C.station==station)*(C.time>stt-2)*(C.time<edt+2); oti=C.time[fp]; oyi=C.elev[fp]; lobs='navd'
        fp = (C2.station == station) * (C2.time > stt-2) * (C2.time < edt+2); oti_tide = C2.time[fp]; oyi_tide = C2.elev[fp];

        if sum(fp)==0:
            fp=(C1.station==station)*(C1.time>stt-2)*(C1.time<edt+2); oti=C1.time[fp]; oyi=C1.elev[fp]; lobs='msl'
            fp = (C3.station == station) * (C3.time > stt-2) * (C3.time < edt+2); oti_tide = C3.time[fp]; oyi_tide = C3.elev[fp];

        #add nan data between oti
        if len(oti)>100:
           ts=find_continuous_sections(oti,1.0); eoti=array([i[-1]+1/24 for i in ts.sections]); eoyi=ones(len(eoti))*nan
           oti=r_[oti,eoti]; oyi=r_[oyi,eoyi]; sind=argsort(oti); oti=oti[sind]; oyi=oyi[sind]
           ts=find_continuous_sections(oti_tide,1.0); eoti=array([i[-1]+1/24 for i in ts.sections]); eoyi=ones(len(eoti))*nan
           oti_tide=r_[oti_tide,eoti]; oyi_tide=r_[oyi_tide,eoyi]; sind=argsort(oti_tide); oti_tide=oti_tide[sind]; oyi_tide=oyi_tide[sind]
        if len(oyi)==0: continue
        if len(oyi)==len(oyi_tide):
            fpn = ~isnan(oyi);
            oti = oti[fpn];
            oyi = oyi[fpn]; oyi_tide=oyi_tide[fpn]
            oyi=oyi-oyi_tide
            foyi = lpfilt(oyi, 1 / 240, 13/ 24)
            # foyi = smooth(oyi, 100)
            # foyi=oyi
            ffpt=(oti>fst)*(oti<fet);
            fptc=(oti>datenum('2016-09-18'))*(oti<datenum('2016-10-4'))
            plot(oti[fptc], foyi[fptc]+offset-nanmean(foyi[ffpt]), 'k',lw=3)
            fptc=(oti>datenum('2016-10-4'))*(oti<datenum('2016-10-10'))
            plot(oti[fptc], foyi[fptc]+offset-nanmean(foyi[ffpt]), color='orange',lw=3)
            fptc=(oti>datenum('2016-10-10 00:00:00'))*(oti<datenum('2016-10-22'))
            plot(oti[fptc], foyi[fptc]+offset-nanmean(foyi[ffpt]), color='orangered',lw=3)
        else:
            foyi = lpfilt(oyi, 1 / 240, 13/ 24)
            plot(oti, foyi+offset, 'r',lw=3)

        # plots and stats for models
        for mn, run in enumerate(runs):
            mti = Model[mn].time
            fp = (Model[mn].bp.station == bp.station[i])
            myi = squeeze(Model[mn].elev[:, fp])
            if any(abs(myi) > 100): continue
            fpn = ~isnan(foyi); oti = oti[fpn]; oyi = foyi[fpn]
            myi = lpfilt(myi, 1 / 24, 13 / 24)
            ffpt = (mti > fst) * (mti < fet);
            plot(mti[::18], myi[::18] + offset - mean(myi[ffpt]), linestyle='None', color=colors[0], marker='o',markerfacecolor='None', ms=4,linewidth=0.5)
            # plot(mti, myi + offset - mean(myi[ffpt]), linestyle='None', color=colors[0], marker='o',markerfacecolor='None', ms=4,linewidth=0.5)

        if ym is not None: setp(gca(), ylim=ym)
        xts, xls = get_xtick(fmt=2, xts=arange(stt,edt+5,5), str='%m-%d')
        setp(gca(), xticks=xts, xticklabels=xls,yticks=arange(0,16,1), xlim=[stt, edt],ylim=[-0.5,15])
        gca().xaxis.grid('on')
        gca().yaxis.grid('on')
        offset = offset +2
ylabel('Water level (m)')
xlabel('Date (2016-)')
tight_layout()
# savefig('{}'.format(sname) + '.png',dpi=900,bbox_inches='tight')
# close()
