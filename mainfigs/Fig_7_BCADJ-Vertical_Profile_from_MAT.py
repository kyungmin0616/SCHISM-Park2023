from pylib import *
import cmocean
import mat73
import matplotlib.colors as colors


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

#import scipy.io as sio
###################################################################################
expname='CTRL-WO-AIR'
transect='NC'
varname='temperature' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

expname2='CTRL-WO-AIR'
transect2='NC'
varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

expname3='CTRL-WO-AIR'
transect3='NC'
varname3='horizontalVelX' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

N=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname,transect,varname))
N2=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname2,transect2,varname2))
N3=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname3,transect3,varname3))

expname='BTADJ'
transect='NC'
varname='temperature' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

expname2='BTADJ'
transect2='NC'
varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

expname3='BTADJ'
transect3='NC'
varname3='horizontalVelX' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

B=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname,transect,varname))
B2=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname2,transect2,varname2))
B3=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname3,transect3,varname3))

N['zv']=N['zv']-B['zv']; N2['zv']=N2['zv']-B2['zv']; N3['zv']=N3['zv']-B3['zv']



###################################################################################
expname='CTRL-WO-AIR'
transect='GA'
varname='temperature' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

expname2='CTRL-WO-AIR'
transect2='GA'
varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

expname3='CTRL-WO-AIR'
transect3='GA'
varname3='horizontalVelX' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

# rho=mat73.loadmat('./ECGOM_GA/outputs/mat/{}_{}_{}.mat'.format(expname,transect,'waterDensity'))
G=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname,transect,varname))
G2=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname2,transect2,varname2))
G3=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname3,transect3,varname3))

expname='BTADJ'
transect='GA'
varname='temperature' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

expname2='BTADJ'
transect2='GA'
varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

expname3='BTADJ'
transect3='GA'
varname3='horizontalVelX' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

# rho=mat73.loadmat('./ECGOM_GA/outputs/mat/{}_{}_{}.mat'.format(expname,transect,'waterDensity'))
B=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname,transect,varname))
B2=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname2,transect2,varname2))
B3=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname3,transect3,varname3))

G['zv']=G['zv']-B['zv']; G2['zv']=G2['zv']-B2['zv']; G3['zv']=G3['zv']-B3['zv'];


###################################################################################
expname='CTRL-WO-AIR'
transect='FL'
varname='temperature' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

expname2='CTRL-WO-AIR'
transect2='FL'
varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

expname3='CTRL-WO-AIR'
transect3='FL'
varname3='horizontalVelX' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

# rho=mat73.loadmat('./ECGOM_GA/outputs/mat/{}_{}_{}.mat'.format(expname,transect,'waterDensity'))
F=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname,transect,varname))
F2=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname2,transect2,varname2))
F3=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname3,transect3,varname3))

expname='BTADJ'
transect='FL'
varname='temperature' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

expname2='BTADJ'
transect2='FL'
varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

expname3='BTADJ'
transect3='FL'
varname3='horizontalVelX' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity
# varname2='horizontalVelY' #'horizontalVelX', turbulentKineticEner, diffusivity, verticalVelocity, viscosity

# rho=mat73.loadmat('./ECGOM_GA/outputs/mat/{}_{}_{}.mat'.format(expname,transect,'waterDensity'))
B=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname,transect,varname))
B2=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname2,transect2,varname2))
B3=mat73.loadmat('./ECGOM_GA/outputs/mat/PaperExp/{}_{}_{}.mat'.format(expname3,transect3,varname3))

F['zv']=F['zv']-B['zv']; F2['zv']=F2['zv']-B2['zv']; F3['zv']=F3['zv']-B3['zv'];
###################################################################################
# EKE
fpt=(datenum('2017-10-12')<=F2['tv'])*(F2['tv']<=datenum('2017-10-22'))

neke=[]; geke=[]; feke=[]
for ekei in arange(len(F2['tv'])):
    neke.append(((N2['zv'][:,:,ekei]-N2['zv'][:,:,fpt].mean(axis=2))**2+(N3['zv'][:,:,ekei]-N3['zv'][:,:,fpt].mean(axis=2))**2)/2)
    geke.append(((G2['zv'][:,:,ekei]-G2['zv'][:,:,fpt].mean(axis=2))**2+(G3['zv'][:,:,ekei]-G3['zv'][:,:,fpt].mean(axis=2))**2)/2)
    feke.append(((F2['zv'][:,:,ekei]-F2['zv'][:,:,fpt].mean(axis=2))**2+(F3['zv'][:,:,ekei]-F3['zv'][:,:,fpt].mean(axis=2))**2)/2)

neke=transpose(neke,axes=[1,2,0]); neke=array(neke)
geke=transpose(geke,axes=[1,2,0]); geke=array(geke)
feke=transpose(feke,axes=[1,2,0]); feke=array(feke)

################################## Plot EKE
fpt=(datenum('2017-10-12')<=F2['tv'])*(F2['tv']<=datenum('2017-10-22'))
cmap='Reds'
# levels= linspace(7, 29,10)
levels= arange(0, 1,2)
levels2= linspace(0, 0.01,40)
ctc='k'
ym=[0,1000]
# norm = colors.TwoSlopeNorm(vmin=levels[0], vcenter=0, vmax=levels[-1])
j=0
figure(1, figsize=[7.3, 4.0])
clf()
for sindt in arange(1,10):
    # if sindt==7: continue
    j=j+1
    sindt2=sindt*24
    sindt=arange(24*sindt,24*(sindt+1))
    time=N['tv'][sindt2]
    print(num2date(time))

    subplot(3,9,j)
    levels = arange(0.5,0.9+0.1,0.1)
    levels2 = linspace(0, 0.1, 40)
    levels2 = linspace(0, 1, 40)
    cbarlabels = arange(levels2.min(), levels2.max() + 0.1, 0.1)

    ctrlx=-130
    xm = [130+ctrlx, 280+ctrlx]
    # title(num2date(time).strftime('%m-%d'))
    CS=contour(N['xv'][:, :, sindt].mean(axis=2)+ctrlx, -N['yv'][:, :, sindt].mean(axis=2), N2['zv'][:, :, sindt].mean(axis=2), levels=levels, colors=ctc, linestyles='-', linewidths=1)
    # clabel(CS,fmt='%2.2f', colors='k', fontsize=8)
    map2=contourf(N2['xv'][:, :, sindt].mean(axis=2)+ctrlx, -N2['yv'][:, :, sindt].mean(axis=2), neke[:,:,sindt].mean(axis=2)/0.1, levels=levels2,cmap=cmap,extend='both')
    gca().invert_yaxis()
    setp(gca(), xticks=arange(xm[0],xm[-1]+50,50),xticklabels=arange(xm[0],xm[-1]+50,50), xlim=xm, ylim=ym)
    gca().invert_yaxis()
    xticks(rotation=45)

    if j > 1: setp(gca(), xticklabels=[], yticklabels=[])
    setp(gca(), xticklabels=[])
    if j == 1: ylabel('Depth (m)')

    # cbar = colorbar(map2, ticks=cbarlabels, orientation='vertical')
    # cbar.ax.tick_params(labelsize=10)
    # setp(gca(), xticklabels=[], yticklabels=[])

    subplot(3,9,j+9)
    levels = arange(0.2,0.6+0.1,0.1)
    levels2 = linspace(0, 0.02, 40)
    levels2 = linspace(0, 1, 40)
    cbarlabels = arange(levels2.min(), levels2.max() + 0.1, 0.1)

    ctrlx=-120
    xm = [120+ctrlx, 270+ctrlx]
    # title(num2date(time).strftime('%m-%d'))
    CS=contour(G['xv'][:, :, sindt].mean(axis=2)+ctrlx, -G['yv'][:, :, sindt].mean(axis=2), G2['zv'][:, :, sindt].mean(axis=2), levels=levels, colors=ctc, linestyles='-', linewidths=1)
    # clabel(CS,fmt='%2.2f', colors='k', fontsize=8)
    contourf(G2['xv'][:, :, sindt].mean(axis=2)+ctrlx, -G2['yv'][:, :, sindt].mean(axis=2), geke[:,:,sindt].mean(axis=2)/0.02, levels=levels2, cmap=cmap,extend='both') #-S['zv'][:,:,tmp].mean(axis=2)
    gca().invert_yaxis()
    setp(gca(), xticks=arange(xm[0],xm[-1]+50,50),xticklabels=arange(xm[0],xm[-1]+50,50), xlim=xm, ylim=ym)
    gca().invert_yaxis()
    xticks(rotation=45)

    if j > 1: setp(gca(), xticklabels=[], yticklabels=[])
    setp(gca(), xticklabels=[])
    if j == 1: ylabel('Depth (m)')

    # cbar = colorbar(map2, ticks=cbarlabels, orientation='vertical')
    # cbar.ax.tick_params(labelsize=10)
    # setp(gca(), xticklabels=[], yticklabels=[])

    subplot(3,9,j+18)
    levels = arange(0.6,1+0.1,0.1)
    levels2 = linspace(0, 0.1, 40)
    levels2 = linspace(0, 1, 40)
    cbarlabels = arange(levels2.min(), levels2.max() + 0.1, 0.1)
    ctrlx=-60
    xm = [60+ctrlx, 210+ctrlx]
    # title(num2date(time).strftime('%m-%d'))
    CS=contour(F['xv'][:, :, sindt].mean(axis=2)+ctrlx, -F['yv'][:, :, sindt].mean(axis=2), F2['zv'][:, :, sindt].mean(axis=2), levels=levels, colors=ctc, linestyles='-', linewidths=1)
    # clabel(CS,fmt='%2.2f', colors='k', fontsize=8)
    contourf(F2['xv'][:, :, sindt].mean(axis=2)+ctrlx, -F2['yv'][:, :, sindt].mean(axis=2), feke[:,:,sindt].mean(axis=2)/0.1, levels=levels2, cmap=cmap,extend='both') #-S['zv'][:,:,tmp].mean(axis=2)
    setp(gca(), xticks=arange(xm[0],xm[-1]+50,50),xticklabels=arange(xm[0],xm[-1]+50,50), xlim=xm, ylim=ym)
    gca().invert_yaxis()
    xticks(rotation=45)
    if j > 1: setp(gca(), xticklabels=[], yticklabels=[])
    if j == 1: xlabel('Distance (km)'); ylabel('Depth (m)')

    # setp(gca(), xticklabels=[])

    # cbar = colorbar(map2, ticks=cbarlabels, orientation='vertical')
    # cbar.ax.tick_params(labelsize=10)
    # setp(gca(), xticklabels=[], yticklabels=[])

    # rcParams['axes.facecolor'] = 'grey'
    # savefig('{}/{}'.format('/Users/kpark350/Downloads/2nd_Paper/ver_tran', j) + '.png')
    # close()
tight_layout()
savefig('/Users/kpark350/Downloads/2nd_Paper/Fig_7.png',dpi=900,bbox_inches='tight')
close()


