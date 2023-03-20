#!/usr/bin/env python3
from pylib import *
import ttide.t_tide as t_tide


from matplotlib.backends.backend_pdf import PdfPages
close("all")
#------------------------------------------------------------------------------
#input
# ------------------------------------------------------------------------------
dir_obs='./obs/USGS/USGS_021989773_WL_2016.csv'
dir_obs2='./obs/USGS/USGS_021989773_RD_2016.csv'
dir_obs3='./obs/USGS/USGS_021989773_PREC_2016.csv'

df = pd.read_csv(dir_obs)
df2 = pd.read_csv(dir_obs2)
df3 = pd.read_csv(dir_obs3)

# print(df)
st=datenum('2016-10-1'); et=datenum('2016-10-15')
otii=array(datenum(df.values[:,2]))
oyii=array(df.values[:,4],dtype='float')*0.3048 # ft to m

otii2=array(datenum(df2.values[:,2]))
oyii2=array(df2.values[:,4],dtype='float')*0.0283168 # cubic ft to cubie m

otii3=array(datenum(df3.values[:,2]))
oyii3=array(df3.values[:,4],dtype='float')*0.0254 # inch to m

#HA
# ts=find_continuous_sections(otii,1)
# fp=(otii>=ts.sections[0,0])*(otii<=ts.sections[0,1])
# HAO = t_tide(oyii[fp],1/4,st,array(32.16555556),synth=0)
foyi = lpfilt(oyii, 1 / 240 * (15/6), 13 / 24)
foyi2 = lpfilt(oyii2, 1 / 240 * (15/6), 13 / 24)


st=datenum('2016-10-1'); et=datenum('2016-10-20')

clf()
plot(otii,foyi,'k') # water level
plot(otii3,oyii3*100,'b') # water level

# setp(gca(),ylim=[-0.5*100, 1.5*100],yticks=linspace(-0.5*100,1.5*100,9))
gca().xaxis.grid('on')
gca().yaxis.grid('on')

twinx()
plot(otii2,foyi2,'g')
# setp(gca(),ylim=[-50, 180],yticks=linspace(-50,180,9))
gca().xaxis.grid('on')
gca().yaxis.grid('on')
xts, xls = get_xtick(fmt=2, xts=linspace(st,et,10), str='%m-%d')
setp(gca(), xticks=xts, xticklabels=xls, xlim=[st, et])
gcf().tight_layout()



SMALL_SIZE = 17
MEDIUM_SIZE = 17
BIGGER_SIZE = 17
rc('font', size=SMALL_SIZE)  # controls default text sizes
rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title