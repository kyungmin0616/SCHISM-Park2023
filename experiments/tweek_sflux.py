from pylib import *

ctlt=34
dir_sflux='../sflux_air1'
fnames=array([i for i in os.listdir(dir_sflux) if i.endswith('.nc') and i.__contains__('air')])
mti=array([array(i.replace('.','_').split('_')[3]).astype('int')for i in fnames])
fpt=(mti>ctlt); fnames=fnames[fpt]; mti=mti[fpt]
sind=argsort(mti); mti=mti[sind]; fnames=fnames[sind]

S=ReadNC('{}/sflux_air_1.00{}.nc'.format(dir_sflux,ctlt))

for i in arange(0,len(S.time.val)):
    S.uwind.val.data[i,:,:]=S.uwind.val.data[i,:,:]*(-(1/(len(S.time.val)-1))*i+1)
    S.vwind.val.data[i,:,:]=S.uwind.val.data[i,:,:]*(-(1/(len(S.time.val)-1))*i+1)
    S.prmsl.val.data[i,:,:]=S.prmsl.val.data[i,:,:]+(101325-S.prmsl.val.data[i,:,:])*(i/(len(S.time.val)-1))
WriteNC('./sflux_air_1.00{}.nc'.format(ctlt),S)

#S=ReadNC('{}/sflux_air_2.00{}.nc'.format(dir_sflux,ctlt))

#for i in arange(0,len(S.time.val)):
#    S.uwind.val.data[i,:,:]=S.uwind.val.data[i,:,:]*(-(1/(len(S.time.val)-1))*i+1)
#    S.uwind.attrs.remove('_FillValue')
#    S.vwind.val.data[i,:,:]=S.uwind.val.data[i,:,:]*(-(1/(len(S.time.val)-1))*i+1)
#    S.prmsl.val.data[i,:,:]=S.prmsl.val.data[i,:,:]+(101325-S.prmsl.val.data[i,:,:])*(i/(len(S.time.val)-1))
#    for j in S.vars:
#        S.{}.attrs.remove('_FillValue')


#WriteNC('./sflux_air_2.00{}.nc'.format(ctlt),S)

#for fname in fnames:
#    S=ReadNC('{}/{}'.format(dir_sflux,fname))
#    for i in arange(0,len(S.time.val)):
#        S.uwind.val.data[i,:,:]=S.uwind.val.data[i,:,:]*0
#        S.vwind.val.data[i,:,:]=S.uwind.val.data[i,:,:]*0
#        S.prmsl.val.data[i,:,:]=S.prmsl.val.data[i,:,:]*0+101325
#    WriteNC('./{}'.format(fname),S)


