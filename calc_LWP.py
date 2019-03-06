""" Script to convert an src=ensmean$INFO.nc file 
    to a tar=LWP$INFO.nc file,
    where $tar denotes the LWP variables
"""

import numpy as np
import os
from netCDF4 import Dataset
from circmeantonc_projectB2 import Rdry
from analy_utils import HHL_creator, nhalo

BASEnowind = '/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/ALLVAR_3D/'
BASEwind05 = '/hymet/adeli/project_B2/512x512_7Kkmwind05_new_1km/postprocessing/composites/ALLVAR_3D/'
BASEwind1 = '/hymet/adeli/project_B2/512x512_7Kkmwind_1km/postprocessing/composites/ALLVAR_3D/'
BASEwind2 = '/hymet/adeli/project_B2/512x512_7Kkmwind2_new_1km/postprocessing/composites/ALLVAR_3D/'
BASEwind4 = '/hymet/adeli/project_B2/512x512_7Kkmwind4_new_1km/postprocessing/composites/ALLVAR_3D/'
BASEA = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/ALLVAR_3D/'

BASE = BASEnowind

BASEU5 = '/net/o3/hymet/adeli/project_B2/512x512_7KkmU5_1km/postprocessing/composites/ALLVAR_3D/'
BASEU10= '/net/o3/hymet/adeli/project_B2/512x512_7KkmU10_1km/postprocessing/composites/ALLVAR_3D/'



ENSN = 'ensmean_day_3d_d0d5.nc'
LWPN = 'LWP_'+ENSN


qtvars = ['Q'+x for x in 'VCIRSG'] #VCIRSG?

offset=0
nhalo=nhalo+offset

#hexps='h250a15 h500a14 h750a14 h1000a14 h1500a14 h500a5 h500a7 h500a10 h50020 h500a25 h500a30'



arange = [5,7,10,14,20,25,30]
hrange = [250,500,750,1000,1500,3000]

na = len(arange)
nh = len(hrange)

expnames=np.zeros([nh,na],dtype=object)
for i in range(nh):
    for j in range(na):
        expnames[i,j] = 'h'+str(hrange[i])+'a'+str(arange[j])


def copyncdims_dimvars(srcnc,tarnc):
    """ copies dimensions from srcnc to the tarnc file."""
    for dim in srcnc.dimensions.values():
        name=dim.name
        size=dim.size
        if size==518:
            size-=2*nhalo
        tarnc.createDimension(name,size=size)

    for var in ['lon','lat']:
        dims = srcnc.variables[var].dimensions
        tarnc.createVariable(var,float,dimensions=dims)
        tarnc.variables[var][:] = srcnc.variables[var][nhalo:-nhalo,nhalo:-nhalo]



if __name__=='__main__':
    
    iparallel = False
    if iparallel:    
        import argparse
        parser = argparse.ArgumentParser(description='Process some integers.')
        parser.add_argument('expnames', metavar='EXPN', type=str, nargs='+',default='h500a20',
                            help='experiment name')
        args = parser.parse_args()
        expns = args.expnames
        #print 'working on a = ' + str(myarange[0])
    else:
        expns = ['h4000a14']
        #myarange = range(na)
    
    expns = ['h1000a14','h2000a14','h3000a14','h4000a14']
    print expns

    OUTPUTPATHS = [BASEU5, BASEU10] 
            
    for BASE in OUTPUTPATHS:
        for expn in expns:
            ah=expn.replace('a','h').split('h')[1:]
            hm,am=int(ah[0]),int(ah[1])
           

            H=HHL_creator(Hm=hm,nx=512-2*offset,ny=512-2*offset,nz=51,a=am)
            H[:] = H[::-1]
            Hhalf = H[1:,:]-H[0:-1,:] 
            Hhalf = Hhalf[::-1,:]


            basepath = BASE+expn+'/60_homo/'
            iensmean = True
            if iensmean:
                srcn = ENSN
            else:
                seedn=filter(lambda x: x.startswith('seed'), os.listdir(basepath)) 
                
                if seedn==[]: continue
                    
                seedn=seedn[0]
                print expn,seedn
                
                srcn = seedn #ENSN
            tarn = 'LWP_'+srcn           
            
            srcfln = basepath+srcn
            tarfln = basepath+tarn

            print tarfln
            
            # Prepare netcdf files
            srcncfl = Dataset(srcfln,'r')
            tarncfl = Dataset(tarfln,'w')
            
            # Error check
            # srcncfl.close()
            # tarncfl.close()
            # continue

            copyncdims_dimvars(srcncfl,tarncfl)

            # Create the new variables

            for qtvar in qtvars:
                tarncfl.createVariable('LWP_'+qtvar,float,dimensions=('time','rlat','rlon'))


            # Copy surface heat fluxes
            tarncfl.createVariable('Esurf',float,dimensions=('time','rlat','rlon'))
            tarncfl.createVariable('Hsurf',float,dimensions=('time','rlat','rlon'))
            tarncfl.createVariable('RHOsurf',float,dimensions=('time','rlat','rlon'))
            
            
            Esurf = srcncfl.variables['EFLUX'][:,-1,nhalo:-nhalo,nhalo:-nhalo]
            Hsurf = srcncfl.variables['HFLUX'][:,-1,nhalo:-nhalo,nhalo:-nhalo]
            
            tarncfl.variables['Esurf'][:] = Esurf
            tarncfl.variables['Hsurf'][:] = Hsurf
            

            # Copy the dimensions from one file to the other
            T = srcncfl.variables['T'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            P = srcncfl.variables['P'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            RHO = P/T/Rdry
            
            # save surface density
            tarncfl.variables['RHOsurf'][:] = RHO[:,-1,:,:]
            
            


            for qtvar in qtvars:
                print qtvar
                for ti in range(48):
                    print ti
                    var = srcncfl.variables[qtvar][:,:,nhalo:-nhalo,nhalo:-nhalo]
                    LWPv = np.sum(RHO[ti,:]*var[ti,:]*Hhalf,axis=0)
                    tarncfl.variables['LWP_'+qtvar][ti,:]=LWPv


            srcncfl.close()
            tarncfl.close()
