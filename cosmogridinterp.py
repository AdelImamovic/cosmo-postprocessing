# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 09:37:39 2016

@author: adeli

""" 
""" Given an nc-file the circular mean around the center is calculated"""


import sys
sys.path.append('/home/adeli/scripts/python/')
from postprocessing_utils import circsym_mean_vec, circsym_mean_2D, circsym_mean_scal
from analy_utils import HHL_creator, nhalo, getha
import numpy as np 
import os

from netCDF4 import Dataset
import matplotlib.pylab as plt
from scipy.interpolate import interp1d


BASEw = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_1km/postprocessing/composites/'
BASEnow = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
BASE = BASEw+'/ALLVAR_3D/'


BASEw05 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_new_1km/postprocessing/composites/ALLVAR_3D/'
BASEw05_old = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_1km/postprocessing/composites/ALLVAR_3D/'
BASEnow = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/ALLVAR_3D/'
BASEw1 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind_1km/postprocessing/composites/ALLVAR_3D/'
BASEw2 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind2_new_1km/postprocessing/composites/ALLVAR_3D/'
BASEw4 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind4_new_1km/postprocessing/composites/ALLVAR_3D/'

BASEU5 = '/net/o3/hymet/adeli/project_B2/512x512_7KkmU5_1km/postprocessing/composites/ALLVAR_3D/'
BASEU10= '/net/o3/hymet/adeli/project_B2/512x512_7KkmU10_1km/postprocessing/composites/ALLVAR_3D/'


Rdry = 287.058 # gas constant of dry air
invRdry = 1. / Rdry

Lv = 2264705. #J/kg
c_p = 1005.7 #J/kg
g_g = 9.81 # m/s2


#
smsall = ['60_30','60_40','60_50','60_70','60_80','60_90',
          '30_homo','40_homo','50_homo','60_homo','70_homo','80_homo','90_homo'] 
smsensi = ['40_10','40_20','40_30','40_50','40_60','40_70']

smsallold = ['60_30','60_40','60_50','60_70','60_80','60_90',
          '40_homo','60_homo','80_homo']
#
orosall = ['flat','125m','250m','500m']
oroheight = {'flat' : 0.0, '125m':125.0, '250m' : 250.0, '500m' : 500.0}


#HHLfile = Dataset('/home/adeli/scripts/python/HHL_fields/HHL_flat.nc')
#HHLhalf = HHLfile.variables['HHL'][0,:] - 500.
#HHLfull = (HHLhalf[0:-1,:] + HHLhalf[1:,:])*0.5
#HHLfile.close()
#
#gridtype='full'
#modgrid = HHLfull if gridtype == 'full' else HHLhalf
#
#intgrid = np.arange(-1000,20000,100)
#truefun = np.random.normal(0,1,size=len(intgrid))
#
#x=modgrid[::-1,0,0]
#y=intgrid
#data=np.random.normal(0,1,(len(x)))
#fun_i = interp1d(x,data,'linear',bounds_error=False,fill_value=0.)
#
#f,ax=plt.subplots(1,1)
#ax.plot(x,data,'o')
#ax.plot(y,fun_i(y),'-')
#assert 0



# Weired behavior
# from calc_LWP import hrange,arange,expnames
arange = [5,7,10,14,20,25,30]
hrange = [250,500,750,1000,2000,4000]

na = len(arange)
nh = len(hrange)

expnames=np.zeros([nh,na],dtype=object)
for i in range(nh):
    for j in range(na):
        expnames[i,j] = 'h'+str(hrange[i])+'a'+str(arange[j])
        
        
def paramextractor(oroname):
    """Extracts oroparameters from an experimental name."""
    
    #oroname=oroname[1:] #clip of 'h'
    
    expnames = {'h500ax10ay20':(500,10,20,'gauss'),
                 'h500ax10ay40':(500,10,40,'gauss'),
                 'h500ax10ay60':(500,10,60,'gauss'),
                 'h500ax20ay10':(500,20,10,'gauss'),
                 'h500ax40ay10':(500,40,10,'gauss'),
                 'h500ax60ay10':(500,60,10,'gauss'),
                 'h1000ax10ay20':(1000,10,20,'gauss'),
                 'h1000ax10ay40':(1000,10,40,'gauss'),
                 'h1000ax10ay60':(1000,10,60,'gauss'),
                 'h1000ax20ay10':(1000,20,10,'gauss'),
                 'h1000ax40ay10':(1000,40,10,'gauss'),
                 'h1000ax60ay10':(1000,60,10,'gauss')}
    
    
    return expnames[oroname]
    
def paramextractor_2(oroname):
    oroname = oroname[1:]
    spl=oroname.split('ax')
    h = spl[0]
    a,ay = spl.split('ay')
    if oroname.endswith('bell'):
        surftopo='bell'
    if oroname.endswith('cos2'):
        surftopo='cos2'
    else:
        surftopo='gauss'
    return h,a,ay,surftopo
    
    

def modgrid2interpgrid(oro,field,intgrid,gridtype='full'):  
    """TODO; bound_error false/true"""
    #filename ='/home/adeli/scripts/python/HHL_fields/512x512/HHL_'+oro+'.nc'
    #HHLfile = Dataset(filename)
    #HHLfile.close() 
     
    
    
    if oro.endswith('bell') or oro.endswith('cos2'):
        oro = oro[:-5]
        surftopo=oro[-4:]
        print oro  
        i,j = np.where(expnames==oro)
        h,a = hrange[i],arange[j]
        ay = a
    if oro.find('ax') > 0:
        h,a,ay,surftopo = paramextractor(oro) 
        print h,a,ay,surftopo
    else:
        surftopo='gauss'
        print oro
        h,a = getha(oro)
        ay=a
        print h,a
    
    
    nz,nx,ny = field.shape
    print field.shape
    HHLhalf = HHL_creator(Hm=h,nx=nx,ny=ny,nz=51,a=a,ay=ay,surftopo=surftopo) 
        
        
        
    #HHLhalf = HHLfile.variables['HHL'][0,:] - 500.
    HHLfull = (HHLhalf[0:-1,:] + HHLhalf[1:,:])*0.5 # determine height of full levs
       

    modgrid = HHLfull if gridtype == 'full' else HHLhalf

    field_interp = np.zeros(intgrid.shape)

    for i in range(nx):
        for j in range(ny):
            f = interp1d(modgrid[:,i,j], field[:,i,j], kind='linear',
                         bounds_error=False,fill_value=np.nan)
            field_interp[:,i,j] = f(intgrid[:,i,j])           

    return field_interp

from scipy.interpolate import griddata    
def modgrid2interpgrid_3d(oro,field,intgrid,gridtype='full'):  
    """ Overcomes shortcomings of modgrid2interpgrid (1d interpolation)."""
    modgrid = HHLfull if gridtype == 'full' else HHLhalf
    
    nz,nx,ny = field.shape
    field_interp = np.zeros(intgrid.shape)

    # fails as intgrid needs to be a meshgrid dim^3
    try:
        field_interp = griddata(modgrid,field,intgrid,method='linear')
    except:
        assert 0
        
    
    return field_interp



if __name__ == '__main__':
   
    OUTPUTPATHS = [BASEnow] 
    interpmethod = modgrid2interpgrid
    # select the files to iterate over, cross product is built oros 'x' sms
    # oros and sms must provide an iterable
    oros = filter(lambda x: x.startswith('h2000a'), os.listdir(BASE))
    oros = ['h250a14','h500a14']
    oros = ['h750a14','h1000a14']
    # oros = ['h1500a14','h2000a14','h3000a14','h4000a14']
    print oros, OUTPUTPATHS
    
    sms = ['60_homo']
    ENSMEAN = 'ensmean_day_3d_d0d5.nc'

    iparallel = 0
    idbg = 0
    
    if iparallel:    
        import argparse
        parser_ = argparse.ArgumentParser(description='Process an expname in expnames')
        parser_.add_argument('expname', action="store")
        myargs = parser_.parse_args()
        expn= myargs.expname
        print expn
        oros = [expn]


    surftopo = 'gauss'
    print "My plan is: "
    print oros
    print sms
    

    # iterate over the selected folder
    for OUTPUTPATH in OUTPUTPATHS:
        for oro in oros:
            for sm in sms:
                # Prepare directory structure information
                EXP = oro + '/' + sm + '/'
                srcpath = OUTPUTPATH + EXP + ENSMEAN
                tarpath = OUTPUTPATH + EXP + 'interp_'+ENSMEAN
                
                print "Source path:"
                print srcpath
                
                print "Target path:"
                print tarpath
                

                # open srnc/tarnc files
                srcnc = Dataset(srcpath,'a')
                tarnc = Dataset(tarpath,'w')
                
                if idbg:
                    srcnc.close()
                    tarnc.close()
                    continue


                # prepare dimensions of the new netcdf file
                # copy dimensions
                for dim in srcnc.dimensions: # of type 'unicode'
                    sz = srcnc.dimensions[dim].size # of type 'int'
                    tarnc.createDimension(dim, sz)  
                
                
                # Define the interpolation grid in the vertical
                lowerlev = np.arange(50,2000,100)
                higherlev = np.arange(2000,10001,1000)
                zinterp = np.concatenate((lowerlev,higherlev))
                
                nzi = len(zinterp)
                nrlat = srcnc.dimensions['rlat'].size
                nrlon = srcnc.dimensions['rlon'].size        
                nt = srcnc.dimensions['time'].size
                                     
                grid_int = np.zeros((nzi,nrlat,nrlon))
                i = 0
                for h in zinterp[::-1]:
                    grid_int[i,:,:] = h
                    i+=1
                    
                     
                # Save grid to tarnc
                tarnc.createDimension('lev_int',nzi)
                
                tarnc.createVariable('zinterp',float,dimensions=('lev_int'))            
                tarnc.createVariable('grid_int',float,
                                     dimensions=('lev_int','rlat','rlon'))
                for d in ('lon','lat'):
                    tarnc.createVariable(d, float, dimensions=('rlat','rlon'))
                    tarnc.variables[d][:] = srcnc.variables[d][:]
                
                tarnc.variables['zinterp'][:] = zinterp[::-1] # cosmo convection for z
                tarnc.variables['grid_int'][:] = grid_int
                
                
                #corrspnd = {'U':'Uinterp','V':'Vinterp','W':'Winterp',
                corrspnd = {'T':'Tinterp','QV':'QVinterp'} #  ,TODO QC, T, and more variables
             
                gridtype = {'U':'full', 'V':'full','W':'half',
                            'T':'full','QV':'full'}
                   
                
                tarnc.createVariable('MSE',float,dimensions=('time','lev_int','rlat','rlon'))     
                tarnc.createVariable('DSE',float,dimensions=('time','lev_int','rlat','rlon'))
                
                
                for var in corrspnd:
                    tarnc.createVariable(corrspnd[var],float,dimensions=('time','lev_int','rlat','rlon'))                
                    for i in range(nt):
                        print 'timestep: ' + str(i)
                        print sm, oro
                        print 'working on '+corrspnd[var]
                        srcfield = srcnc.variables[var][i,:]
                        tarfield = interpmethod(oro,srcfield,grid_int,gridtype=gridtype[var])
                        tarnc.variables[corrspnd[var]][i,:] = tarfield


                for i in range(nt):
                    tarnc.variables['DSE'][i,:] = tarnc.variables['Tinterp'][i,:] + g_g/c_p*tarnc.variables['grid_int'][:]
                    tarnc.variables['MSE'][i,:] = tarnc.variables['DSE'][i,:] + Lv/c_p*tarnc.variables['QVinterp'][i,:]
                        
                

                
                    
                
                
                srcnc.close()
                tarnc.close()


    
