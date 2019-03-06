# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 09:37:39 2016

@author: adeli

""" 
""" Given an nc-file the circular mean around the center is calculated"""


import sys
sys.path.append('/home/adeli/scripts/python/')
from postprocessing_utils import circsym_mean_vec, circsym_mean_2D, circsym_mean_scal
from analy_utils import HHL_creator, nhalo
import numpy as np 
from netCDF4 import Dataset
import matplotlib.pylab as plt
import os

BASEB2 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/ALLVAR_3D/'
BASEB2_o3 = '/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/ALLVAR_3D/'
BASEA = '/hymet/adeli/project_A/256x256_7Kkmnowind_1km/postprocessing/composites/'


Rdry = 287.058 # gas constant of dry air
invRdry = 1. / Rdry

oroheight = {'flat' : 0.0, '125m':125.0, '250m' : 250.0, '500m' : 500.0}

hrange = [250,500,750,1000,1500,3000,4000] #range(250,1001,250)
arange = [5,7,10,14,20,25,30]
nh = len(hrange)
na = len(arange)

expnames=np.zeros((nh,na),dtype=object)

for i in range(nh):
    for j in range(na):
        expnames[i,j] = 'h'+str(hrange[i])+'a'+str(arange[j])


# cut off large parte of model boundaries to savememory?
isavemem = False      
if isavemem: nhalo = nhalo+128



if __name__ == '__main__':


    iparallel = False
    if iparallel:    
        import argparse
        parser = argparse.ArgumentParser(description='Process some integers.')
        parser.add_argument('integers', metavar='N', type=int, nargs='+',
                            help='an integer for the accumulator')
        args = parser.parse_args()
        myarange = args.integers
        print 'working on a = ' + str(myarange[0])
    else:
        myarange=range(na)

    myhrange = range(nh)


    iensmean = True #True #else seed
    
    if iensmean:
        srcflpref = 'LWP_ensmean_day_3d_d0d5.nc'
        tarflname = 'circmean_'+srcflpref
    else:
        srcflpref = 'seed'
        tarflname = 'circmean_seed_d0.nc'


    #myhrange,myarange=np.where(expnames=='h1000a7')
    print myarange,myhrange
    print expnames[:,myarange[0]]



    OUTPUTPATH = BASEB2_o3
    
    #select the files to iterate over, cross product is built oros X sms
    oros = expnames
    sms = ['60_homo']

    print OUTPUTPATH
   
    surftopo = 'gauss'
    


    for hi in [6]: #range(nh): #[3]: #myhrange:
        for hj in [3]: #range(na):#[4]:# myarange:
            
            # Prepare structure and look for the desired file
            # syntax: filename endswith d0.nc
            if surftopo in ['bell','cos2']: 
                EXP = expnames[hi,hj] + '_'+surftopo+'/60_homo/'
            else: #gauss
                EXP = expnames[hi,hj] + '/60_homo/'

            fls = os.listdir(OUTPUTPATH+EXP)
            flname = filter(lambda t: t.startswith(srcflpref),fls)
            # check if the desired file is available
            if flname == []:
                continue
            flname = flname[0]
	     
            # Define source and target file paths            
            srcpath = OUTPUTPATH+EXP+flname
            tarpath = OUTPUTPATH + EXP + tarflname #'circmean_day_d0d5.nc'
            print "Source path:"
            print srcpath
            print "Target path:"
            print tarpath
            
            srcnc = Dataset(srcpath,'r')
            tarnc = Dataset(tarpath,'w')
            
            # prepare dimensions of the new netcdf file
            # copy dimensions
            for dim in srcnc.dimensions:
                sz = srcnc.dimensions[dim].size
                tarnc.createDimension(dim, sz)

               
            # FIELDS TO POSTPROCESS    
            # DYNAMICS
            # 3D (nt x nz x nx x ny)
            print 'loading vars'
            
            varnames = filter(lambda x: x.startswith('LWP'),srcnc.variables.keys())
            print varnames

            TWP = np.zeros(srcnc.variables['LWP_QV'][:].shape)
            for var in varnames:
                data = srcnc.variables[var][:]
                TWP += data
                
            # NUMBER OF DIMENSIONS        
            nt,nx,ny=TWP.shape   
            
            # coordinates 
            print 'copying dimensions'
            tarnc.createDimension('r',nx/2) 
            x = np.arange(nx/2)
          
            # create r coordinate
            r = tarnc.createVariable('r',float,dimensions=('r'))
            r[:] = x
                 
                 
            # Circular means
            print 'create new vars'


            #create new variable names and assign circmean
            for varn in varnames:
                print varn
                circvarn = varn+'r'
                oldvar = srcnc.variables[varn][:]               
                circmeanvar = tarnc.createVariable(circvarn,float,dimensions=('time','r'))
                for i in range(nt):
                    circmeanvar[i,:] = circsym_mean_2D(oldvar[i,:])
                    print i
            
            # add TWP too:
            circmeanTWP = tarnc.createVariable('TWPr',float,dimensions=('time','r'))
            for i in range(nt):
                circmeanTWP[i,:] = circsym_mean_2D(TWP[i,:])

            
              
            # add E and H as well
            srcnc.close()
            srcpath =  OUTPUTPATH+EXP+'ensmean_day_3d_d0d5.nc'
            
            
            srcnc = Dataset(srcpath,'r')
            oldnewvar = {'EFLUX':'EFLUXr','HFLUX':'HFLUXr'}


            for oldvarn in oldnewvar:
                print oldvarn
                newvarn = oldnewvar[oldvarn]
                
                nhalob= nhalo
                olddata = srcnc.variables[oldvarn][:,-1,nhalob:-nhalob,nhalob:-nhalob]
                newdata = tarnc.createVariable(newvarn,float,dimensions=('time','r'))
                
                
                for i in range(nt):
                    print i
                    newdata[i,:] = circsym_mean_2D(olddata[i,:])

            
            


            tarnc.close()
            
    
