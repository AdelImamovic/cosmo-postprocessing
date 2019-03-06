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

# Define the path that contains the ensemble fields
#BASEA_OLD = '/hymet/adeli/project_A/300x300_8Kkmwind_1km/postprocessing/composites/'
#BASEA2_OLD = '/hymet/adeli/project_A/384x384_8Kkmnowind_1km/postprocessing/composites/'
BASEA = '/net/o3/hymet/adeli/project_A/256x256_7Kkmnowind_1km/postprocessing/composites/'

BASEwjan = '/net/o3/hymet/wjan/for_wjan/Projecta/'

Rdry = 287.058 # gas constant of dry air
invRdry = 1. / Rdry



#
smsall = ['60_30','60_40','60_50','60_70','60_80','60_90',
          '40_homo','60_homo','80_homo'] 
smsensi = ['40_10','40_20','40_30','40_50','40_60','40_70']

smprojb_wjan = ['60_40_radius_'+str(n) for n in range(5,41,5)]



#
orosall = ['flat','125m','250m','500m']
oroheight = {'flat' : 0.0, '125m':125.0, '250m' : 250.0, '500m' : 500.0}



if __name__ == '__main__':
    OUTPUTPATH = BASEA
    
    #select the files to iterate over, cross product is built oros X sms
    oros = ['flat']
    sms = ['10_homo','20_homo','30_homo','50_homo','70_homo']

    print OUTPUTPATH
    print sms
    sms = ['60_30']

    for oro in oros:
        for sm in sms:
            
            # Prepare directory structure information
            EXP = oro + '/' + sm + '/'
            srcpath = OUTPUTPATH + EXP + 'ensmean_day_d1d5.nc'#'ensmean_day_d1d5.nc'
            tarpath = OUTPUTPATH + EXP + 'circmean_day_d1d5.nc'
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
            U = srcnc.variables['U'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            V = srcnc.variables['V'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            W = srcnc.variables['W'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            T = srcnc.variables['T'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            P = srcnc.variables['P'][:,:,nhalo:-nhalo,nhalo:-nhalo]     
            
            # PHYSICS
            # 3D (nt x nz x nx x ny)        
            QV = srcnc.variables['QV'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            QR = srcnc.variables['QR'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            QC = srcnc.variables['QC'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            QI = srcnc.variables['QI'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            QS = srcnc.variables['QS'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            QG = srcnc.variables['QG'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            # 2D (nt x nx x ny)
            TOT_PREC = srcnc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo]
            
            # DIAGNOSTICS
            # 3D (nt x nz x nx x ny)
            HFLUX = srcnc.variables['HFLUX'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            EFLUX = srcnc.variables['EFLUX'][:,:,nhalo:-nhalo,nhalo:-nhalo]
            # 2D (nt x nx x ny)
            HPBL = srcnc.variables['HPBL'][:,nhalo:-nhalo,nhalo:-nhalo]
            TQC = srcnc.variables['TQC'][:,nhalo:-nhalo,nhalo:-nhalo]
            CAPE_ML = srcnc.variables['CAPE_ML'][:,nhalo:-nhalo,nhalo:-nhalo]
            CIN_ML = srcnc.variables['CIN_ML'][:,nhalo:-nhalo,nhalo:-nhalo]  
            CAPE_MU = srcnc.variables['CAPE_MU'][:,nhalo:-nhalo,nhalo:-nhalo]
            CIN_MU = srcnc.variables['CIN_MU'][:,nhalo:-nhalo,nhalo:-nhalo]
            LCL_ML = srcnc.variables['LCL_ML'][:,nhalo:-nhalo,nhalo:-nhalo]
            LFC_ML = srcnc.variables['LFC_ML'][:,nhalo:-nhalo,nhalo:-nhalo]            
        

            VAR_3D_vec = {}
            VAR_3D_scal = {}
            VAR_2D_scal = {}


            # NUMBER OF DIMENSIONS        
            nt,nz,nx,ny=U.shape        

            # coordinates 
            tarnc.createDimension('r',nx/2) 
            x = np.arange(nx/2)
            z = np.arange(nz)
            X, Z0 = np.meshgrid(x,z)
            Z = HHL_creator(oroheight[oro],nx,ny,nz+1)
            Z=Z[nx/2:,nx/2,1:]
            Z=np.transpose(Z)            
            Xnc = tarnc.createVariable('X',float,dimensions=('lev','r'))
            Znc = tarnc.createVariable('Z',float,dimensions=('lev','r'))                   
            Xnc[:] = X
            Znc[:] = Z
        
        
        
        
            # Circular means
            Urz = tarnc.createVariable('Urz',float,dimensions=('time','lev','r'))
            Wrz = tarnc.createVariable('Wrz',float,dimensions=('time','lev','r'))
            Trz = tarnc.createVariable('Trz',float,dimensions=('time','lev','r'))
            Prz = tarnc.createVariable('Prz',float,dimensions=('time','lev','r'))
            RHOrz = tarnc.createVariable('RHOrz',float,dimensions=('time','lev','r'))
            speedrz = tarnc.createVariable('Speedrz',float,dimensions=('time','lev','r'))

            QVrz = tarnc.createVariable('QVrz',float,dimensions=('time','lev','r'))    
            QRrz = tarnc.createVariable('QRrz',float,dimensions=('time','lev','r'))
            QSrz = tarnc.createVariable('QSrz',float,dimensions=('time','lev','r'))
            QIrz = tarnc.createVariable('QIrz',float,dimensions=('time','lev','r'))
            QCrz = tarnc.createVariable('QCrz',float,dimensions=('time','lev','r'))
            QGrz = tarnc.createVariable('QGrz',float,dimensions=('time','lev','r')) 
            
            # TODO : DESTAGGER
            HFLUXrz = tarnc.createVariable('HFLUXrz',float,dimensions=('time','lev_2','r')) 
            EFLUXrz = tarnc.createVariable('EFLUXrz',float,dimensions=('time','lev_2','r')) 
            
            TOT_PRECr = tarnc.createVariable('TOT_PRECr',float,dimensions=('time','r'))
            HPBLr = tarnc.createVariable('HPBLr',float,dimensions=('time','r'))
            TQCr = tarnc.createVariable('TQCr',float,dimensions=('time','r'))
            CAPE_MLr = tarnc.createVariable('CAPE_MLr',float,dimensions=('time','r'))
            CIN_MLr = tarnc.createVariable('CIN_MLr',float,dimensions=('time','r'))          
            CAPE_MUr = tarnc.createVariable('CAPE_MUr',float,dimensions=('time','r'))
            CIN_MUr = tarnc.createVariable('CIN_MUr',float,dimensions=('time','r')) 
            LCL_MLr = tarnc.createVariable('LCL_MLr',float,dimensions=('time','r'))
            LFC_MLr = tarnc.createVariable('LFC_MLr',float,dimensions=('time','r'))




            #_3Drzvars = {name+'rz': tarnc.createVariable(name,float,dimensions=('time','lev','r')) for name in _3Dvarnames}
            #_2Drzvars = {name+'r': tarnc.createVariable(name,float,dimensions=('time','r')) for name in _2Dvarnames}            
            
            
            # TODO: Optimise loops for locality exploitation
            for i in range(nt): 
                print EXP + ' timestep=' + str(i)
                # 3D_vec
                Urz[i,:],Wrz[i,:] = circsym_mean_vec((U[i,:],V[i,:],W[i,:]))
                speedrz[i,:] = np.sqrt(Urz[i,:]**2+Wrz[i,:]**2)
                
                # 3D_scal
                Trz[i,:] = circsym_mean_scal(T[i,:])
                Prz[i,:] = circsym_mean_scal(P[i,:])
                RHOrz[i,:] = (Prz[i,:] / Trz[i,:]) * invRdry 
                QVrz[i,:] = circsym_mean_scal(QV[i,:])
                QSrz[i,:] = circsym_mean_scal(QS[i,:])
                QIrz[i,:] = circsym_mean_scal(QI[i,:])
                QCrz[i,:] = circsym_mean_scal(QC[i,:])
                QGrz[i,:] = circsym_mean_scal(QG[i,:])
                QRrz[i,:] = circsym_mean_scal(QR[i,:])
                EFLUXrz[i,:] = circsym_mean_scal(EFLUX[i,:])
                HFLUXrz[i,:] = circsym_mean_scal(HFLUX[i,:])
                
                # 2D
                TOT_PRECr[i,:] = circsym_mean_2D(TOT_PREC[i,:])
                HPBLr[i,:] = circsym_mean_2D(HPBL[i,:])
                TQCr[i,:] = circsym_mean_2D(TQC[i,:])
                CAPE_MLr[i,:] = circsym_mean_2D(CAPE_ML[i,:])
                CIN_MLr[i,:] = circsym_mean_2D(CIN_ML[i,:])
                CAPE_MUr[i,:] = circsym_mean_2D(CAPE_MU[i,:])
                CIN_MUr[i,:] = circsym_mean_2D(CIN_MU[i,:])
                LCL_MLr[i,:] = circsym_mean_2D(LCL_ML[i,:])
                LFC_MLr[i,:] = circsym_mean_2D(LFC_ML[i,:])
                
                
                
                         
                
            tarnc.close()
            srcnc.close()
    
