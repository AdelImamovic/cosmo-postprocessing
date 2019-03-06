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
OUTPUTPATH = BASEB2_o3

Rdry = 287.058 # gas constant of dry air
invRdry = 1. / Rdry

oroheight = {'flat' : 0.0, '125m':125.0, '250m' : 250.0, '500m' : 500.0}

hrange = [250,500,750,1000,1500,2000, 3000,4000] #range(250,1001,250)
arange = [5,7,10,14,20,25,30]
nh = len(hrange)
na = len(arange)

expnames=np.zeros((nh,na),dtype=object)

for i in range(nh):
    for j in range(na):
        expnames[i,j] = 'h'+str(hrange[i])+'a'+str(arange[j])


isavemem = True       
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
        srcflpref = 'ensmean'
        tarflname = 'circmean_day_d0d5.nc'
    else:
        srcflpref = 'seed'
        tarflname = 'circmean_seed_d0.nc'


    #myhrange,myarange=np.where(expnames=='h1000a7')
    print myarange,myhrange
    print expnames[:,myarange[0]]



    
    
    #select the files to iterate over, cross product is built oros X sms
    oros = expnames
    sms = ['60_homo']

    print OUTPUTPATH
   
    surftopo = 'gauss'
    


    for hi in [5,7]: #myhrange:
        for hj in [3]: #myarange:
            
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
            print flname
	     
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
            print U.shape
            # coordinates 
            print 'copying dimensions'
            tarnc.createDimension('r',nx/2) 
            x = np.arange(nx/2)
            z = np.arange(nz)
            X, Z0 = np.meshgrid(x,z)
            #Z = HHL_creator(oroheight[oro],nx,ny,nz+1)
            Z=HHL_creator(hrange[hi],nx,ny,nz+1,arange[hj],surftopo=surftopo)
            

            print Z.shape
            #Z=Z[nx/2:,nx/2,1:]
            Z = Z[1:,nx/2,nx/2:]
            print Z.shape
            #Z=np.transpose(Z)
          

            Xnc = tarnc.createVariable('X',float,dimensions=('level','r'))
            Znc = tarnc.createVariable('Z',float,dimensions=('level','r'))                   
            Xnc[:] = X
            Znc[:] = Z
            print Z.shape
                  
            # Circular means
            print 'create new vars'
            Urz = tarnc.createVariable('Urz',float,dimensions=('time','level','r'))
            Wrz = tarnc.createVariable('Wrz',float,dimensions=('time','level','r'))
            Trz = tarnc.createVariable('Trz',float,dimensions=('time','level','r'))
            Prz = tarnc.createVariable('Prz',float,dimensions=('time','level','r'))
            RHOrz = tarnc.createVariable('RHOrz',float,dimensions=('time','level','r'))
            speedrz = tarnc.createVariable('Speedrz',float,dimensions=('time','level','r'))

            QVrz = tarnc.createVariable('QVrz',float,dimensions=('time','level','r'))    
            QRrz = tarnc.createVariable('QRrz',float,dimensions=('time','level','r'))
            QSrz = tarnc.createVariable('QSrz',float,dimensions=('time','level','r'))
            QIrz = tarnc.createVariable('QIrz',float,dimensions=('time','level','r'))
            QCrz = tarnc.createVariable('QCrz',float,dimensions=('time','level','r'))
            QGrz = tarnc.createVariable('QGrz',float,dimensions=('time','level','r')) 
            
            # TODO : DESTAGGER
            HFLUXrz = tarnc.createVariable('HFLUXrz',float,dimensions=('time','level1','r')) 
            EFLUXrz = tarnc.createVariable('EFLUXrz',float,dimensions=('time','level1','r')) 
            
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
            print 'calculating circular means'
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
    
