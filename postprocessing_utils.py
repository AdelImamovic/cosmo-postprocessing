# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 14:44:17 2016

@author: adeli
"""
import sys
sys.path.append('/home/adeli/scripts/')
from analy_utils import ncvar2pyvar, nhalo, HHL_creator
from alt_colormaps import viridis, inferno, plasma, magma

import os
import matplotlib.pylab as plt
import numpy as np
import re

from math import sqrt
from joblib import Parallel, delayed




FIGPATH = 'testfigs/'


from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

def anchor_name(name):
    return AnchoredText(name,
          prop=dict(size=5), frameon=True,
          loc=2)





def ENT_flux(s1,s2,vecfield,vcoord=0):
    """ W -> E flux, S -> N, B -> T flux
    through simplex boundaries with origin at s1 and diagonal point at s2
    Arguments:
    -- s1: simplex origin coordinate
    -- s2: simplex diagonal coordinate
    -- vecfield:3D vector field
    """
    x1,y1,z1 = s1
    x2,y2,z2 = s2
    U, V, W = vecfield
    
    # simulate density profile to be exponential
    rho_dummy = 10.**(-vcoord / 10000.)


        
    Efl = (np.sum(U[x1,y1:y2,z1:z2],axis=0)-np.sum(U[x2,y1:y2,z1:z2],axis=0))*rho_dummy
    Efl1 = np.sum(U[x1,y1:y2,z1:z2],axis=0)
    Efl2 = np.sum(U[x2,y1:y2,z1:z2],axis=0)   
    plt.contourf(V[x1:x2,y2,z1:z2])
    assert 0    
    
    Nfl = (np.sum(V[x1:x2,y1,z1:z2],axis=0)-np.sum(V[x1:x2,y2,z1:z2],axis=0))*rho_dummy
    Tfl = np.sum(W[x1:x2,y1:y2,z1])-np.sum(W[x1:x2,y1:y2,z2])
    
    return Efl, Nfl, Tfl,Efl1,Efl2
    
    
def COSMOgridtoCartesiangrid(field):
    """ COSMO scalar fields are z',rlon,rlat coordinates
    z' indicates that highest model level has index 0.

    transpose axes to get rlon,rlat,z object. 
    highest model level has maximum index now
    index 0 for ground
    """
    field = np.swapaxes(field,0,2) # z <-> y ; field =field (y,x,z)
    field = np.swapaxes(field,0,1) # y <-> z ; field =field(x,y,z)
    field = field[:,:,::-1] # z <-> z'
    return field



def mask_field(field,threshold):
    """ Sets all field values fijk=field[i,j,k,..] to 0 if fijk < threshold.

    Arguments
     -- field: 
             np.array
     -- threshold:
             float
    """
    field[np.where(field)<np.mean(field)] = 0.
    return field

def plot_fluxes(datapath,time,f=False,ax=False):
    # extract dynamic variables & clip halo points
    U = ncvar2pyvar(datapath,'U')[time,:,nhalo:-nhalo,nhalo:-nhalo]
    V = ncvar2pyvar(datapath,'V')[time,:,nhalo:-nhalo,nhalo:-nhalo]    
    W = ncvar2pyvar(datapath,'W')[time,:,nhalo:-nhalo,nhalo:-nhalo]
    vcoord =ncvar2pyvar(datapath,'vcoord')[:1:-1]
    U = COSMOgridtoCartesiangrid(U)
    V = COSMOgridtoCartesiangrid(V)
    W = COSMOgridtoCartesiangrid(W)

        
    if not f:
        f, ax = plt.subplots(1,1)

    Efl,Nfl,Tfl,Efl1,Efl2 = ENT_flux((100,100,0),(200,200,49),(U,V,W),vcoord)
    print np.shape(Efl)
    ax.plot(Efl,vcoord,label='EW')
    ax.set_ylim([0,22000])
    ax.legend()
    
    
def circsym_mean_vec(field):
    """ Circular mean of field field(z',x,y) -> f(z,r)
    
    Paramters
    ----
        fxyz:
        c:
    Returns:
    ----
        frz
    TODO:
    ----
        translate grid so mountain top is staggered - center at 149.5?
    
    """
    U,V,W = field
    nz,nx,ny=U.shape
    c = nx/2.+0.5
    
    Urz=np.zeros((nz,nx/2))
    Wrz=np.copy(Urz)
    x = np.arange(nx)-c
    y = np.arange(ny)-c

    X,Y=np.meshgrid(x,y)
    R=np.sqrt(X**2+Y**2)
    PHI=np.arctan(Y/X)
    Uradial=np.zeros(U.shape)

    for k in range(nz):
        Uradial[k,:]=U[k,:]*np.cos(PHI)+V[k,:]*np.sin(PHI)
        Uradial[k,:,:nx/2]*=-1. #circular symmetry
    
    if 0: # debug
        lev=38
        cf=plt.contourf(Uradial[lev,:],cmap='coolwarm')
        cbar=plt.colorbar(cf)
        assert 0
        
    for k in range(nz):
        levelu = Uradial[k,:,:]
        levelw = W[k,:,:]
        for r in range(nx/2):
                Urz[k,r] = np.mean(levelu[np.logical_and(R>=r*1.,R<r+1.)])
                Wrz[k,r] = np.mean(levelw[np.logical_and(R>=r*1.,R<r+1.)])
                
    return Urz,Wrz
    
    
def circsym_mean_scal(field):
    """ cylindiric mean of a scalar field"""
    nz,nx,ny = field.shape
    frz = np.zeros((nz,nx/2))
    c = nx/2.+0.5
    x = np.arange(nx)-c
    y = np.arange(ny)-c
    X,Y=np.meshgrid(x,y)
    R=np.sqrt(X**2+Y**2)
    for k in range(nz):
        level = field[k,:,:]
        for r in range(nx/2):
            frz[k,r] = np.mean(level[np.logical_and(R>=r*1.,R<r+1.)])
    return frz
    
def circsym_mean_2D(field):
    """ cylindiric mean of a scalar field"""
    nx,ny = field.shape
    frz = np.zeros(nx/2)
    c = nx/2+0.5
    x = np.arange(nx)-c
    y = np.arange(ny)-c
    X,Y=np.meshgrid(x,y)
    R=np.sqrt(X**2+Y**2)
    for r in range(nx/2):
            frz[r] = np.mean(field[np.logical_and(R>=r*1.,R<r+1.)])
    return frz
    
    
    

def plot_uvw(datapath, time, f=False, ax=False, var = 'uw',comvarlev=42):
    """ Plots streamlines of an experiment
    
    var (string)
    --- string of concat variable names to plot
    
    testdata (string)
    --- absolute path of data 
    
    time (int)
    --- time slice under consideration
    
    TODO:
    --- generalize beyond var = 'UV'
    --- timeslice, rather time
    --- interpolation between staggered / collocated grids
    --- quiver low resolution output
    --- 
    """     

    imaskfield = True

    # center constants, window size
    c = 150
    sz = 100 / 2    
    
    # Prepare coordinate fields
    x = np.arange(c-sz,c+sz)
    y = np.arange(c-sz,c+sz)
    z = ncvar2pyvar(path,'vcoord')[-1::-1]        
    z = z[1:]
    
    Xy, Yx = np.meshgrid(x,y)
    Xz, Zx = np.meshgrid(x,z)  
    Xz = np.transpose(Xz)
    Zx = np.transpose(Zx)

    heights = {
        'flat':0.,
        '100m':100.,
        '250m':250.,
        '500m':500.
        }    
    
    # grid point heights
    height = [heights[h] for h in heights if re.search(h,datapath)][0]
    HHL = HHL_creator(Hm=height)
    Zx = HHL[c-sz:c+sz,153, 1:]
    Zx = Zx[:,::-1]

    
    # extract dynamic variables & clip halo points
    # TODO: interpolate to collocated grid
    U = ncvar2pyvar(datapath,'U')[time,:,nhalo:-nhalo,nhalo:-nhalo]
    V = ncvar2pyvar(datapath,'V')[time,:,nhalo:-nhalo,nhalo:-nhalo]    
    W = ncvar2pyvar(datapath,'W')[time,1:,nhalo:-nhalo,nhalo:-nhalo]
       
     # interpolate to grid
    U = interpolate2grid(U)
    V = interpolate2grid(V)
    W = interpolate2grid(W)
    
    # Convert to Cartesian grid (zlev',rlon,rlat) --> (rlon,rlat,zlev)
    U = COSMOgridtoCartesiangrid(U)
    V = COSMOgridtoCartesiangrid(V)
    W = COSMOgridtoCartesiangrid(W)
    
    # idebug
    print np.shape(U)
    print np.shape(V)
    print np.shape(W)
    print np.shape(Xz)
    print np.shape(Zx)   
    print np.shape(Xy)
    print np.shape(Yx) 
    
    # if f,ax hasn't been passed (default), then create it
    if (not f):
        f, ax = plt.subplots(1,1)   

    if var == 'uw': 
        imean = True

        # mean of U, W
        if imean:
            #not correct, need HHL as weights
            U = np.mean(U[:,148:152,:],axis=1) 
            W = np.mean(W[:,148:152,:],axis=1) 
                
        # window selection    
        U = U[c-sz:c+sz,:]
        W = W[c-sz:c+sz,:]
        speeduw=np.sqrt(U**2+W**2)   
        
        # mask vectors with magnitude smaller than threshold        
        if imaskfield:
            mask = np.zeros(speeduw.shape)
            mask[speeduw < 0.2]=True
            U = np.ma.masked_array(U,mask=mask)
            W = np.ma.masked_array(W,mask=mask)

 

        # contour fields of the (u,w)-speed       
        cf = ax.contourf(Xz,Zx,speeduw,cmap=viridis,vmin=0,vmax=10,alpha=0.5)   
        
        # vector fields
        Q = ax.quiver(Xz,Zx,U,W)
        qk = plt.quiverkey(Q, 2, 1.05, 0.5, r'$1 \frac{m}{s}$', labelpos='W',
                       fontproperties={'weight': 'bold'})
                       
        # cosmetics
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylim([0,12000])
        ax.set_xlim([100,201])
        ax.set_yticks(np.arange(0,12001,1000))
        #ax.set_ylabel('AGL height in m')
    
        f.savefig(FIGPATH+'test.pdf',bbox_inch = 'tight', dpi = 300)
        return cf
        
    if var == 'uv':
        imean = 0
        if imean:
            U = np.mean(U[c-sz:c+sz,c-sz:c+sz,comvarlev-2:comvarlev+2],axis=2)
            V = np.mean(V[c-sz:c+sz,c-sz:c+sz,comvarlev-2:comvarlev+2],axis=2)
            W = np.mean(W[c-sz:c+sz,c-sz:c+sz,comvarlev-2:comvarlev+2],axis=2)
        else:
            U = U[c-sz:c+sz,c-sz:c+sz,comvarlev]
            V = V[c-sz:c+sz,c-sz:c+sz,comvarlev]
            W = W[c-sz:c+sz,c-sz:c+sz,comvarlev]
            
        
        speeduv = np.sqrt(U**2+V**2)
        if imaskfield:
            mask = np.zeros(speeduv.shape)
            mask[speeduv < 0.2] = True
            U = np.ma.masked_array(U,mask=mask)
            V = np.ma.masked_array(V,mask=mask)
        
        cf =ax.contourf(Xy,Yx,Urz, cmap = viridis, alpha=0.5, vmin=-2, vmax=10)
        Q = ax.quiver(Xy,Yx,U,V)
        qk = plt.quiverkey(Q, 2, 1.05, 0.5, r'$1 \frac{m}{s}$', labelpos='W',
                       fontproperties={'weight': 'bold'})
        ax.set_aspect('equal')
        return cf
        

def plot_rzfields(oro):
    #OUTPUTPATH = '/lhome/adeli/Project_B/ensemblemean/'
    OUTPUTPATH='/hymet/adeli/project_B/300x300_7Kkmnowind_projectB/postprocessing/composites/'

    oro='250m'
    #runsel=['60_40_'+oro+'_7KKm_nowind_d15tod20.nc'
    runsel=['60_40','60_homo']
    
    x = np.arange(150)
    z = np.arange(50)    
    X,Z=np.meshgrid(x,z)
    Z=HHL_creator(Hm=250.0)
    Z=Z[150:,153,1:]
    Z=np.transpose(Z)
    Z=Z[::-1,:]

    xst = 1
    zst = 1
    
    timesel = [23,24]#25,26,27]#,22,23]
    nt=len(timesel)
    nr=len(runsel)
    
    f,ax=plt.subplots(nt,nr)
    
    for i in range(nt):
        for j in range(nr):
            datapath=OUTPUTPATH+oro+'/'+runsel[j]+'/'+'ensmean.nc'

            U = ncvar2pyvar(datapath,'U')[timesel[i],:,nhalo:-nhalo,nhalo:-nhalo]
            V = ncvar2pyvar(datapath,'V')[timesel[i],:,nhalo:-nhalo,nhalo:-nhalo]
            W = ncvar2pyvar(datapath,'W')[timesel[i],:,nhalo:-nhalo,nhalo:-nhalo]
            Urz,Wrz=circsym_mean_vec((U,V,W))
            Urz=Urz[::-1,:]
            Wrz=Wrz[::-1,:]
            #ax[i,0].plot(np.sum(Urz[:,20]*np.exp(-np.log(10)*Z[:,20]/10000))*Z[:,20],Z[:,20],label=runsel[j]+'im')
            #ax[i,1].plot(Urz[:,20]*np.exp(-np.log(10)*Z[:,20]/10000),Z[:,20],label=runsel[j])
            #ax[i,j].set_ylim([0,12000])
            speedrz=np.sqrt(Urz**2+Wrz**2)
            cf=ax[i,j].contourf(X,Z,speedrz,cmap=inferno, vmin=0,vmax=5)
            ax[i,j].set_xlim([0,80])
            ax[i,j].set_ylim([0,5000])
            Q=ax[i,j].quiver(X[::zst,::xst],Z[::zst,::xst],
                            Urz[::zst,::xst],Wrz[::zst,::xst],color='white')
            qk = plt.quiverkey(Q, 1, 1.05, 0.5, r'$1 \frac{m}{s}$', labelpos='W',
                       fontproperties={'weight': 'bold'})
    return f,ax


def interpolate2grid():
    x = np.arange(150)
    z = np.arange(50)    
    X,Z=np.meshgrid(x,z)
    Z=HHL_creator(Hm=500.0)
    Z=Z[150:,153,1:]
    Z=np.transpose(Z)
    Z=Z[::-1,:]
    
    OUTPUTPATH = '/lhome/adeli/Project_B/ensemblemean/'
    
    runsel = '60_homo_500m_7KKm_nowind_d15tod20.nc'
    datapath = OUTPUTPATH+runsel
    
    time = 20   
    
    U = ncvar2pyvar(datapath,'U')[time,:,nhalo:-nhalo,nhalo:-nhalo]
    #U = COSMOgridtoCartesiangrid(U)
    V = ncvar2pyvar(datapath,'V')[time,:,nhalo:-nhalo,nhalo:-nhalo]
    W = ncvar2pyvar(datapath,'W')[time,:,nhalo:-nhalo,nhalo:-nhalo]
    #V = COSMOgridtoCartesiangrid(V)

    Urz,Wrz=circsym_mean((U,V,W))
    Urz=Urz[::-1,:]

    xinterp = np.arange(0,100.,2)
    zinterp = np.arange(600,2000,100)  
    
    Xint, Zint = np.meshgrid(xinterp,zinterp)
    from scipy.interpolate import  griddata
    
    datapoints = np.array([X.flatten(),Z.flatten()])
    datapoints = np.transpose(datapoints)
    gridUrz = griddata(datapoints,Urz.flatten(),(Xint,Zint))
    gridWrz = griddata(datapoints,Wrz.flatten(),(Xint,Zint))
    plt.quiver(Xint,Zint,gridUrz,gridWrz)

def T_test():
    x = np.arange(150)
    z = np.arange(50)    
    X,Z=np.meshgrid(x,z)
    Z=HHL_creator(Hm=500.0)
    Z=Z[150:,153,1:]
    Z=np.transpose(Z)
    Z=Z[::-1,:]
    
    OUTPUTPATH = '/lhome/adeli/Project_B/ensemblemean/'    
    runsel = '60_homo_500m_7KKm_nowind_d15tod20.nc'
    datapath = OUTPUTPATH+runsel
    time =24
    Tref = ncvar2pyvar(datapath,'T')[time,:,nhalo:-nhalo,nhalo:-nhalo]
    Tpert = ncvar2pyvar(OUTPUTPATH+'60_80_500m_7KKm_nowind_d15tod20.nc','T')[time,:,nhalo:-nhalo,nhalo:-nhalo]
    Trzref = circsym_mean_scal(Tref) 
    Trzpert = circsym_mean_scal(Tpert)
    Trzpert = Trzpert[::-1,:]
    Trzref =Trzref[::-1,:]
    cf=plt.contourf(X,Z,Trzpert-Trzref,cmap='coolwarm')
    cbar=plt.colorbar(cf)


if __name__ == '__main__':
    f,ax=plot_rzfields(2)
    assert 0
    
    f,ax,cf = plot_rzfields(2)    
    
    OUTPUTPATH = '/lhome/adeli/Project_B/ensemblemean/'
    
    runsel = '60_homo_500m_7KKm_nowind_d15tod20.nc'
    datapath = OUTPUTPATH+runsel
    
    time = 20   
    
    U = ncvar2pyvar(datapath,'U')[time,:,nhalo:-nhalo,nhalo:-nhalo]
    #U = COSMOgridtoCartesiangrid(U)
    V = ncvar2pyvar(datapath,'V')[time,:,nhalo:-nhalo,nhalo:-nhalo]
    W = ncvar2pyvar(datapath,'W')[time,:,nhalo:-nhalo,nhalo:-nhalo]
    #V = COSMOgridtoCartesiangrid(V)

    Urz,Wrz=circsym_mean((U,V,W))
    Urz=Urz[::-1,:]
    Wrz=Wrz[::-1,:]


    x = np.arange(150)
    z = np.arange(50)    
    X,Z=np.meshgrid(x,z)
    Z=HHL_creator(Hm=500.0)
    Z=Z[150:,153,1:]
    Z=np.transpose(Z)
    Z=Z[::-1,:]
    print Z.shape,X.shape,Urz.shape
    f,ax=plt.subplots(1,1)
    #cf=plt.contourf(X,Z,Urz,cmap='coolwarm')
    cf=ax.contourf(X,Z,np.sqrt(Urz**2+Wrz**2),cmap='coolwarm')
    ax.quiver(X,Z,Urz,Wrz,alpha=0.5)
    cbar=plt.colorbar(cf)
    ax.set_xlim([0,100])
    ax.set_ylim([0,8000])
    
      
    
    
        
if __name__ == '__maain__':
    OUTPUTPATH = '/lhome/adeli/Project_B/ensemblemean/'
    
        
    timesl = (20,21,22,23,24)#23,24)#,24,26)#,24)
    ntimesl = len(timesl)
    
    allruns = os.listdir(OUTPUTPATH)        
    runsel = [run for run in allruns if (re.search('_250m_',run) and not re.search('homo',run))]
    #runsel = [runsel[0],runsel[-1]]
    runsel.sort()
    print runsel    
    runsel = [runsel[0],runsel[-1]]
    nruns = len(runsel)    
    idebug = False

    
    oro='250m'
    runsel=['60_40_'+oro+'_7KKm_nowind_d15tod20.nc',
            '60_50_'+oro+'_7KKm_nowind_d15tod20.nc',
    '60_homo_'+oro+'_7KKm_nowind_d15tod20.nc',
    '60_70_'+oro+'_nowind_d15tod20.nc',
    '60_80_'+oro+'_7KKm_nowind_d15tod20.nc',    
    ]
    
    
    
    runsel = '60_homo_500m_7KKm_nowind_d15tod20.nc'
    time = 22
    plot_fluxes(OUTPUTPATH+runsel,time)

    assert 0
    
    
    
    
    
    if idebug:  
        runsel = '60_30_250m_7KKm_nowind_d15tod20.nc'
        time = 22
        f, ax =plt.subplots(1,1)
        path = OUTPUTPATH+runsel
        cf = plot_uvw(path,time,f=f,ax=ax,var='uv',comvarlev=16)
        at=anchor_name('time = '+str(time / 2))            
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        ax.add_artist(at)   

        
    else:        
        path = OUTPUTPATH+runsel[0]
        nruns = len(runsel) 
    
        f, ax = plt.subplots(ntimesl,nruns,sharex=True,sharey=True)
        #Parallel(n_jobs=8)(delayed(plot_uvw)(OUTPUTPATH+runsel[j],timesl[i],f=f,ax=ax[i,j]) for j in range(nruns) for i in range(ntimesl))

        for j in range(nruns):
            path = OUTPUTPATH + runsel[j]  
            for i in range(ntimesl):
                cf = plot_uvw(path,timesl[i],f=f,ax=ax[i,j],var='uw')
                if i == 0:
                        at=anchor_name('time = '+str(timesl[i]/2))            
                        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
                        ax[i,j].add_artist(at)   
                #ax[i,j].set_title(timesl[i])
                
        vcoordvec=ncvar2pyvar(path,'vcoord')[:]    
        f.tight_layout()
        f.savefig(FIGPATH+oro+'_uw.pdf',bboxes_inch='tight')
        
    #cbaxes = f.add_axes([0.95,0.05,0.05,0.9]) 
    #cb=plt.colorbar(cf,orientation='vertical',cax=cbaxes)
    #cb.set_label('wind speed in m/s')
    #f.savefig('test.pdf',bboxes_inch='tight')