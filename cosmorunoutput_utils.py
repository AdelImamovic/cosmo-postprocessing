# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:03:15 2016

@author: adeli
"""

from ensemble_utils import ExperimentEnv
from analy_utils import ncvar2pyvar
import os
import numpy as np
import matplotlib.pylab as plt
import re
import scipy
import style_sheet as s_s 



def interpolate2coll(vari, vert=1):
    levels=np.arange(1,len(vari))
    levelsp1=np.arange(1,len(vari)+1)
    vari=scipy.interp(levels,levelsp1,vari)  
    return vari
    
def interpolate2stag(self,vert=1):
    pass


class COSMORunOutput(ExperimentEnv):
    
    def __init__(self,path,runname):
        self.runname=runname
        self.datapath=path+runname+'/'
        self.files=os.listdir(self.datapath)
        self.refpath='/lhome/adeli/Project_A/output/references/'
        

        for fl in self.files:
            if not (fl.startswith('lfff') and fl.endswith('.nc')): 
                self.files.remove(fl)
                continue

        (self.files).sort()

        
        
        
        #determine the timeslice indices
        self.timeslices=[]
        self.hourssincestart=[]
        self.gen_timeslice_indices()
        #determine the ensemble member
        self.tag=()
        self.gettag()
        self.refpath+=self.tag[0]+'/'
        
        self.plotcol=ExperimentEnv._smcols[self.tag[1]]
        
        #get initial profiles of run
        self.initprofiles={'U':[],'T':[],'RELHUM':[],'HHL':[]}
        #self.getinitprofiles()
        
        
    def gettag(self):
        """Determines the tag of the Runoutput. Tag=(oro,sm)."""        
        for oro in ExperimentEnv._oros:
            if re.search(oro, self.runname):
                curoro=oro
                break
        for sm in ExperimentEnv._sms:
            if re.search(sm, self.runname):
                cursm=sm
                break
        self.tag=curoro,cursm
      
        
    def getinitprofiles(self,iplot=False):
        """Generates init properties from non-changing initfile."""
        reffilepath=''#TODO
        self.height_levels=ncvar2pyvar(reffilepath+'lfff00000000.nc','HHL')
        reffilename='lfff00000000.nc'
        
        for vari in self.initprofiles:
            profile=ncvar2pyvar(self.datapath+reffilename,vari)[0,:,0,0]
            self.initprofiles[vari]=profile
        #interpolate to collocated grid
        self.initprofiles['HHL']=interpolate2coll(self.initprofiles['HHL'])
                    
        if iplot:
            f,ax=plt.subplots(1,3)
            i=0
            for vari in self.initprofiles:
                if vari=='HHL':
                    continue
                ax[i].plot(self.initprofiles[vari],self.initprofiles['HHL'],label=vari)
                ax[i].set_xlabel(vari)
                ax[i].set_ylim([500,12000])
                #ax[i].set_aspect()
                if i==0:                
                    ax[i].set_yticks(np.arange(1000,12000,1000))
                else:
                    ax[i].set_yticks(np.arange(1000,12000,1000),[])
                
                i+=1
            f.tight_layout()

    @staticmethod            
    def time2COSMOfilename(time):
        return 'lfff00'+str(time)+'00.nc'
        
    @staticmethod
    def time2outputindex(time):
        """ Assume format."""        
        return 26 #for testing
        
    def gen_timeslice_indices(self):
        """Loops over all the filenames to generate the extract timeindices.
        and to create hours since start."""
        for filename in self.files:
            days=int(filename[4:6])
            hours=int(filename[6:8])
            minutes=int(filename[8:10])
            (self.timeslices).append(str(filename[4:10]))
            #(self.hourssincestart).append(24.*days+hours+(minutes/30)*0.5)
        self.hourssincestart=np.linspace(0,24,len(self.timeslices))            
            
    def plot_uv_w(self,time=1300,level=40,f=False,ax=False,figname='test.pdf',save=False):
        if (not f):
            f,ax=plt.subplots(1,1)#self.illustrate_domain_properties()
            ax.set_aspect('equal')
        
        filename=self.time2COSMOfilename(time)
        
        try:
            HHL=ncvar2pyvar(self.datapath+filename,'HHL')[0,level,1,1]
        except:
            print 'no HHL in file'
        background={'U0':[]}
        variables={'U':[],'V':[],'W':[]}
        dims={'lat':[],'lon':[]}
        print 'here'
        for dim in dims:
            print self.datapath+filename
            dims[dim]=ncvar2pyvar(self.datapath+filename,dim)[nhalo:-nhalo,nhalo:-nhalo]*np.pi/180.*6371.
         
        for vari in variables:
            
            variables[vari]=np.transpose(ncvar2pyvar(self.datapath+filename,vari)[0,level,nhalo:-nhalo,nhalo:-nhalo])
            
        for back in background:

            background[back]=ncvar2pyvar(self.datapath+'lfff00000000.nc','U')[0,level,nhalo:-nhalo,nhalo:-nhalo]
            
        #f,ax=self.illustrate_domain_properties(save=False)
       
        U=variables['U']
        V=variables['V']
        W=variables['W']
        #U=U-background['U0']
        speed=np.sqrt(U**2+V**2)
        
        UN=U/speed
        VN=V/speed    
        
        
        #X,Y=np.meshgrid(dims['rlon'],dims['rlat'])
        X,Y=dims.values()              
        #ax.streamplot(X,Y,U,V,color='Red',arrowsize=0.1)
        #ax[1].quiver(X,Y,UN,VN,color='Red',headlength=7)
        ax.contourf(X,Y,W,cmap=my_green_cm,vmin=1.0)
        #ax.set_title(str(self.height_levels[level])+' m')
        #ax.text(150,250,'max speed='+str(np.max(speed))+' m/s')
        ax.set_xlim([75,225])
        ax.set_ylim([75,225])
        #for i in range(2): 
        #ax.set_aspect('equal')    
        f.savefig(figname,bboxes_inch='tight')
        return f,ax
        
        
    def plot_uw_v(self,time,rlat):
        pass
        
    def plot_vw_u(self,time,rlon):
        pass
    
    #this functions could be realised as decorators (do the looping over the files)
    def dc_uhmax2(self,vari='U',hourslim=False):
        """plots dc of vari in hourslim interval."""
        backgr=ncvar2pyvar(self.datapath+'lfff00000000.nc',vari)[0,:,:,:]
     
        i=0
        u_max=self.hourssincestart*1.
        u_prime=u_max*1.
        v_prime=u_max*1.
        speed=u_max*1.
        for filename in self.files:
            print filename
            U_prime=ncvar2pyvar(self.datapath+filename,'U')[0,:,:,:]-backgr
            V_prime=ncvar2pyvar(self.datapath+filename,'V')[0,:,:,:]
            speed[i]=np.max(np.sqrt(U_prime[:,103:203,103:203]**2+V_prime[:,103:203,103:203]**2))
            u_prime[i]=np.max(U_prime[:,103:203,103:203])
            v_prime[i]=np.max(V_prime[:,103:203,103:203])
            i+=1
        
        #plt.plot(self.hourssincestart,u_max)
        plt.xlabel('hours since start')
        plt.ylabel('horizontal u in m/s')
        #plt.xlim([0,24])
        return u_prime,v_prime,speed
    
    def dc_wmax(self):
        w_max=np.array(self.hourssincestart)*1.
        i=0        
        for filename in self.files:
            print self.datapath+filename
            w=ncvar2pyvar(self.datapath+filename,'W')[:,:,:,:]
            print np.shape(w)
        
            w_max[i]=np.max(w)
            i+=1
            #assert 0
        return w_max
        
    def dc_uhmax(self):
        uh_max=np.array(self.hourssincestart)*1.
        i=0        
        for filename in self.files:
            print self.datapath+filename
            u=ncvar2pyvar(self.datapath+filename,'U')[:,:,:,:]-ncvar2pyvar(self.refpath+'lfff00000000.nc','U')
            v=ncvar2pyvar(self.datapath+filename,'V')[:,:,:,:]
            uh=np.sqrt(u**2+v**2)
        
            uh_max[i]=np.max(uh)
            i+=1
            #assert 0
        return uh_max
        
    
    def calc_dc(var):
        """....

        TODO:
        ----
            implement dispatch method.
            
        """
        pass
        
    def dc_2D(self,vari='ALHFL_S'):
        
        dc_sh=(self.hourssincestart)*1.
        i=0        
        for filename in self.files:
            #print filename
            dc_sh[i]=-1*np.mean(ncvar2pyvar(self.datapath+filename,'ALHFL_S')[0,103:203,103:203])
            i+=1
        return dc_sh



        
def plot_dc(path,vari='ALHFL_S',ttl='flat'):
    fig,axs=plt.subplots(1,1)
    runs=os.listdir(path)
    for run in runs:
        run+='/'
        temp=COSMORunOutput(path,run)
        data=temp.dc_2D(vari=vari)
        axs.plot(temp.hourssincestart,data,color=temp.plotcol,label=temp.tag[1])
        axs.text(temp.hourssincestart[np.where(data==np.max(data[0:45]))],np.max(data[0:45]),temp.tag[1],fontsize=8,color=temp.plotcol)
    axs.set_ylabel('LH surface flux in W/m')
    axs.set_xlabel('local time in hours')
    axs.set_xticks(np.arange(5)*6)
    y0,y1=axs.get_ylim()
    axs.set_title(ttl)
       
    axs.plot((12,12),(y0,y1),color='red')
    axs.set_xlim([0,24])
    return fig,axs
    
if __name__=='__main__':    
    idbg=True
    oro='100m'
    path='/lhome/adeli/Project_A/output/8Kkm_wind/relaxoff/output_full/'+oro+'/'    
    runs=os.listdir(path)
    run=runs[12]
    #for run in runs:
        
    temp=COSMORunOutput(path,run) 
#    f,ax=temp.plot_uv_w(time=1000,level=33)
#    ax.set_title(temp.tag[0]+' '+temp.tag[1])
#    f.show()
    #temp.plot_uv_w()
    

    
#    if not idbg:
#        fig,ax=plot_dc(path,ttl=oro)
#        figname=oro+'_LH.pdf'
#        fig.savefig(s_s.FIGPATH+figname)


#    runs=os.listdir(path)
#    run=runs[0]
#    f,ax=plt.subplots(1,1)
#    for run in runs:
#        
#        temp=COSMORunOutput(path,run+'/')
#        print temp.hourssincestart
#        ax.plot(temp.hourssincestart,temp.dc_2D(),color=temp.plotcol)

  
    #fig,axs=plot_dc(path,vari='ALHFL_S')
    
    #fig.savefig(FIGPATH+'LH.pdf')
    
    #Figure 4
        
        
#    #Figure 10 - updraf velocities
#    f10,ax10=plt.subplots(1,1)
#    for run in runs:
#        temp=COSMORunOutput(path,run)
#        data=temp.dc_wmax()
#        ax10.plot(temp.hourssincestart,data,color=temp.plotcol)
#        
#    ax10.set_xlabel('local time in hours')
#    ax10.set_ylabel('maximum w in m/s')
#    ax10.set_title(oro)    
#    f10.savefig(s_s.FIGPATH+'wmax'+oro+'.pdf')
    
        #Figure 10 - updraf velocities
    f10,ax10=plt.subplots(1,1)
    for run in runs:
        temp=COSMORunOutput(path,run)
        data=temp.dc_uhmax()
        ax10.plot(temp.hourssincestart,data,color=temp.plotcol)
        
    ax10.set_xlabel('local time in hours')
    ax10.set_ylabel('maximum u_horiz in m/s')
    ax10.set_title(oro)    
    f10.savefig(s_s.FIGPATH+'uhmax'+oro+'.pdf')
    
    
        
        
