# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 14:42:43 2016

@author: adeli
"""

#import COSMOoutputanalysis_utils as utls
import os,re
from analy_utils import ncvar2pyvar
import matplotlib.pylab as plt
import numpy as np
from collections import OrderedDict

nhalo=3
center=150

fig=plt.figure()
ax=plt.axes()

#plt.set_aspect('equal')
#matplotlib.style.use('ggplot')


ATM='7Kkm_wind'
ANALYSISPATH='/lhome/adeli/Project_A/output/'+ATM+'/'
FIGPATH='/lhome/adeli/Project_A/scripts/figures/'

oro='500m'
sm='60_homo'
smhet='radius'
timeindex=72
listwithruns=os.listdir(ANALYSISPATH)
allruns=()
allruns=[run for run in listwithruns if (re.search(oro,run) and run.endswith('.nc'))]

curruns=[]

for run in allruns:
    if re.search(sm,run) or re.search(smhet,run):    
        curruns.append(run)

size=['min','max','median']

minmaxmean={}
plt.hold(True)

for currun in curruns:
    totprecfields=ncvar2pyvar(ANALYSISPATH+currun,'TOT_PREC')
    cumtotprecfield=np.cumsum(totprecfields,axis=0)
    cumtotprecfieldday1=cumtotprecfield[timeindex][nhalo:-nhalo,nhalo:-nhalo]
    
    
    distances=np.arange(2,151,2)
    
    distanceindices=distances[::-1]
    
    totprec_dist=[]
    for dist in distances:
        totprec_dist.append(np.mean(cumtotprecfieldday1[(center-dist):(center+dist),(center-dist):(center+dist)]))
        
    if re.search('60_homo',currun):   
        plt.plot(2*distances,totprec_dist,color='green',alpha=0.8,label='HOM60')#,label=currun[:20])
        #plt.hlines(totprec_dist[25],0,140,color='green')        
    if re.search('40_homo',currun):   
        plt.plot(2*distances,totprec_dist,color='brown',alpha=0.8,label='HOM40')
        #plt.hlines(totprec_dist[25],0,140,color='brown') 
    if re.search('80_homo',currun):   
        plt.plot(2*distances,totprec_dist,color='blue',alpha=0.8,label='HOM80')    
        #plt.hlines(totprec_dist[25],0,140,color='blue') 
    if re.search('radius',currun):
        if (re.search('60_70',currun) or re.search('60_80',currun) or re.search('60_90',currun)):
            plt.plot(2*distances,totprec_dist,color='black',alpha=0.8,label='NEGA')
            continue
        plt.plot(2*distances,totprec_dist,color='red',alpha=0.8,label='POSA')
        #plt.hlines(totprec_dist[25],0,140,color='red')


handles, labels1 = plt.gca().get_legend_handles_labels()

labelsutf8=[l.encode('UTF8') for l in labels1]

by_label = OrderedDict(zip(labelsutf8, handles))

plt.ylim([0,15])
plt.title(ATM+' and '+oro)
plt.xlabel('averaging window length in km')
plt.ylabel('rain amount in mm after '+str(timeindex/2)+'h')    
plt.savefig(FIGPATH+'time'+str(timeindex/2)+smhet+'_'+ATM+'_'+sm+'_'+oro+'.pdf',dpi=300,bboxes_inch='tight')
