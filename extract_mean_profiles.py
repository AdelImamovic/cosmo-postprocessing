"""
    Script extracts hourly domain-mean profiles of
    
    T,P,QV,QC,U,V,W
    
    from a COSMO run.
   
    and saves it into a new file.
    
"""


from netCDF4 import Dataset
import numpy as np
import os



if __name__=='__main__':
    
    varnames = ['T','P','QV','QC','U','V']
    
    tarvars = {varn : np.zeros(48,50) for varn in varnames}
    
    #create target netcdf file atmos.nc
    basepath = '/users/adeli/'
    tarnc = Dataset(basepath+'atmos.nc','w')
    
    
    #open src nc
    srcpath = '/scratch/snx3000/adeli/project_B2/512x512_7Kkmnowind_1km'
    srcpath +='/rawfiles/h0a0_test/60_homo/seed76996/output/'
    
    coutputfiles = os.listdir(srcpath)
    
    #extract hourly mean fields
    # extract first day only
    filter(lambda x: x.startswith('lfff01'),outputfiles)
    
    coutputfiles.sort()
    
    srcnc0 =Dataset(srcpath+coutputfiles[0])
    
    #create dimensions
    for dim in srcnc0.dimensions: # of type 'unicode'
        sz = srcnc0.dimensions[dim].size # of type 'int'
        tarnc.createDimension(dim, sz)  
    
    for var in varnames:
        tarnc.createVariable(var,float,dimensions=('time','lev'))
    
    ti = 0 
    for coutfl in coutputfiles:
        # iterate over variables
        srcnc = Dataset(srcpath+coutfl,'r')
        for varn in varnames:
            # extract domain mean
            var = srcnc.variables[varn][:]
            dommeanprof = np.mean(var,axis=(2,3))
            
            # write to target
            tarnc.variables[var][ti,:] = dommeanprof
            ti+=0
        srcnc.close()
        
    # close target netcdf file    
    tarnc.close()
    
    
    