
# coding: utf-8

# In[ ]:


import numpy as np
import pyJHTDB
from pyJHTDB import libJHTDB


# In[ ]:


import json
with open("parameters-getCutout.txt","r") as pf:
    params=json.load(pf)

auth_token=params["token"]
tstart=int(params.get("tstart"))
tend=int(params.get("tend"))
tstep=int(params.get("tstep"))
xstart=int(params.get("xstart"))
ystart=int(params.get("ystart"))
zstart=int(params.get("zstart"))
xend=int(params.get("xend"))
yend=int(params.get("yend"))
zend=int(params.get("zend"))
xstep=int(params.get("xstep",1))
ystep=int(params.get("ystep",1))
zstep=int(params.get("zstep",1))
Filter_Width=int(params.get("Filter_Width",1))
time_step=int(params.get("time_step",0))
fields=params.get("fields","u")
data_set=params.get("dataset","isotropic1024coarse")
output_filename=params.get("output_filename",data_set)

#if fields == 'u':
#    VarName="Velocity"
#    dim = 3
#elif fields == 'a':
#    VarName="Vector Potential"
#    dim = 3
#elif fields == 'b':
#    VarName="Magnetic Field"
#    dim = 3
#elif fields == 'p':
#    VarName="Pressure"
#    dim = 1
#elif fields == 'd':
#    VarName="Density"
#    dim = 1
#elif fields == 't':
#    VarName="Temperature"
#    dim = 1


# In[ ]:


lJHTDB = libJHTDB()
lJHTDB.initialize()
lJHTDB.add_token(auth_token)

## "filename" parameter is the file names of output files, if filename='N/A', no files will be written. 
##             For example, if filename='results', the function will write "results.h5" and "results.xmf".
## The function only returns the data at the last time step within [t_start:t_step:t_end]
## The function only returns the data in the last field. For example, result=p if field=[up].
result = lJHTDB.getbigCutout(
        data_set=data_set, fields=fields, t_start=tstart, t_end=tend, t_step=tstep,
        start=np.array([xstart, ystart, zstart], dtype = np.int),
        end=np.array([xend, yend, zend], dtype = np.int),
        step=np.array([xstep, ystep, zstep], dtype = np.int),
        filter_width=Filter_Width,
        filename=output_filename)

lJHTDB.finalize()
print(result.shape)

