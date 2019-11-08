#import packages
from numpy import *
from scipy.misc import *
from geopy.distance import vincenty
import matplotlib.pyplot as plt
from IPython import get_ipython
import pickle
from netCDF4 import Dataset
import time
from joblib import Parallel, delayed
import multiprocessing



tic = time.time()

#clear figures, console and variables
plt.close("all")
print(chr(27) + "[2J")
#get_ipython().magic('reset -sf')
#from matplotlib import cm


# Define paths and filenames

# Load data
path = '/Users/JennaPalmer/Google Drive/Baylor & Jenna/Gulf of Maine/Data/'
source = Dataset(path+'GLAD.nc','r')


u = source.variables['u'][:]
u = transpose(u)
v = source.variables['v'][:]
v = transpose(v)
lat = source.variables['lat'][:]
lat = transpose(lat)
lon = source.variables['lon'][:]
lon = transpose(lon)

source.close()

[nd,nt] = u.shape

oceantime = arange(nt)


# initialize variables
combos = int(comb(nd,2))

# output variables
l = zeros((combos,nt))
t = zeros((combos,nt))
r = zeros((combos,nt))


#for tt in range(0,nt): 
def timeloop(tt):   
    # temp variables
#    ds = zeros((nd-1,nd))
#    du2 = zeros((nd-1,nd))
#    dv2 = zeros((nd-1,nd))

    d = 0;
    for ii in range(0,nd-1):
        for jj in range(ii+1,nd):           
            du = u[jj,tt]-u[ii,tt]
            dv = v[jj,tt]-v[ii,tt]
            
            if isnan(du) == False:
                dx = vincenty((lat[ii,tt],lon[ii,tt]),(lat[ii,tt],lon[jj,tt])).km
                dy = vincenty((lat[ii,tt],lon[ii,tt]),(lat[jj,tt],lon[ii,tt])).km
                r[d,tt] = vincenty((lat[ii,tt],lon[ii,tt]),(lat[jj,tt],lon[jj,tt])).km
            
                # get correct sign as vincenty returns absolute distance
                dx = sign(lon[jj,tt]-lon[ii,tt])*dx
                dy = sign(lat[jj,tt]-lat[ii,tt])*dy
            
            
                dl = (du*dx + dv*dy)/r[d,tt]
                dt = (-du*dy + dv*dx)/r[d,tt]
            
                l[d,tt] = power(dl,2)
                t[d,tt] = power(dt,2)
                
            else:
                l[d,tt] = nan
                t[d,tt] = nan
                r[d,tt] = nan
            d = d+1
            
    
    prev = int(round(float(tt)/nt,2)*100)
    current = int(round(float(tt+1)/nt,2)*100)
    if current != prev:
        print('Progress Calculating: ' + str(current) + '%')
        
    return l,t,r
        
# Find number of workers
num_cores = multiprocessing.cpu_count()

#Calculate structure functions 
results = Parallel(n_jobs=num_cores)(delayed(timeloop)(tt) for tt in range(0,nt))


############## Binning ########################

# declare vars
Bins = []

kk=0
Bins.extend([0])
while power(10,((kk+1)*0.25)) <= ceil(nanmax(r)):
        Bins.extend([power(10,((kk+1)*0.25))])
        kk+=1
        
lsf = zeros((len(Bins)-1,nt))
tsf = zeros((len(Bins)-1,nt))
count = zeros((len(Bins)-1,nt));

for tt in range(0,nt): 
    for ii in range(0,len(Bins)-1):
        ltemp = l[:,tt]
        ttemp = t[:,tt]
        
        ind = where(logical_and(r[:,tt] <= Bins[ii+1],r[:,tt] >Bins[ii]))
        
        lsf[ii,tt] = nanmean(ltemp[ind])
        tsf[ii,tt] = nanmean(ttemp[ind])
        count[ii,tt] = len(ltemp[~isnan(ltemp)])
        
        
    prev = int(round(float(tt)/nt,2)*100)
    current = int(round(float(tt+1)/nt,2)*100)
    if current != prev:
        print('Progress Binning: ' + str(current) + '%')

# Get ensemble SFs

sumtsf = nansum(multiply(tsf,count), axis=1)
tottsf = nansum(count,axis = 1)
enstsf = divide(sumtsf,tottsf)

sumlsf = nansum(multiply(lsf,count), axis=1)
totlsf = nansum(count,axis = 1)
enslsf = divide(sumlsf,totlsf)

         
# center bins
rmid = []
#rmid.extend([.5])

for bb in range(0,len(Bins)-1):
     rmid.extend([(Bins[bb]+Bins[bb+1])/2])
     


# save variables to dictionary
VSF = {'l_unbinned': l ,'t_unbinned':t,'r_unbinned':r, 'l_binned': lsf, 't_binned': tsf, 'r_binned': rmid, 'freq':count, 'enstsf': enstsf, 'enslsf':enslsf}
# save dictionary to file
pickle.dump( VSF, open( path+"GLAD_VSF.p", "wb" ) )
    
#data = pickle.load( open( path+"GLAD_VSF.p", "rb" ) )

toc = time.time()-tic
print('Seconds elapsed: ' + str(toc))

loglog(rmid,enstsf)


