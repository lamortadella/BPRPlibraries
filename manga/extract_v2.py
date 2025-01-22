from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import pandas as pd
import pathlib
import numpy as np

#read catalog of SDSS-Gaia data points provided by Michele
pdir = pathlib.Path('.').resolve()
ff1 = pdir / 'mastarall-v3_1_1-v1_7_7_goodstars_inDR3.fits' #22144 entries
df1 = Table.read(ff1)
print(df1.columns)
df1=df1['MANGAID','source_id','RAdeg','DEdeg','INPUT_LOGG','INPUT_TEFF','INPUT_FE_H','INPUT_ALPHA_M']
df1=df1['MANGAID','source_id']
df1=df1.to_pandas()
df1['source_id']=df1['source_id'].astype(str)
df1['sMANGAID']=df1['MANGAID'].astype(str)

#read the calspec library
dir="./"
calsp=pd.read_csv(dir+'CALSPEC-resultCOORDn.csv')
calsp['source_id']=calsp['source_id'].astype(str)
print(calsp.columns)

#MERGE
df2=df1.merge(calsp, how='inner', on=['source_id','source_id'])
df2['MANGAID']=df2.MANGAID.str.decode('utf-8') 
print(df2)

# READ SPECTRA =read in the HDUs using a FITS reader (astropy, pyfits, or fitsio)
hdulist = fits.open('mastar-goodspec-v3_1_1-v1_7_7.fits')
hdulist.info()
lala = hdulist[1].header
data = hdulist[1].data

#create pandas array of spectral parameters
#print(data.columns)
lala=data['MANGAID'].tolist()
lala2=np.array(lala)
spID = pd.DataFrame(lala2.T, columns=['1'])
spID2=spID.rename(columns={'1':'MANGAID'})

lala=data['PLATE'].tolist()
lala2=np.array(lala)
plate=lala2.T
spID2['PLATE']=plate

lala=data['MJD'].tolist()
lala2=np.array(lala)
mjd=lala2.T
spID2['MJD']=mjd

lala=data['EXPTIME'].tolist()
lala2=np.array(lala)
mjd=lala2.T
spID2['EXPTIME']=mjd

lala=data['NEXP_VISIT'].tolist()
lala2=np.array(lala)
mjd=lala2.T
spID2['NEXP_VISIT']=mjd


#add spectral parameters into the table
df3=df2.merge(spID2, how='left', on=['MANGAID','MANGAID'])

#save the table
df3.to_csv('calspec_sdss_sp.csv',index=0)
print(df3)

crossid=df3['MANGAID'].astype('str')

lala=data['MANGAID'].tolist()
lala2=np.array(lala)
spID = pd.DataFrame(lala2.T, columns=['1'])
datamanga=spID.rename(columns={'1':'MANGAID'})

nlen=len(crossid)
for i in range(0,nlen):
  print(crossid.values[i])
  idx=datamanga.loc[datamanga['MANGAID']==crossid.values[i]].index
  print(crossid.values[i])
  
  print(idx,idx.size)
  np=idx.size
  if(np > 0):
   for j in range(0,np):
    plt.plot(data['WAVE'][idx[j]],data['FLUX'][idx[j]])
    #plt.show()
    plt.savefig(crossid.values[i]+"_"+str(j)+'.eps',dpi=300,bbox_inches="tight")
  
    print(j, crossid.values[i])
    wave=data['WAVE'][idx[j]]
    spectrum = pd.DataFrame(wave.T, columns=['1'])
    flux=data['FLUX'][idx[j]]
    spectrum['flux']=flux
    spectrum2=spectrum.rename(columns={'1':'wave'})
    spectrum2.to_csv(crossid.values[i]+"_"+str(j)+'.csv',index=0)
    
 
    
