from pathlib import Path
import os,sys
R=Path.cwd();P=R;O=R/'experiments/E06_sampling_spatial_sensitivity';O.mkdir(exist_ok=True);sys.path.insert(0,str(P/'src'));os.environ['AQUATIC_PANEL']=os.environ.get('AQUATIC_PANEL',str(R/'data/raw/china_aquatic_buffering_panel_v2.parquet'))
import numpy as np,pandas as pd
import statsmodels.api as sm
from scipy.special import expit
from scipy.stats import t
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from audit_pipeline import load_panel
source=(P/'src/harmonized_confounding.py').read_text();source=source[source.index('_,d=load_panel()'):source.index('rows=[];frags=[];diags=[]')];source=source.replace("for name,s in [('M0',d),('M1',c)]:\n r,a=group_effects(s);r['model']=name;a['model']=name;save(r,out+name+'_fragility.csv');save(a,out+name+'_AFB.csv')",'')
inject="ctx=pd.read_csv(ROOT/'experiments/E01_spatial_support/national_corrected_1km_contexts.csv');d=d.merge(ctx,on='grid_id');d['context']=np.where(d.near_perm_1km_corrected==1,'permanent',np.where(d.near_seasonal_1km_corrected==1,'seasonal','other'));"
source=source.replace("fields=['elevation'",inject+"fields=['elevation'");ROOT=P;exec(source)
rows=[];mor=[];res=[]
for ev in ['EVI','LAI','combined']:
 fit=sm.GLM(s[ev].to_numpy(),x,family=sm.families.Binomial()).fit();b=fit.params;p=fit.fittedvalues;bread=np.linalg.inv((x*(p*(1-p))[:,None]).T@x);score=x*(s[ev].to_numpy()-p)[:,None];mus={};gr={}
 for ctx in ['permanent','seasonal','other']:
  for hh in [0,1]:
   xp=x.copy();xp[:,1]=hh;xp[:,2]=ctx=='seasonal';xp[:,3]=ctx=='other';xp[:,4]=hh*(ctx=='seasonal');xp[:,5]=hh*(ctx=='other');pp=expit(xp@b);mus[ctx,hh]=pp.mean();gr[ctx,hh]=(xp*(pp*(1-pp))[:,None]).mean(0)
 for width in [0,1,2,5]:
  labels=s.grid_id if width==0 else np.floor(s.lon/width).astype(str)+'_'+np.floor(s.lat/width).astype(str);codes,levels=pd.factorize(labels);G=len(levels);scores=np.zeros((G,x.shape[1]));np.add.at(scores,codes,score);V=bread@(scores.T@scores)@bread*(G/(G-1))*((len(s)-1)/(len(s)-x.shape[1]));crit=t.ppf(.975,G-1)
  for ctx in ['seasonal','other']:
   est=mus[ctx,1]-mus[ctx,0]-mus['permanent',1]+mus['permanent',0];grad=gr[ctx,1]-gr[ctx,0]-gr['permanent',1]+gr['permanent',0];se=np.sqrt(grad@V@grad);rows.append(dict(response=ev,comparison=ctx,block_degrees=width,n_clusters=G,AFB=est,se=se,ci_low=est-crit*se,ci_high=est+crit*se))
 z=s[['grid_id','lon','lat']].copy();z['residual']=s[ev].to_numpy()-p;g=z.groupby('grid_id').agg({'lon':'first','lat':'first','residual':'mean'});xy=np.column_stack([g.lon*np.cos(np.deg2rad(g.lat.mean())),g.lat])*111.32;idx=cKDTree(xy).query(xy,k=9)[1][:,1:];N=len(g);A=csr_matrix((np.ones(N*8),(np.repeat(np.arange(N),8),idx.ravel())),shape=(N,N));A=A.maximum(A.T);v=g.residual.to_numpy().copy();v-=v.mean();I=N/A.sum()*(v@(A@v))/(v@v);rng=np.random.default_rng(20261125);per=[]
 for k in range(199):
  vv=rng.permutation(v);per.append(N/A.sum()*(vv@(A@vv))/(vv@vv))
 mor.append(dict(response=ev,Morans_I=I,p=(1+np.sum(np.abs(np.array(per)+1/(N-1))>=abs(I+1/(N-1))))/200,n_grid=N,k=8));print(ev,'complete',flush=True)
pd.DataFrame(rows).to_csv(O/'adjusted_spatial_contrasts.csv',index=False);pd.DataFrame(mor).to_csv(O/'adjusted_residual_moran.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(mor).to_string(index=False))

