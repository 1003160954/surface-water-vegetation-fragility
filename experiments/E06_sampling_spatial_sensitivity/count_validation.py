from pathlib import Path
import sys,os,json
R=Path.cwd();P=R;O=R/'experiments/E06_sampling_spatial_sensitivity';sys.path.insert(0,str(P/'src'));os.environ['AQUATIC_PANEL']=os.environ.get('AQUATIC_PANEL',str(R/'data/raw/china_aquatic_buffering_panel_v2.parquet'))
from audit_pipeline import load_panel,group_effects
import numpy as np,pandas as pd
_,d=load_panel();ctx=pd.read_csv(P/'experiments/E01_spatial_support/national_corrected_1km_contexts.csv');d=d.merge(ctx,on='grid_id');d['context']=np.where(d.near_perm_1km_corrected==1,'permanent',np.where(d.near_seasonal_1km_corrected==1,'seasonal','other'));cv=pd.read_csv(P/'experiments/E02_confounding/harmonized_landcover_climate_covariates.csv');d=d.merge(cv,on='grid_id');q=d.HDSI.quantile([.25,.75]);vs=['EVI','LAI','VPD','LST_Day','SM_rootproxy','Precip'];source=Path(os.environ.get('RAW_CORE_PANEL','data/raw/china_panel_16day_core_patched.parquet'));p=pd.read_parquet(source,columns=['grid_id','year','step16','doy_mid']+vs);p.loc[~p.EVI.between(-.2,1.2),'EVI']=np.nan;p.loc[~p.LAI.between(0,12),'LAI']=np.nan
counts=[];merged=None
for v in vs:
 g=p.groupby(['grid_id','step16'])[v];n=g.transform('count');mu=g.transform('mean');sd=g.transform('std');xx=p[v];n0=n-1;mu0=(n*mu-xx)/n0;ss=(n-1)*sd**2-(xx-mu)*(xx-mu0);sd0=np.sqrt(np.maximum(ss,0)/(n0-1));valid=xx.notna()&(sd>1e-6);held=xx.notna()&(sd0>1e-6)&(n0>=3);mask=p.doy_mid.between(*( (213,273) if v in ['EVI','LAI'] else (152,243)))
 z=p.loc[mask,['grid_id','year']].copy();z['nstep']=valid[mask].astype(int);z['held']=held[mask].astype(int);z['ref']=n[mask].where(valid[mask]);a=z.groupby(['grid_id','year']).agg(nstep=('nstep','sum'),held=('held','sum'),min_ref=('ref','min')).reset_index();a=a.rename(columns={c:v+'_'+c for c in ['nstep','held','min_ref']});merged=a if merged is None else merged.merge(a,on=['grid_id','year']);print(v,flush=True)
s=d.merge(merged,on=['grid_id','year']);s.to_parquet(O/'count_validation_panel.parquet',index=False)
for v in vs:
 for ctx,z in s.groupby('context'):
  for field in ['nstep','min_ref']:
   for value,n in z[v+'_'+field].value_counts(dropna=False).items():counts.append(dict(variable=v,context=ctx,quantity=field,value=value,n_records=n))
pd.DataFrame(counts).to_csv(O/'observation_counts_by_context.csv',index=False)
rows=[]
for minimum in [2,3]:
 z=s[s[[v+'_nstep' for v in vs]].ge(minimum).all(axis=1)];r,a=group_effects(z,B=2000,q_fixed=q);a['minimum_steps']=minimum;rows.append(a);print('min',minimum,len(z),flush=True)
pd.concat(rows).to_csv(O/'minimum_steps_AFB.csv',index=False)
s['retained']=s[[v+'_held' for v in vs]].ge(1).all(axis=1);assert s.retained.sum()==209985
s['region']='Other';conds=[(s.elevation.ge(3000)&s.lon.between(75,105)&s.lat.between(25,42),'Tibetan Plateau'),(s.lon.ge(118)&s.lat.ge(40),'Northeast'),(s.lon.lt(105)&s.lat.ge(35),'Northwest arid'),(s.lon.ge(105)&s.lon.lt(122)&s.lat.ge(34)&s.lat.lt(40),'North China'),(s.lon.ge(110)&s.lon.lt(123)&s.lat.ge(27)&s.lat.lt(34),'East China'),(s.lon.ge(105)&s.lat.lt(27),'South China'),(s.lon.lt(110)&s.lat.lt(35),'Southwest')];assigned=np.zeros(len(s),bool)
for mask,label in conds:s.loc[mask&~assigned,'region']=label;assigned|=mask
rows=[]
for field in ['context','LC2018_grid_mode','region']:
 for cat,z in s.groupby(field):rows.append(dict(grouping=field,category=cat,eligible=len(z),retained=int(z.retained.sum()),excluded=int((~z.retained).sum()),retention_fraction=z.retained.mean()))
pd.DataFrame(rows).to_csv(O/'heldout_exclusions.csv',index=False);print('ALL DONE',flush=True)
