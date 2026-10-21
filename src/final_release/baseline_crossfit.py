from pathlib import Path
import sys,json,gc
import numpy as np,pandas as pd
import os
ROOT=Path(__file__).resolve().parents[2];os.chdir(ROOT)
sys.path.insert(0,str(ROOT/'src'))
from audit_pipeline import load_panel,group_effects
out=Path('final_results/data/baseline');out.mkdir(parents=True,exist_ok=True)
source=Path(os.environ.get('RAW_CORE_PANEL','data/raw/china_panel_16day_core_patched.parquet'))
vs=['EVI','LAI','VPD','LST_Day','SM_rootproxy','Precip']
p=pd.read_parquet(source,columns=['grid_id','year','step16','doy_mid']+vs)
assert not p.duplicated(['grid_id','year','step16']).any()
p.loc[(p.EVI<-.2)|(p.EVI>1.2),'EVI']=np.nan;p.loc[(p.LAI<0)|(p.LAI>12),'LAI']=np.nan
_,d=load_panel();ctx=pd.read_csv('experiments/E01_spatial_support/national_corrected_1km_contexts.csv')
d=d.merge(ctx,on='grid_id',validate='many_to_one');d['context']=np.where(d.near_perm_1km_corrected==1,'permanent',np.where(d.near_seasonal_1km_corrected==1,'seasonal','other'))
merged=None;countrows=[]
for v in vs:
 g=p.groupby(['grid_id','step16'])[v];n=g.transform('count').to_numpy();mu=g.transform('mean').to_numpy();sd=g.transform('std').to_numpy();x=p[v].to_numpy();finite=np.isfinite(x)
 # Delete the target year's single value from the sufficient statistics.
 n0=n-1;mu0=(n*mu-x)/np.where(n0>0,n0,np.nan)
 ss=(n-1)*sd**2-(x-mu)*(x-mu0)
 sd0=np.sqrt(np.maximum(ss,0)/np.where(n0>1,n0-1,np.nan))
 z=np.clip((x-mu)/np.where(sd>1e-6,sd,np.nan),-5,5)
 z0=np.clip((x-mu0)/np.where((sd0>1e-6)&(n0>=3),sd0,np.nan),-5,5)
 season=(213,273) if v in ['EVI','LAI'] else (152,243);mask=p.doy_mid.between(*season).to_numpy()
 s=p.loc[mask,['grid_id','year']].copy();s[v+'_reproduced']=z[mask];s[v+'_heldout']=z0[mask];s[v+'_count']=np.isfinite(z0[mask]).astype(int)
 a=s.groupby(['grid_id','year']).agg({v+'_reproduced':'mean',v+'_heldout':'mean',v+'_count':'sum'}).reset_index()
 merged=a if merged is None else merged.merge(a,on=['grid_id','year'],how='outer',validate='one_to_one')
 for c in range(7):countrows.append(dict(variable=v,reference_count=c,n_target_steps=int((mask&finite&(n==c)).sum())))
 print('finished',v,flush=True)
del p;gc.collect()
s=d.merge(merged,on=['grid_id','year'],validate='one_to_one');checks=[]
for v in vs:
 col=('late_summer_' if v in ['EVI','LAI'] else 'summer_')+v+'_z'
 err=abs(s[col]-s[v+'_reproduced']);checks.append(dict(variable=v,max_abs_difference=err.max(),n_compared=int(err.notna().sum())))
 assert err.max()<1e-8,(v,err.max())
pd.DataFrame(checks).to_csv(out/'reference_reproduction.csv',index=False);pd.DataFrame(countrows).to_csv(out/'reference_sample_counts.csv',index=False)
ok=s[[v+'_heldout' for v in vs]].notna().all(axis=1);z=s[ok].copy();q=d.HDSI.quantile([.25,.75])
results=[];agreement=[]
for ev in ['EVI','LAI']:
 old=z[ev].to_numpy();new=(z[ev+'_heldout']<-.5).astype(int).to_numpy()
 agreement.append(dict(response=ev,n=len(z),agreement=float((old==new).mean()),original_event_rate=old.mean(),heldout_event_rate=new.mean(),changed=int((old!=new).sum())))
newcombined=((z.EVI_heldout<-.5)|(z.LAI_heldout<-.5)).astype(int)
agreement.append(dict(response='combined',n=len(z),agreement=float((z.combined==newcombined).mean()),original_event_rate=z.combined.mean(),heldout_event_rate=newcombined.mean(),changed=int((z.combined!=newcombined).sum())))
for name in ['reference_same_rows','heldout_response_fixed_exposure','heldout_response_and_exposure']:
 t=z.copy()
 if name!='reference_same_rows':
  t['EVI']=(t.EVI_heldout<-.5).astype(int);t['LAI']=(t.LAI_heldout<-.5).astype(int);t['combined']=t[['EVI','LAI']].max(axis=1)
 if name=='heldout_response_and_exposure':t['HDSI']=(t.VPD_heldout+t.LST_Day_heldout-t.SM_rootproxy_heldout-t.Precip_heldout)/4
 r,a=group_effects(t,q_fixed=None if name=='heldout_response_and_exposure' else q);a['scenario']=name;r['scenario']=name;results.append(a);r.to_csv(out/(name+'_risks.csv'),index=False)
pd.concat(results).to_csv(out/'baseline_AFB.csv',index=False);pd.DataFrame(agreement).to_csv(out/'event_agreement.csv',index=False)
cfg=dict(source=str(source),source_rows=5619038,reference='exclude target year at each grid/16-day position',min_reference_observations=3,sample_sd_ddof=1,sd_mask=1e-6,clip=[-5,5],eligible_main_rows=len(d),paired_valid_rows=len(z),paired_valid_grids=z.grid_id.nunique(),excluded_rows=int((~ok).sum()),bootstrap=300,seed=20261020,long_independent_baseline=False)
(out/'config.json').write_text(json.dumps(cfg,indent=2));print(json.dumps(cfg));print(pd.DataFrame(agreement).to_string(index=False));print(pd.concat(results).to_string(index=False))
