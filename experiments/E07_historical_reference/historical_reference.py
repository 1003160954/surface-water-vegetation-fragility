from pathlib import Path
import os,sys,json,gc,time
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[2];O=Path(os.environ.get('HISTORICAL_OUTPUT',str(Path(__file__).resolve().parent)));O.mkdir(parents=True,exist_ok=True);sys.path.insert(0,str(P/'src'));from audit_pipeline import load_panel,group_effects
_,d=load_panel();ctx=pd.read_csv(P/'experiments/E01_spatial_support/national_corrected_1km_contexts.csv');d=d.merge(ctx,on='grid_id');d['context']=np.where(d.near_perm_1km_corrected==1,'permanent',np.where(d.near_seasonal_1km_corrected==1,'seasonal','other'));grids=set(d.grid_id);vs=['EVI','LAI','VPD','LST_Day','SM_rootproxy','Precip'];root=Path(os.environ['HISTORICAL_ARCHIVE']);regions=['R_5_2','R_5_3','R_6_1','R_6_2','R_6_3','R_6_4'];manifest=[]
def clean(z):
 for v in vs:z[v]=pd.to_numeric(z[v],errors='coerce').mask(lambda x:x<=-9990)
 z.loc[~z.EVI.between(-.2,1.2),'EVI']=np.nan;z.loc[~z.LAI.between(0,12),'LAI']=np.nan
 return z
def build_region(reg):
 manifest=[]
 target=O/(reg+'_reference.parquet')
 if target.exists():print('cached',reg,flush=True);return
 acc=None
 for year in range(2001,2018):
  src=root/f'year={year}'/f'region_id={reg}'/f'dynamic_step16_y{year}_{reg}.parquet';z=pd.read_parquet(src,columns=['grid_id','year','step16']+vs);z=z[z.grid_id.isin(grids)].copy();assert not z.duplicated(['grid_id','year','step16']).any();z=clean(z);z=z.set_index(['grid_id','step16']);v=z[vs];a=pd.concat([v.fillna(0).add_suffix('_sum'),(v.fillna(0)**2).add_suffix('_ss'),v.notna().astype(int).add_suffix('_n')],axis=1);acc=a if acc is None else acc.add(a,fill_value=0);manifest.append({'file':str(src.relative_to(root)).replace('\\','/'),'year':year,'region':reg,'matched_rows':len(z),'bytes':src.stat().st_size});print(reg,year,len(z),flush=True)
 for v in vs:
  n=acc[v+'_n'];mu=acc[v+'_sum']/n;var=(acc[v+'_ss']-n*mu**2)/(n-1);acc[v+'_mean']=mu;acc[v+'_sd']=np.sqrt(var.clip(lower=0))
 acc.reset_index().to_parquet(target,index=False);pd.DataFrame(manifest).to_csv(O/(reg+'_input_inventory.csv'),index=False);del acc;gc.collect()

from concurrent.futures import ThreadPoolExecutor
with ThreadPoolExecutor(max_workers=3) as pool:
 list(pool.map(build_region,regions))
ref=pd.concat([pd.read_parquet(O/(r+'_reference.parquet')) for r in regions],ignore_index=True);assert not ref.duplicated(['grid_id','step16']).any()
source=Path(os.environ['RAW_CORE_PANEL']);p=pd.read_parquet(source,columns=['grid_id','year','step16','doy_mid']+vs);p=clean(p[p.grid_id.isin(grids)].copy());p=p.merge(ref,on=['grid_id','step16'],validate='many_to_one');summary=[];results=[];agreements=[];supports=[]
for minimum in [10,15]:
 merged=None
 for v in vs:
  season=p.doy_mid.between(*( (213,273) if v in ['EVI','LAI'] else (152,243)));valid=(p[v+'_n']>=minimum)&(p[v+'_sd']>1e-6);zz=((p[v]-p[v+'_mean'])/p[v+'_sd']).where(valid).clip(-5,5);z=p.loc[season,['grid_id','year']].copy();z[v]=zz[season];z[v+'_count']=z[v].notna().astype(int);a=z.groupby(['grid_id','year']).agg({v:'mean',v+'_count':'sum'}).reset_index();merged=a if merged is None else merged.merge(a,on=['grid_id','year']);summary.append({'min_reference_years':minimum,'variable':v,'valid_target_steps':int(z[v].notna().sum())})
 s=d.merge(merged,on=['grid_id','year'],how='left',suffixes=('','_long'));ok=s[[('EVI_long' if v=='EVI' else 'LAI_long' if v=='LAI' else v) for v in vs]].notna().all(axis=1);z=s[ok].copy();supports.append({'min_reference_years':minimum,'n_rows':len(z),'n_grids':z.grid_id.nunique(),'excluded':len(d)-len(z)});s[['grid_id','year','context']].assign(retained=ok).to_csv(O/f'sample_retention_min{minimum}.csv',index=False)
 for name in ['main_reference_paired','independent_response_fixed_exposure','independent_response_and_exposure']:
  t=z.copy()
  if name!='main_reference_paired':t['EVI']=(t.EVI_long<-.5).astype(int);t['LAI']=(t.LAI_long<-.5).astype(int);t['combined']=t[['EVI','LAI']].max(axis=1)
  if name=='independent_response_and_exposure':t['HDSI']=(t.VPD+t.LST_Day-t.SM_rootproxy-t.Precip)/4
  r,a=group_effects(t,B=2000,q_fixed=None if name=='independent_response_and_exposure' else d.HDSI.quantile([.25,.75]));r.assign(scenario=name,min_reference_years=minimum).to_csv(O/f'context_risks_min{minimum}_{name}.csv',index=False);a['scenario']=name;a['min_reference_years']=minimum;results.append(a)
  for ev in ['EVI','LAI','combined']:agreements.append({'minimum':minimum,'scenario':name,'event':ev,'event_rate':t[ev].mean(),'agreement_with_main':(t[ev].to_numpy()==z[ev].to_numpy()).mean()})
  print('RESULT',minimum,name,a[['comparison','response','AFB','ci_low','ci_high']].to_dict('records'),flush=True)
 pd.concat(results).to_csv(O/'independent_reference_AFB.csv',index=False);pd.DataFrame(agreements).to_csv(O/'event_agreement.csv',index=False);pd.DataFrame(supports).to_csv(O/'sample_support.csv',index=False)
pd.DataFrame(summary).to_csv(O/'valid_steps.csv',index=False);print('COMPLETE',flush=True)
