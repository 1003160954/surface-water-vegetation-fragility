from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime,sys,os
ROOT=Path(__file__).resolve().parents[1]; ROOT.mkdir(parents=True,exist_ok=True)
for p in ['gee','src','figures','supplementary','tests','configs','docs','logs','experiments/E00_original_reproduction','experiments/E01_spatial_support','experiments/E02_confounding','experiments/E03_independent_baseline','experiments/E04_robustness','experiments/E05_figure_audit']:(ROOT/p).mkdir(parents=True,exist_ok=True)
DATA=Path(os.environ.get('AQUATIC_PANEL',str(ROOT/'data/raw/china_aquatic_buffering_panel_v2.parquet'))); SEED=20261020

def save(df,path):
 p=ROOT/path;p.parent.mkdir(parents=True,exist_ok=True);df.to_csv(p,index=False)

def load_panel():
 raw=pd.read_parquet(DATA); d=raw.loc[raw.valid_model_both_y].copy()
 d['context']=d.water_context.replace({'near_permanent_water':'permanent','near_seasonal_water':'seasonal','far_or_no_water':'other'})
 d['EVI']=(d.late_summer_EVI_z<-.5).astype(int);d['LAI']=(d.late_summer_LAI_z<-.5).astype(int);d['combined']=((d.EVI==1)|(d.LAI==1)).astype(int)
 return raw,d

def group_effects(d,group='context',stress='HDSI',B=300,cluster='grid_id',seed=SEED,q_fixed=None):
 # Quantiles are fixed at original full analysis-sample values within each specified experiment.
 q=d[stress].quantile([.25,.75]) if q_fixed is None else q_fixed; d=d.copy(); d['hi']=d[stress]>q.iloc[1]; d['lo']=d[stress]<q.iloc[0]
 cats=list(d[group].dropna().unique()); rows=[];draws={};rng=np.random.default_rng(seed)
 clusters=pd.Index(d[cluster].drop_duplicates());ci=clusters.get_indexer(d[cluster]);G=len(clusters)
 # Common resampling weights preserve covariance between exposure groups and responses.
 weights=np.array([np.bincount(rng.integers(G,size=G),minlength=G) for _ in range(B)],dtype=np.float32) if B else None
 for c in cats:
  m=d[group].eq(c).to_numpy(); sums=[]
  for event in ['EVI','LAI','combined']:
   mats=np.column_stack([np.bincount(ci,weights=m*d['lo'].to_numpy(),minlength=G),np.bincount(ci,weights=m*d['hi'].to_numpy(),minlength=G),np.bincount(ci,weights=m*d['lo'].to_numpy()*d[event].to_numpy(),minlength=G),np.bincount(ci,weights=m*d['hi'].to_numpy()*d[event].to_numpy(),minlength=G)])
   nlo,nhi,slo,shi=mats.sum(axis=0); effect=shi/nhi-slo/nlo if nlo*nhi else np.nan
   boot=weights@mats if B else np.zeros((0,4)); vals=boot[:,3]/np.where(boot[:,1]>0,boot[:,1],np.nan)-boot[:,2]/np.where(boot[:,0]>0,boot[:,0],np.nan)
   draws[c,event]=vals
   lo,hi=np.nanquantile(vals,[.025,.975]) if np.isfinite(vals).any() else (np.nan,np.nan)
   rows.append(dict(group=c,response=event,n_total=int(m.sum()),n_grid=d.loc[m,'grid_id'].nunique(),n_low=int(nlo),n_high=int(nhi),risk_low=slo/nlo if nlo else np.nan,risk_high=shi/nhi if nhi else np.nan,fragility=effect,ci_low=lo,ci_high=hi,q25=q.iloc[0],q75=q.iloc[1],cluster=cluster,n_clusters=G,bootstrap=B,seed=seed,valid_draws=int(np.isfinite(vals).sum())))
 out=pd.DataFrame(rows);contrasts=[]
 if 'permanent' in cats:
  for c in cats:
   if c=='permanent':continue
   for e in ['EVI','LAI','combined']:
    a=out[(out.group==c)&(out.response==e)].iloc[0];b=out[(out.group=='permanent')&(out.response==e)].iloc[0];v=draws[c,e]-draws['permanent',e];lo,hi=np.nanquantile(v,[.025,.975]) if np.isfinite(v).any() else (np.nan,np.nan)
    contrasts.append(dict(comparison=c,response=e,AFB=a.fragility-b.fragility,ci_low=lo,ci_high=hi,n_total=len(d),n_grid=d.grid_id.nunique(),n_low=int(d.lo.sum()),n_high=int(d.hi.sum()),cluster=cluster,n_clusters=G,bootstrap=B,seed=seed,valid_draws=int(np.isfinite(v).sum())))
 return out,pd.DataFrame(contrasts)

def main():
 raw,d=load_panel();o='experiments/E00_original_reproduction/'
 save(pd.DataFrame({'field':raw.columns,'dtype':raw.dtypes.astype(str).values,'missing_n':raw.isna().sum().values,'missing_fraction':raw.isna().mean().values}),o+'data_dictionary.csv')
 stages=[('raw',raw),('valid_x',raw[raw.valid_x]),('valid_any_response',raw[raw.valid_model_any_y]),('valid_both_responses',d)]
 save(pd.DataFrame([dict(stage=k,n_rows=len(v),n_grids=v.grid_id.nunique(),years=','.join(map(str,sorted(v.year.unique())))) for k,v in stages]),o+'sample_flow.csv')
 vals={'rows':len(d),'unique_grids':d.grid_id.nunique(),'HDSI_Q25':d.HDSI.quantile(.25),'HDSI_Q75':d.HDSI.quantile(.75),'EVI_neg_rate':d.EVI.mean(),'LAI_neg_rate':d.LAI.mean(),'joint_failure_rate':d.combined.mean(),'dist_perm_mean_km':d.dist_perm_water_km.mean(),'dist_perm_median_km':d.dist_perm_water_km.median(),'dist_perm_p90_km':d.dist_perm_water_km.quantile(.9),'dist_perm_max_km':d.dist_perm_water_km.max()}
 hist=pd.read_csv(ROOT/'docs/historical/main_overview.csv').set_index('metric').value
 save(pd.DataFrame([dict(metric=k,recomputed=v,historical=float(hist[k]),difference=v-float(hist[k])) for k,v in vals.items()]),o+'overview_comparison.csv')
 save(d.groupby('year').agg(n=('grid_id','size'),n_grid=('grid_id','nunique'),EVI_event=('EVI','mean'),LAI_event=('LAI','mean'),combined_event=('combined','mean')).reset_index(),o+'year_summary.csv')
 save(d.groupby('region_id').agg(n=('grid_id','size'),n_grid=('grid_id','nunique')).reset_index(),o+'source_tile_counts.csv')
 save(d.groupby('water_dist_band',observed=True).dist_perm_water_km.agg(['min','max','mean','count']).reset_index(),o+'stored_band_audit.csv')
 flags={'duplicate_grid_year':int(raw.duplicated(['grid_id','year']).sum()),'fraction_noninteger_landcover':float((abs(d.LC_Type1-d.LC_Type1.round())>1e-6).mean()),'HDSI_max_abs_difference':float(abs(d.HDSI-(d.summer_VPD_z+d.summer_LST_Day_z-d.summer_SM_rootproxy_z-d.summer_Precip_z)/4).max()),'complete_case_covariates':int(d[['elevation','slope','LC_Type1','lon','lat','year']].notna().all(axis=1).sum())}
 (ROOT/o/'audit_flags.json').write_text(json.dumps(flags,indent=2))
 risk,afb=group_effects(d);save(risk,o+'context_fragility.csv');save(afb,o+'AFB.csv')
 h=pd.read_csv(ROOT/'docs/historical/aquatic_fragility_buffering_AFB_with_CI.csv');h['response']=h.response.replace({'Joint':'combined'});h['comparison']=h.contrast.str.startswith('near_seasonal').map({True:'seasonal',False:'other'});cmp=afb.merge(h[['response','comparison','AFB']],on=['response','comparison'],suffixes=('_recomputed','_historical'));cmp['difference']=cmp.AFB_recomputed-cmp.AFB_historical;save(cmp,o+'AFB_point_comparison.csv')
 d['distance_band']=pd.cut(d.dist_perm_water_km,[-1e-9,.5,1,2,5,10,np.inf],labels=['0-0.5','0.5-1','1-2','2-5','5-10','>10']).astype(str)
 band,_=group_effects(d,'distance_band');save(band,o+'distance_fragility.csv')
 d['fine_band']=pd.cut(d.dist_perm_water_km,[-1e-9,.25,.5,1,2,3,5,7.5,10,15,np.inf]).astype(str);fine,_=group_effects(d,'fine_band');save(fine,o+'fine_distance_fragility.csv')
 manifest={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input':str(DATA),'sha256':hashlib.file_digest(DATA.open('rb'),'sha256').hexdigest(),'bytes':DATA.stat().st_size,'rows_raw':len(raw),'rows_analysis':len(d),'seed':SEED,'bootstrap':300,'original_notebook_found':(ROOT/'docs/original_sources/cell_59.py').exists(),'method':'reconstruction; recovered original code validated separately','python':sys.version}
 (ROOT/'logs/E00_manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(vals,indent=2));print(afb.to_string(index=False));print('E00 FINISHED',flush=True)
if __name__=='__main__':main()
