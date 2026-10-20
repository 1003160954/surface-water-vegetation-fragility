from pathlib import Path
import sys,ast,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
_,d=load_panel();q=d.HDSI.quantile([.25,.75]);d['exposure']=np.where(d.HDSI<=q.iloc[0],'low',np.where(d.HDSI>=q.iloc[1],'high','mid'))
src=(ROOT/'docs/original_sources/cell_59.py').read_text(encoding='utf8');node=next(n for n in ast.parse(src).body if isinstance(n,ast.FunctionDef) and n.name=='assign_macro_region');exec(compile(ast.Module(body=[node],type_ignores=[]),'<recovered-cell59-function>','exec'))
d['macro_region']=d.apply(assign_macro_region,axis=1);o='experiments/E00_original_reproduction/';rows=[]
for e in ['EVI','LAI','combined']:
 arrays={}
 for ctx in ['permanent','other','seasonal']:
  for lev in ['low','high']:arrays[ctx,lev]=d.loc[(d.context==ctx)&(d.exposure==lev),e].to_numpy()
 rng=np.random.default_rng(42);boots={c:[] for c in ['other','seasonal']}
 for _ in range(500):
  a=[rng.choice(arrays[c,l],size=len(arrays[c,l]),replace=True).mean() for c in ['permanent','other','seasonal'] for l in ['low','high']]
  boots['other'].append(a[3]-a[2]-a[1]+a[0]);boots['seasonal'].append(a[5]-a[4]-a[1]+a[0])
 for c,v in boots.items():
  lo,hi=np.percentile(v,[2.5,97.5]);rows.append(dict(response=e,comparison=c,AFB=np.mean(v),ci_low=lo,ci_high=hi,n_total=len(d),n_grid=d.grid_id.nunique(),bootstrap=500,seed=42,unit='independent grid-year within context/exposure',point_estimator='mean of bootstrap contrasts'))
a=pd.DataFrame(rows);hist=pd.read_csv(ROOT/'docs/historical/aquatic_fragility_buffering_AFB_with_CI.csv');hist['response']=hist.response.replace({'Joint':'combined'});hist['comparison']=hist.contrast.str.startswith('near_seasonal').map({True:'seasonal',False:'other'});c=a.merge(hist[['response','comparison','AFB','ci_low','ci_high']],on=['response','comparison'],suffixes=('_rerun','_historical'))
for k in ['AFB','ci_low','ci_high']:c[k+'_difference']=c[k+'_rerun']-c[k+'_historical']
save(a,o+'original_cell59_AFB_rerun.csv');save(c,o+'original_cell59_exact_comparison.csv');assert c[['AFB_difference','ci_low_difference','ci_high_difference']].abs().to_numpy().max()<1e-12
# Reconstruct all region/context/exposure supports; do not omit failed estimands.
save(d.groupby(['macro_region','context','exposure']).agg(n_rows=('grid_id','size'),n_grid=('grid_id','nunique')).reset_index(),o+'regional_context_exposure_counts.csv');save(d[['grid_id','year','macro_region','context','exposure']],o+'sample_region_map.csv')
regions=[]
for region,s in d.groupby('macro_region'):
 near=s[s.context=='permanent'];other=s[s.context=='other'];counts={c+'_'+lev:int(((s.context==c)&(s.exposure==lev)).sum()) for c in ['permanent','other'] for lev in ['low','high']};reason=''
 if len(near)<100 or len(other)<100:reason='original total context count below 100'
 elif min(counts.values())<30:reason='original low/high cell count below 30'
 for e in ['EVI','LAI','combined']:
  row=dict(region=region,response=e,n_near=len(near),n_other=len(other),n_grids_region=s.grid_id.nunique(),**counts,omitted_original=bool(reason),reason=reason)
  if not reason:
   arr=[s.loc[(s.context==c)&(s.exposure==l),e].to_numpy() for c in ['permanent','other'] for l in ['low','high']];rng=np.random.default_rng(42);vals=[]
   for _ in range(400):
    v=[rng.choice(x,size=len(x),replace=True).mean() for x in arr];vals.append(v[3]-v[2]-v[1]+v[0])
   row.update(AFB=np.mean(vals),ci_low=np.percentile(vals,2.5),ci_high=np.percentile(vals,97.5),bootstrap=400,seed=42)
  regions.append(row)
save(pd.DataFrame(regions),o+'original_regional_AFB_rerun_with_omissions.csv');print(c.to_string(index=False));print(pd.DataFrame(regions).query("response == 'combined'").to_string(index=False));print('EXACT ORIGINAL AFB REPRODUCED',flush=True)
