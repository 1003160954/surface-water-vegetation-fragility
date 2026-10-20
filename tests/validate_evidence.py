from pathlib import Path
import sys,json,datetime,ast
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));from build_independent_baseline import build_reference,apply_reference,VARS
checks=[]
def check(name,condition):
 checks.append({'check':name,'pass':bool(condition)})
 if not condition:raise AssertionError(name)
# Numerical provenance checks use executed evidence, not expected direction of ecological effects.
e=pd.read_csv(ROOT/'experiments/E00_original_reproduction/original_cell59_exact_comparison.csv');check('original AFB and95%CI exact to1e-12',e[['AFB_difference','ci_low_difference','ci_high_difference']].abs().to_numpy().max()<1e-12)
e=pd.read_csv(ROOT/'experiments/E00_original_reproduction/risk_exact_comparison.csv');check('all original context risks exact',e.risk_difference.abs().max()==0);check('all original risk sample counts exact',(e.recomputed_n==e.historical_n).all())
f=pd.read_csv(ROOT/'experiments/E00_original_reproduction/sample_flow.csv');check('main sample236995 and48843',int(f.iloc[-1].n_rows)==236995 and int(f.iloc[-1].n_grids)==48843)
e=pd.read_csv(ROOT/'experiments/E01_spatial_support/original_GEE_algorithm_exact_rerun.csv');check('original GEE distance algorithm rerun matches panel',e.difference_km.abs().max()<1e-6)
e=pd.read_csv(ROOT/'experiments/E01_spatial_support/independent_native_distance_check_original_threshold.csv');check('native geodesic contradicts original scale labels',(e.native_pixel_geodesic_min_km>e.original_table_km*2).all());check('metric projection matches independent geodesic within50m',abs(e.native_pixel_geodesic_min_km-e.UTM_fast_distance_km).max()<.05)
e=pd.read_csv(ROOT/'experiments/E02_confounding/harmonized/model_diagnostics.csv');check('harmonized models converge and full rank',e.converged.all() and (e['rank']==e.n_parameters).all())
e=pd.read_csv(ROOT/'experiments/E01_spatial_support/national_corrected_1km_contexts.csv');check('corrected 1km classification covers48843 unique grids',len(e)==48843 and e.grid_id.nunique()==48843)
check('corrected flags binary',all(set(e[c].unique())<={0,1} for c in ['near_perm_1km_corrected','near_seasonal_1km_corrected']))
e=pd.read_csv(ROOT/'experiments/E01_spatial_support/lightweight_vs_full_distance_QC.csv');check('lightweight and full distance classification agree17056 grids',(e.n_grids_compared==17056).all() and (e.mismatches_between64_and2048_pixel_search==0).all())
e=pd.read_csv(ROOT/'experiments/E01_spatial_support/corrected_1km/harmonized/model_diagnostics.csv');check('corrected models converge and full rank',e.converged.all() and (e['rank']==e.n_parameters).all())
check('heavy task pause flag preserved',(ROOT/'experiments/E01_spatial_support/heavy_pause.flag').exists())
for p in ROOT.rglob('*.py'):ast.parse(p.read_text(encoding='utf-8-sig'))
check('all delivered Python sources parse',True)
# Synthetic fixtures below test independence and missingness rules only; they are never scientific results.
h=pd.DataFrame([dict(grid_id='test',year=y,step16=1,**{v:float(y-2000) for v in VARS}) for y in range(2001,2018)]);h['EVI']=h.EVI/20;h['LAI']=h.LAI/2;ref=build_reference(h);check('17 independent historical values counted',int(ref.EVI_ref_count.iloc[0])==17)
cur=h.iloc[[0]].copy();cur['year']=2018;x=apply_reference(cur,ref);check('reference application produces finite z',np.isfinite(x.EVI_z_independent.iloc[0]))
rejected=False
try:build_reference(h,end=2018)
except ValueError:rejected=True
check('overlapping baseline prohibited',rejected)
rejected=False
try:build_reference(pd.concat([h,h.iloc[[0]]]))
except ValueError:rejected=True
check('duplicate historical grid-year-step prohibited',rejected)
check('figures have nonempty PDF outputs',all(p.stat().st_size>1000 and p.open('rb').read(4)==b'%PDF' for p in (ROOT/'figures').glob('*.pdf')))
(ROOT/'logs/validation.json').write_text(json.dumps({'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks},indent=2));print(json.dumps(checks,indent=2))
