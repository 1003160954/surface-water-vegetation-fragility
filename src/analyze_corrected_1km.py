from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
_,d=load_panel();c=pd.read_csv(ROOT/'experiments/E01_spatial_support/national_corrected_1km_contexts.csv');assert c.grid_id.nunique()==d.grid_id.nunique() and not c.grid_id.duplicated().any();assert set(c.near_perm_1km_corrected.unique())<={0,1};assert set(c.near_seasonal_1km_corrected.unique())<={0,1};z=d.merge(c,on='grid_id',validate='many_to_one');z['original_context']=z.context;z['context']=np.where(z.near_perm_1km_corrected==1,'permanent',np.where(z.near_seasonal_1km_corrected==1,'seasonal','other'));out='experiments/E01_spatial_support/corrected_1km/';r,a=group_effects(z);save(r,out+'context_fragility.csv');save(a,out+'AFB.csv');g=z.drop_duplicates('grid_id');save(pd.crosstab(g.original_context,g.context).rename_axis('original_context').reset_index(),out+'context_reclassification_grid.csv');save(pd.crosstab(z.original_context,z.context).rename_axis('original_context').reset_index(),out+'context_reclassification_grid_year.csv');print('CORRECTED1KM AFB',a.to_string(index=False),flush=True)
inject="\n_ctx=pd.read_csv(ROOT/'experiments/E01_spatial_support/national_corrected_1km_contexts.csv');d=d.merge(_ctx,on='grid_id',validate='many_to_one');d['context']=np.where(d.near_perm_1km_corrected==1,'permanent',np.where(d.near_seasonal_1km_corrected==1,'seasonal','other'))\n"
source=(Path(__file__).resolve().parent/'harmonized_confounding.py').read_text(encoding='utf-8-sig').replace("out='experiments/E02_confounding/harmonized/'","out='experiments/E01_spatial_support/corrected_1km/harmonized/'").replace("fields=['elevation'",inject+"fields=['elevation'").replace('original manuscript panel; corrected-distance sensitivity separate','corrected30m UTM 1km centroid contexts, actual original threshold80/10')
exec(compile(source,'<harmonized_corrected1km>','exec'),{'__name__':'__main__','__file__':__file__})
source=(Path(__file__).resolve().parent/'common_landscape.py').read_text(encoding='utf-8-sig').replace("_,d=load_panel();cv=","_,d=load_panel();"+inject+"cv=").replace('experiments/E02_confounding/harmonized/','experiments/E01_spatial_support/corrected_1km/harmonized/').replace('original distance grouping','corrected30m UTM1km contexts')
try:
 exec(compile(source,'<landscape_corrected1km>','exec'),{'__name__':'__main__','__file__':__file__})
except Exception as _error:
 (ROOT/'experiments/E01_spatial_support/corrected_1km/landscape_error.json').write_text(json.dumps({'status':'not_estimable_or_failed','reason':str(_error)},indent=2))
 print('Landscape step:',str(_error),flush=True)
source=(Path(__file__).resolve().parent/'robustness.py').read_text(encoding='utf-8-sig');source=source.split('# Symmetric binary')[0];source=source.replace("raw,d=load_panel();out=", "raw,d=load_panel();"+inject+"out=").replace('experiments/E04_robustness/','experiments/E01_spatial_support/corrected_1km/robustness/')
exec(compile(source,'<robustness_corrected1km>','exec'),{'__name__':'__main__','__file__':__file__})
print('CORRECTED1KM ANALYSES COMPLETE',flush=True)
