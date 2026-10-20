from pathlib import Path
import sys,json
BASE=Path(__file__).resolve().parents[1];ROOT=BASE/'outputs/surface-water-vegetation-fragility';sys.path.insert(0,str(ROOT/'src'))
from audit_pipeline import *
source=(ROOT/'src/harmonized_confounding.py').read_text(encoding='utf8')
inject="""
_full_q=d.HDSI.quantile([.25,.75])
_ctx=pd.read_csv(ROOT/'experiments/E01_spatial_support/national_corrected_1km_contexts.csv')
d=d.merge(_ctx,on='grid_id',validate='many_to_one')
d['context']=np.where(d.near_perm_1km_corrected==1,'permanent',np.where(d.near_seasonal_1km_corrected==1,'seasonal','other'))
_before=dict(n_rows=len(d),n_grids=int(d.grid_id.nunique()))
d=d[~d.LC2018_grid_mode.isin([12,14])].copy()
(ROOT/'experiments/E02_confounding/corrected_noncrop').mkdir(exist_ok=True)
(ROOT/'experiments/E02_confounding/corrected_noncrop/exclusion_config.json').write_text(json.dumps({'excluded_landcover_classes':[12,14],'landcover_product':'2018 MODIS/061/MCD12Q1 grid spatial mode','before':_before,'after':{'n_rows':len(d),'n_grids':int(d.grid_id.nunique())},'HDSI_cutpoints_fixed_to_full_sample':_full_q.tolist(),'purpose':'agricultural-composition sensitivity, not identification of irrigation-free vegetation','remaining_confounding':['irrigation outside mapped classes','groundwater','soil','water management'],'GEE_jobs_submitted':0},indent=2))
"""
source=source.replace("out='experiments/E02_confounding/harmonized/'","out='experiments/E02_confounding/corrected_noncrop/'")
source=source.replace("fields=['elevation'",inject+"\nfields=['elevation'")
source=source.replace("q=d.HDSI.quantile([.25,.75]);print", "q=_full_q;print")
source=source.replace('r,a=group_effects(s);','r,a=group_effects(s,q_fixed=_full_q);')
source=source.replace('original manuscript panel; corrected-distance sensitivity separate','corrected 30m UTM1km contexts; excluded2018 modal classes12and14')
exec(compile(source,'<corrected_noncrop_models>','exec'),{'__name__':'__main__','__file__':str(ROOT/'src/harmonized_confounding.py')})
print('CORRECTED NONCROP SENSITIVITY COMPLETE',flush=True)
