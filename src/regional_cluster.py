from pathlib import Path
import sys,json,ast
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
_,d=load_panel();mapping=pd.read_csv(ROOT/'experiments/E00_original_reproduction/sample_region_map.csv');d=d.merge(mapping[['grid_id','year','macro_region']],on=['grid_id','year'],validate='one_to_one');q=d.HDSI.quantile([.25,.75]);rows=[];frag=[]
for region,s in d.groupby('macro_region'):
 r,a=group_effects(s,q_fixed=q);a['region']=region;r['region']=region;counts=[]
 for ctx in ['permanent','other']:
  sub=s[s.context==ctx];counts += [sum(sub.HDSI<q.iloc[0]),sum(sub.HDSI>q.iloc[1])]
 a['insufficient_other_vs_permanent']=min(counts)<30;a['minimum_context_exposure_n']=min(counts);rows.append(a);frag.append(r)
save(pd.concat(rows),'experiments/E05_figure_audit/regional_AFB_grid_cluster.csv');save(pd.concat(frag),'experiments/E05_figure_audit/regional_fragility_grid_cluster.csv');print('regional bootstrap complete',flush=True)
