from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
m=pd.read_csv(ROOT/'experiments/E01_spatial_support/pilot_within_cell_metrics.csv');g=pd.read_csv(ROOT/'experiments/E01_spatial_support/pilot_selected_grids.csv');a=m[m.alignment=='original'].merge(g[['grid_id','dist_perm_water_km']],on='grid_id');save(a,'experiments/E01_spatial_support/pilot_centroid_comparison.csv');summary=[]
for scenario,t in a.groupby('scenario'):
 summary.append(dict(scenario=scenario,n_grid=len(t),median_water_fraction=t.water_fraction.median(),median_centroid_km=t.dist_perm_water_km.median(),median_mean_vegetation_distance_km=t.mean_vegetation_distance_km.median(),median_distance_difference_km=(t.mean_vegetation_distance_km-t.dist_perm_water_km).median(),spearman_centroid_vs_within=t[['dist_perm_water_km','mean_vegetation_distance_km']].corr(method='spearman').iloc[0,1],min_distance_coverage=t.vegetation_distance_coverage.min()))
save(pd.DataFrame(summary),'experiments/E01_spatial_support/pilot_summary.csv');orig=m[m.alignment=='original'];shift=m[m.alignment!='original'];align=orig.merge(shift,on=['grid_id','scenario'],suffixes=('_original','_shift'));align['mean_distance_change_km']=align.mean_vegetation_distance_km_shift-align.mean_vegetation_distance_km_original;align['water_fraction_change']=align.water_fraction_shift-align.water_fraction_original;save(align,'experiments/E01_spatial_support/pilot_alignment_metric_comparison.csv')
_,d=load_panel();s=d[d.grid_id.isin(g.grid_id)].copy();q=d.HDSI.quantile([.25,.75]);results=[];frags=[]
for scenario,t in a.groupby('scenario'):
 z=s.merge(t[['grid_id','mean_vegetation_distance_km','water_fraction']],on='grid_id',how='inner');
 for metric in ['dist_perm_water_km','mean_vegetation_distance_km']:
  # Binary pilot distance contrast only: comparator >1 km versus near <=1 km. Seasonal water not redefined here.
  z['pilot_context']=np.where(z[metric]<=1,'permanent','comparison_gt1km');r,af=group_effects(z,group='pilot_context');r['distance_metric']=metric;r['scenario']=scenario;af['distance_metric']=metric;af['scenario']=scenario;af['warning']='24-grid regional pilot; local HDSI quartiles, broad finite-sample uncertainty; not national inference';results.append(af);frags.append(r)
save(pd.concat(results),'experiments/E01_spatial_support/pilot_distance_AFB.csv');save(pd.concat(frags),'experiments/E01_spatial_support/pilot_distance_fragility.csv')
print(pd.DataFrame(summary).to_string(index=False));print(m.vegetation_distance_coverage.describe().to_string());print('pilot n rows',len(s),'grids',s.grid_id.nunique())
