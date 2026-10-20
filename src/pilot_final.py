from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
out='experiments/E01_spatial_support/'
m=pd.read_csv(ROOT/out/'pilot_within_cell_metrics.csv');g=pd.read_csv(ROOT/out/'pilot_selected_grids.csv');cent=pd.read_csv(ROOT/out/'pilot_same_raster_centroid_distances.csv');a=m[m.alignment=='original'].merge(g[['grid_id','dist_perm_water_km']],on='grid_id').merge(cent,on=['grid_id','scenario']);save(a,out+'pilot_centroid_comparison.csv')
summary=[]
for scenario,t in a.groupby('scenario'):
 summary.append(dict(scenario=scenario,n_grid=len(t),median_water_fraction=t.water_fraction.median(),median_original_table_centroid_km=t.dist_perm_water_km.median(),median_same_raster_centroid_km=t.recomputed_centroid_km.median(),median_mean_vegetation_distance_km=t.mean_vegetation_distance_km.median(),median_within_minus_same_raster_centroid_km=(t.mean_vegetation_distance_km-t.recomputed_centroid_km).median(),spearman_same_raster_centroid_vs_within=t[['recomputed_centroid_km','mean_vegetation_distance_km']].corr(method='spearman').iloc[0,1],min_distance_coverage=t.vegetation_distance_coverage.min()))
save(pd.DataFrame(summary),out+'pilot_summary.csv');_,d=load_panel();s=d[d.grid_id.isin(g.grid_id)];results=[];frags=[]
for scenario,t in a.groupby('scenario'):
 z=s.merge(t[['grid_id','recomputed_centroid_km','mean_vegetation_distance_km','water_fraction']],on='grid_id');metrics=['recomputed_centroid_km','mean_vegetation_distance_km']+(['dist_perm_water_km'] if scenario=='main' else [])
 for metric in metrics:
  z['pilot_context']=np.where(z[metric]<=1,'permanent','comparison_gt1km');r,af=group_effects(z,group='pilot_context');r['distance_metric']=metric;r['scenario']=scenario;frags.append(r)
  if af.empty:
   af=pd.DataFrame([dict(comparison='comparison_gt1km',response=e,AFB=np.nan,ci_low=np.nan,ci_high=np.nan,n_total=len(z),n_grid=z.grid_id.nunique(),valid_draws=0,status='not_estimable_no_near_group') for e in ['EVI','LAI','combined']])
  else:af['status']='pilot_estimate'
  af['distance_metric']=metric;af['scenario']=scenario;af['warning']='24-grid regional pilot; within-pilot HDSI quartiles; does not validate national distance conclusions';results.append(af)
save(pd.concat(results),out+'pilot_distance_AFB.csv');save(pd.concat(frags),out+'pilot_distance_fragility.csv');print(pd.DataFrame(summary).to_string(index=False))
