from pathlib import Path
import sys,json,datetime
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
cloud=Path(__import__('os').environ.get('AQUATIC_EXPORT_DIR','data/exports/aquatic_spatial_metrics'));out=ROOT/'experiments/E01_spatial_support';_,d=load_panel();need=set(d.grid_id.unique());summary=[]
for prefix,kind in [('E01_corrected_centroids_20261020_UTM','centroid'),('E01_within_cell_20261020_UTM','within_cell')]:
 files=sorted(cloud.glob(prefix+'*.csv'));a=pd.concat([pd.read_csv(p) for p in files],ignore_index=True) if files else pd.DataFrame();n=len(a);rec={'type':kind,'n_files':len(files),'n_rows':n,'expected_grids':len(need),'complete':False}
 if len(a):
  assert not a.grid_id.duplicated().any();assert set(a.grid_id).issubset(need);a.to_csv(out/(kind+'_export_received.csv'),index=False);rec['complete']=set(a.grid_id)==need
  if rec['complete']:
   if kind=='centroid':a.to_csv(out/'national_corrected_centroids.csv',index=False)
   else:
    a['permanent_water_fraction']=a.permanent_area/a.cell_area;a['vegetation_distance_coverage']=a.distance_supported_area/a.vegetation_area;a['vegetation_mean_perm_distance_km']=a.distance_area_sum/a.distance_supported_area/1000;a['valid_support']=(a.vegetation_area>0)&(a.vegetation_distance_coverage>=.999);a.to_csv(out/'national_within_cell_metrics.csv',index=False)
 summary.append(rec)
(out/'export_collection_status.json').write_text(json.dumps({'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'exports':summary},indent=2));print(json.dumps(summary,indent=2))
