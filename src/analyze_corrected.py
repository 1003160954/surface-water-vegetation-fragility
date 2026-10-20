from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
_,d=load_panel();base=ROOT/'experiments/E01_spatial_support';cp=base/'national_corrected_centroids.csv'
if cp.exists():
 c=pd.read_csv(cp);assert c.grid_id.nunique()==d.grid_id.nunique() and not c.grid_id.duplicated().any();z=d.merge(c,on='grid_id',validate='many_to_one');z['original_context']=z.context;z['context']=np.where(z.dist_perm_corrected_km<=1,'permanent',np.where(z.dist_seasonal_corrected_km<=1,'seasonal','other'));out='experiments/E01_spatial_support/corrected_centroid/';r,a=group_effects(z);save(r,out+'context_fragility.csv');save(a,out+'AFB.csv');save(pd.crosstab(z.original_context,z.context).rename_axis('original_context').reset_index(),out+'context_reclassification_grid_year.csv');g=z.drop_duplicates('grid_id');save(pd.crosstab(g.original_context,g.context).rename_axis('original_context').reset_index(),out+'context_reclassification_grid.csv');save(g[['grid_id','lon','lat','dist_perm_water_km','dist_perm_corrected_km','dist_seasonal_water_km','dist_seasonal_corrected_km','original_context','context']],out+'distance_comparison.csv')
 for deg in [1,2,5]:
  zz=z.copy();zz['block']=(np.floor(zz.lon/deg).astype(int).astype(str)+'_'+np.floor(zz.lat/deg).astype(int).astype(str));r,a=group_effects(zz,cluster='block');r['block_degrees']=deg;a['block_degrees']=deg;save(a,out+'AFB_spatial_'+str(deg)+'deg.csv')
 # Covariate model sensitivity is fitted by the same declared harmonized model, only distance grouping changes.
 source=(Path(__file__).resolve().parent/'harmonized_confounding.py').read_text(encoding='utf-8-sig');source=source.replace("out='experiments/E02_confounding/harmonized/'","out='experiments/E01_spatial_support/corrected_centroid/harmonized/'")
 inject="\n_c=pd.read_csv(ROOT/'experiments/E01_spatial_support/national_corrected_centroids.csv');d=d.merge(_c,on='grid_id',validate='many_to_one');d['context']=np.where(d.dist_perm_corrected_km<=1,'permanent',np.where(d.dist_seasonal_corrected_km<=1,'seasonal','other'))\n"
 source=source.replace("fields=['elevation'",inject+"fields=['elevation'")
 source=source.replace("original manuscript panel; corrected-distance sensitivity separate","corrected 30m UTM centroid distances, actual original water threshold80/10")
 exec(compile(source,'<harmonized_corrected_distance>','exec'),{'__name__':'__main__'})
 print('Corrected centroid results analyzed',flush=True)
wp=base/'national_within_cell_metrics.csv'
if wp.exists():
 w=pd.read_csv(wp);s=d.merge(w,on='grid_id',validate='many_to_one');s=s[s.valid_support].copy();out='experiments/E01_spatial_support/within_cell/';sample=[]
 for metric in ['dist_perm_water_km','vegetation_mean_perm_distance_km']:
  z=s.copy();z['context']=np.where(z[metric]<=1,'permanent','far_mean_distance');r,a=group_effects(z);r['distance_metric']=metric;a['distance_metric']=metric;save(r,out+metric+'_fragility.csv');save(a,out+metric+'_AFB.csv');sample.append(dict(metric=metric,n_rows=len(z),n_grid=z.grid_id.nunique(),n_near_rows=int((z.context=='permanent').sum()),n_near_grid=z[z.context=='permanent'].grid_id.nunique()))
 z=s.copy();z['distance_band']=pd.cut(z.vegetation_mean_perm_distance_km,[-1e-9,.25,.5,1,2,5,10,20,40,float('inf')]).astype(str);r,_=group_effects(z,group='distance_band');save(r,out+'vegetation_distance_band_fragility.csv');z['water_fraction_band']=pd.cut(z.permanent_water_fraction,[-1e-9,0,.01,.05,1],labels=['0','0-1%','1-5%','>5%']).astype(str);r,_=group_effects(z,group='water_fraction_band');save(r,out+'water_fraction_fragility.csv');save(pd.DataFrame(sample),out+'same_support_sample_counts.csv');save(pd.DataFrame([dict(raw_grids=len(w),valid_support_grids=int(w.valid_support.sum()),excluded_grids=int((~w.valid_support).sum()),coverage_min=float(w.vegetation_distance_coverage.min()))]),out+'support_qc.csv');print('Within-cell results analyzed',flush=True)

# Refresh user-facing evidence only after a successful executed analysis.
import subprocess as _subprocess
_subprocess.run([sys.executable,"-X","utf8",str(ROOT/"src/build_report.py")],check=True)
