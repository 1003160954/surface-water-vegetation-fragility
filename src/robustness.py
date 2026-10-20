from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
raw,d=load_panel();out='experiments/E04_robustness/';results=[];risks=[]
components={'VPD':d.summer_VPD_z,'LST':d.summer_LST_Day_z,'SM_deficit':-d.summer_SM_rootproxy_z,'Precip_deficit':-d.summer_Precip_z};C=pd.DataFrame(components,index=d.index)
variants={'equal':C.mean(axis=1)}
Z=(C-C.mean())/C.std(ddof=1);u,sv,vt=np.linalg.svd(Z.to_numpy(),full_matrices=False);loadings=vt[0];pc=Z.to_numpy()@loadings
if np.corrcoef(pc,d.HDSI)[0,1]<0:pc=-pc;loadings=-loadings
variants['PCA']=pc;save(pd.DataFrame({'component':list(C),'PC1_loading':loadings,'input_mean':C.mean().values,'input_sd':C.std(ddof=1).values}),out+'PCA_loadings.csv');save(C.corr().rename_axis('component').reset_index(),out+'component_correlations.csv')
for col in C:variants['omit_'+col]=C.drop(columns=col).mean(axis=1);variants['only_'+col]=C[col]
for name,v in variants.items():
 s=d.copy();s['variant']=v;r,a=group_effects(s,stress='variant');r['scenario']=name;a['scenario']=name;results.append(a);risks.append(r)
save(pd.concat(results),out+'HDSI_variant_AFB.csv');save(pd.concat(risks),out+'HDSI_variant_fragility.csv');print('variants complete',flush=True)
results=[];risks=[]
for t in [-.25,-.5,-.75,-1.]:
 s=d.copy();s['EVI']=(s.late_summer_EVI_z<t).astype(int);s['LAI']=(s.late_summer_LAI_z<t).astype(int);s['combined']=s[['EVI','LAI']].max(axis=1);r,a=group_effects(s);r['threshold']=t;a['threshold']=t;results.append(a);risks.append(r)
save(pd.concat(results),out+'event_threshold_AFB.csv');save(pd.concat(risks),out+'event_threshold_fragility.csv');print('thresholds complete',flush=True)
results=[];risks=[]
for y in sorted(d.year.unique()):
 # LOO here excludes an analysis year; it does not recalculate step-wise anomaly baselines.
 s=d[d.year!=y].copy();r,a=group_effects(s);r['excluded_year']=y;a['excluded_year']=y;results.append(a);risks.append(r)
save(pd.concat(results),out+'leave_analysis_year_out_AFB.csv');save(pd.concat(risks),out+'leave_analysis_year_out_fragility.csv');print('LOO complete',flush=True)
results=[];risks=[]
for deg in [1,2,5]:
 s=d.copy();s['spatial_block']=(np.floor(s.lon/deg).astype(int).astype(str)+'_'+np.floor(s.lat/deg).astype(int).astype(str));r,a=group_effects(s,cluster='spatial_block');r['block_degrees']=deg;a['block_degrees']=deg;results.append(a);risks.append(r)
save(pd.concat(results),out+'spatial_block_AFB.csv');save(pd.concat(risks),out+'spatial_block_fragility.csv');print('blocks complete',flush=True)
# Symmetric binary 8-nearest-neighbour graph; coordinates approximated in km at study mean latitude.
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
g=d.groupby('grid_id').agg(lon=('lon','first'),lat=('lat','first'),HDSI=('HDSI','mean'),EVI=('EVI','mean'),LAI=('LAI','mean'),combined=('combined','mean')).reset_index();xy=np.column_stack([g.lon*111.32*np.cos(np.deg2rad(g.lat.mean())),g.lat*111.32]);idx=cKDTree(xy).query(xy,k=9)[1][:,1:];N=len(g);W=csr_matrix((np.ones(N*8),(np.repeat(np.arange(N),8),idx.ravel())),shape=(N,N));W=W.maximum(W.T);S0=W.sum();rng=np.random.default_rng(SEED);rows=[]
for col in ['HDSI','EVI','LAI','combined']:
 z=g[col].to_numpy()-g[col].mean();den=z@z;I=float(N/S0*(z@(W@z))/den);perms=[]
 for _ in range(199):
  zz=rng.permutation(z);perms.append(float(N/S0*(zz@(W@zz))/den))
 p=(1+np.sum(np.abs(np.asarray(perms)+1/(N-1))>=abs(I+1/(N-1))))/200
 rows.append(dict(variable=col,Morans_I=I,p_permutation_two_sided=p,n_grid=N,k=8,graph='symmetric_binary_knn',permutations=199,seed=SEED,interpretation='grid means; HDSI mean near zero by construction'))
save(pd.DataFrame(rows),out+'Morans_I.csv')
(ROOT/out/'run_config.json').write_text(json.dumps({'seed':SEED,'bootstrap':300,'block_sizes_degrees':[1,2,5],'block_origin_lon_lat':[0,0],'block_grid':'floor(lon/size), floor(lat/size)','block_size_warning':'degrees imply latitude-varying physical size; sensitivity, not estimated decorrelation range','quantiles':'recomputed per variant and per leave-year-out sample; held fixed across bootstrap draws','PCA':'PC1 of sample-standardized seasonal stress components, sign aligned with equal index','LOO':'analysis-year exclusion only, original anomaly baseline retained','permanence_threshold_AFB':'not run: alternative raster-derived distance surfaces unavailable'},indent=2))
print('E04 FINISHED',flush=True)
