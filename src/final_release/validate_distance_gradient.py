from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.sandwich_covariance import cov_cluster
from scipy.special import expit
from patsy import dmatrix,build_design_matrices
sys.path.insert(0,str(Path('src').resolve()))
from audit_pipeline import load_panel,group_effects,ROOT
out=Path('outputs/distance_validation');out.mkdir(parents=True,exist_ok=True)
_,d=load_panel();q=d.HDSI.quantile([.25,.75]);cpath=ROOT/'experiments/E01_spatial_support/centroid_export_received.csv';c=pd.read_csv(cpath)
assert not c.grid_id.duplicated().any()
cov=pd.read_csv(ROOT/'experiments/E02_confounding/harmonized_landcover_climate_covariates.csv')
d=d.merge(cov,on='grid_id',validate='many_to_one');g=d.drop_duplicates('grid_id').copy();g['available']=g.grid_id.isin(c.grid_id)
g['zone']=np.floor((g.lon+180)/6).astype(int)+1
s=d.merge(c.drop(columns=['system:index','.geo'],errors='ignore'),on='grid_id',validate='many_to_one')
n_before=len(s);n_grid_before=s.grid_id.nunique();s=s[s.dist_perm_corrected_km_in_range.eq(1)&s.dist_perm_corrected_km.notna()].copy()
meta={'national_rows':len(d),'national_grids':d.grid_id.nunique(),'received_grids':len(c),'received_zones':sorted(c.zone.unique().tolist()),'rows_before_range_filter':n_before,'grids_before_range_filter':n_grid_before,'rows_in_range':len(s),'grids_in_range':s.grid_id.nunique(),'q25':q.iloc[0],'q75':q.iloc[1],'search_limit_km':61.44,'bootstrap':300,'seed':20261020,'heavy_tasks':'remain paused; no remote computation','partial_sample':'completed UTM exports, not random national sample','input_sha256':hashlib.sha256(cpath.read_bytes()).hexdigest()}
print(json.dumps(meta),flush=True)
g.groupby(['zone','available']).size().unstack(fill_value=0).to_csv(out/'coverage_by_utm_zone.csv')
fields=['lon','lat','elevation','slope','climate_annual_temp_C','climate_annual_precip_mm']
rep=[]
for col in fields:
 a=g.loc[g.available,col].dropna();b=g.loc[~g.available,col].dropna();rep.append(dict(variable=col,available_mean=a.mean(),unavailable_mean=b.mean(),national_mean=g[col].mean(),available_sd=a.std(),unavailable_sd=b.std(),SMD=(a.mean()-b.mean())/np.sqrt((a.var()+b.var())/2)))
pd.DataFrame(rep).to_csv(out/'sample_representativeness.csv',index=False)
pd.crosstab(g.LC2018_grid_mode,g.available,normalize='columns').to_csv(out/'landcover_composition.csv')
g[['grid_id','lon','lat','available','zone']].to_csv(out/'sample_coverage.csv',index=False)
s.drop_duplicates('grid_id')[['grid_id','lon','lat','zone','dist_perm_water_km','dist_perm_corrected_km']].to_csv(out/'paired_distances.csv',index=False)
# Same cells/years, fixed national thresholds, identical resampling seed across distance definitions.
schemes={'fine':[0,.25,.5,1,2,3,5,7.5,10,15,30,61.44],'coarse':[0,1,5,10,20,40,61.44]}
quant=np.quantile(s.drop_duplicates('grid_id').dist_perm_corrected_km,[0,.2,.4,.6,.8,1]);quant[0]=0;quant[-1]=61.44;schemes['grid_quintiles']=np.unique(quant).tolist()
allbins=[]
for metric in ['dist_perm_water_km','dist_perm_corrected_km']:
 for name,edges in schemes.items():
  if name=='grid_quintiles' and metric=='dist_perm_water_km':continue
  z=s.copy();labels=[f'{a:g}-{b:g}' for a,b in zip(edges[:-1],edges[1:])];z['band']=pd.cut(z[metric],[-1e-9]+edges[1:],labels=labels).astype('string');r,_=group_effects(z,group='band',q_fixed=q)
  r['metric']=metric;r['scheme']=name;r['bin_index']=r.group.map(dict(zip(labels,range(len(labels)))));r['supported']=(r.n_low>=30)&(r.n_high>=30)&(r.n_grid>=30);allbins.append(r)
  print('bins',metric,name,flush=True)
bins=pd.concat(allbins);bins.to_csv(out/'distance_band_fragility.csv',index=False)
# Near/far comparisons use paired resampling; retain middle distances in the eligible sample but not contrast.
cr=[]
for near,far in [(1,5),(1,10),(5,20)]:
 z=s.copy();z['context']=np.where(z.dist_perm_corrected_km<=near,'permanent',np.where(z.dist_perm_corrected_km>=far,'far','middle'))
 for cluster,deg in [('grid_id',0),('block',1),('block',2),('block',5)]:
  if deg:z['block']=(np.floor(z.lon/deg).astype(int).astype(str)+'_'+np.floor(z.lat/deg).astype(int).astype(str))
  r,a=group_effects(z,cluster=cluster,q_fixed=q);a=a[a.comparison.eq('far')].copy();a['near_max_km']=near;a['far_min_km']=far;a['block_degrees']=deg;cr.append(a)
pd.concat(cr).to_csv(out/'unadjusted_far_minus_near.csv',index=False)
print('contrasts complete',flush=True)
# Adjusted nonlinear distance association. Fixed national HDSI cutpoints; common empirical standardization.
cc=fields+['year','LC2018_grid_mode'];m=s[s[cc].notna().all(axis=1)&s.LC2018_grid_mode.between(1,17)&((s.HDSI<q.iloc[0])|(s.HDSI>q.iloc[1]))].copy()
h=(m.HDSI>q.iloc[1]).to_numpy(float);dist=m.dist_perm_corrected_km.to_numpy();zz=np.log1p(dist);upper=np.log1p(61.44)
B=dmatrix('0 + bs(z, df=5, degree=3, include_intercept=False, lower_bound=0, upper_bound=upper)',{'z':zz,'upper':upper});info=B.design_info;basis=np.asarray(B)
core=np.column_stack([np.ones(len(m)),h,basis,basis*h[:,None]]);controls=[];names=[]
for col in fields:
 v=m[col].to_numpy();v=np.log1p(v) if col=='climate_annual_precip_mm' else v;controls.append((v-v.mean())/v.std(ddof=1));names.append(col)
controls.extend([controls[0]**2,controls[1]**2,controls[0]*controls[1]])
for col in ['year','LC2018_grid_mode']:
 for val in sorted(m[col].unique())[1:]:controls.append(m[col].eq(val).to_numpy(float))
controls=np.column_stack(controls);x=np.column_stack([core,controls]);rank=np.linalg.matrix_rank(x)
assert rank==x.shape[1],(rank,x.shape)
clusters={'grid':m.grid_id.to_numpy()}
for deg in [1,2,5]:clusters[f'block_{deg}deg']=(np.floor(m.lon/deg).astype(int).astype(str)+'_'+np.floor(m.lat/deg).astype(int).astype(str)).to_numpy()
sg=m.drop_duplicates('grid_id');lo,hi=np.quantile(sg.dist_perm_corrected_km,[.05,.95]);points=np.unique(np.r_[np.linspace(lo,hi,35),1,5,10,20]);points=points[(points>=dist.min())&(points<=dist.max())]
pred=[];contrast=[];diag=[];tests=[]
for event in ['EVI','LAI','combined']:
 fit=sm.GLM(m[event].to_numpy(),x,family=sm.families.Binomial()).fit(maxiter=100)
 Vs={k:cov_cluster(fit,v) for k,v in clusters.items()};b=fit.params;cache={}
 for dk in points:
  bb=np.asarray(build_design_matrices([info],{'z':np.array([np.log1p(dk)]),'upper':upper})[0])[0];mu=[];grad=[]
  for hh in [0,1]:
   vcore=np.r_[1,hh,bb,bb*hh];p=expit(vcore@b[:12]+controls@b[12:]);w=p*(1-p);mu.append(p.mean());grad.append(np.r_[vcore*w.mean(),(controls*w[:,None]).mean(axis=0)])
  est=mu[1]-mu[0];gr=grad[1]-grad[0];cache[dk]=(est,gr)
  for kind,V in Vs.items():
   se=np.sqrt(max(0,gr@V@gr));pred.append(dict(response=event,distance_km=dk,fragility=est,risk_low=mu[0],risk_high=mu[1],ci_low=est-1.96*se,ci_high=est+1.96*se,uncertainty=kind,n_local_grids=int((abs(sg.dist_perm_corrected_km-dk)<=1).sum())))
 for a,bk in [(1,5),(1,10),(5,20)]:
  if a not in cache or bk not in cache:continue
  delta=cache[bk][0]-cache[a][0];gr=cache[bk][1]-cache[a][1]
  for kind,V in Vs.items():
   se=np.sqrt(max(0,gr@V@gr));contrast.append(dict(response=event,near_km=a,far_km=bk,fragility_difference=delta,ci_low=delta-1.96*se,ci_high=delta+1.96*se,uncertainty=kind))
 for kind,V in Vs.items():
  R=np.zeros((5,x.shape[1]));R[:,7:12]=np.eye(5);wt=fit.wald_test(R,cov_p=V,scalar=True);tests.append(dict(response=event,uncertainty=kind,interaction_wald=float(wt.statistic),p_value=float(wt.pvalue),df=5))
 diag.append(dict(response=event,converged=bool(fit.converged),n_rows=len(m),n_grids=m.grid_id.nunique(),rank=int(rank),parameters=x.shape[1],condition=float(np.linalg.cond(x))))
 print('model',event,diag[-1],flush=True)
pred=pd.DataFrame(pred);contrast=pd.DataFrame(contrast);pred.to_csv(out/'adjusted_spline_predictions.csv',index=False);contrast.to_csv(out/'adjusted_distance_contrasts.csv',index=False);pd.DataFrame(diag).to_csv(out/'model_diagnostics.csv',index=False);pd.DataFrame(tests).to_csv(out/'interaction_tests.csv',index=False)
meta.update(schemes=schemes,model='high * bs(log1p(corrected distance), df=5) + scaled lon + lat + elevation + slope + MAT + log1p(MAP) + lon^2 + lat^2 + lon*lat + C(year) + C(2018 modal landcover)',standardization='common empirical covariate distribution in complete-case low/high subset',curve_range_5_95=[lo,hi],model_block_counts={k:len(np.unique(v)) for k,v in clusters.items()},CI='pointwise normal delta method, clustered sandwich; conditional on empirical standardization distribution',limitations=['partial nonrandom UTM coverage','centroid distances, not within-cell vegetation distances','no grid realignment or resolution reaggregation','spatial block sizes are sensitivity choices, not estimated correlation lengths','multiple exploratory contrasts, no multiplicity correction'])
(out/'analysis_config.json').write_text(json.dumps(meta,indent=2))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
f,axs=plt.subplots(1,3,figsize=(14,4.6),layout='constrained')
for ax,e in zip(axs,['EVI','LAI','combined']):
 for metric,color,label in [('dist_perm_water_km','#d48b30','Original distance'),('dist_perm_corrected_km','#217b9c','Corrected distance')]:
  z=bins[(bins.response==e)&(bins.metric==metric)&(bins.scheme=='coarse')].sort_values('bin_index'); xx=z.bin_index.to_numpy()+(-.08 if metric=='dist_perm_water_km' else .08);ax.errorbar(xx,z.fragility,yerr=[z.fragility-z.ci_low,z.ci_high-z.fragility],fmt='o-',color=color,label=label,capsize=3)
 ax.set(xticks=range(6),xticklabels=['0–1','1–5','5–10','10–20','20–40','40–61.44'],xlabel='Centroid distance class (km)',ylabel='Fragility increase',title=e);ax.tick_params(axis='x',rotation=35);ax.axhline(0,color='grey',lw=.6)
axs[0].legend(fontsize=8);f.suptitle('Same partial sample; fixed national HDSI thresholds; grid-cluster 95% intervals')
for ext in ['png','pdf']:f.savefig(out/f'Fig_C2_paired_distance_bands.{ext}',dpi=200)
plt.close(f)
f,axs=plt.subplots(1,3,figsize=(14,4.6),layout='constrained')
for ax,e in zip(axs,['EVI','LAI','combined']):
 z=pred[(pred.response==e)&pred.uncertainty.eq('grid')].sort_values('distance_km');ax.plot(z.distance_km,z.fragility,color='#217b9c');ax.fill_between(z.distance_km,z.ci_low,z.ci_high,color='#217b9c',alpha=.18,label='Grid-cluster 95% CI');v=pred[(pred.response==e)&pred.uncertainty.eq('block_2deg')].sort_values('distance_km');ax.plot(v.distance_km,v.ci_low,'--',color='#b47836',lw=1);ax.plot(v.distance_km,v.ci_high,'--',color='#b47836',lw=1,label='2° block 95% CI');ax.set(title=e,xlabel='Corrected centroid distance (km)',ylabel='Adjusted fragility');ax.axhline(0,color='grey',lw=.6)
axs[0].legend(fontsize=8);f.suptitle('Adjusted distance association in completed UTM zones; pointwise intervals')
for ext in ['png','pdf']:f.savefig(out/f'Fig_C2_adjusted_distance.{ext}',dpi=200)
plt.close(f)
f,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained');axs[0].scatter(g.lon,g.lat,s=2,c=np.where(g.available,'#217b9c','#d8d8d8'),rasterized=True);axs[0].set(xlabel='Longitude',ylabel='Latitude',title='Completed distance coverage (blue)');axs[0].set_aspect('equal');rep=pd.DataFrame(rep);axs[1].barh(rep.variable,rep.SMD,color='#217b9c');axs[1].axvline(0,color='grey');axs[1].set(xlabel='Standardized mean difference\nAvailable minus unavailable grids',title='Sample composition')
for ext in ['png','pdf']:f.savefig(out/f'Fig_C2_sample_coverage.{ext}',dpi=200)
plt.close(f)
print(contrast.to_string(index=False),flush=True)
print('DONE',flush=True)
