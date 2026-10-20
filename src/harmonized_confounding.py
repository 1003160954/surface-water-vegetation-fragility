from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
import statsmodels.api as sm
from scipy.special import expit
_,d=load_panel();out='experiments/E02_confounding/harmonized/';new=pd.read_csv(ROOT/'experiments/E02_confounding/harmonized_landcover_climate_covariates.csv');d=d.merge(new,on='grid_id',validate='many_to_one')
fields=['elevation','slope','lon','lat','year','LC2018_grid_mode','climate_annual_temp_C','climate_annual_precip_mm'];cc=d[fields].notna().all(axis=1)&d.LC2018_grid_mode.between(1,17);c=d[cc].copy();q=d.HDSI.quantile([.25,.75]);print('eligible',len(d),'complete_case',len(c),flush=True)
for name,s in [('M0',d),('M1',c)]:
 r,a=group_effects(s);r['model']=name;a['model']=name;save(r,out+name+'_fragility.csv');save(a,out+name+'_AFB.csv')
s=c[(c.HDSI<q.iloc[0])|(c.HDSI>q.iloc[1])].copy();h=(s.HDSI>q.iloc[1]).to_numpy(float);seas=s.context.eq('seasonal').to_numpy(float);other=s.context.eq('other').to_numpy(float);x=np.column_stack([np.ones(len(s)),h,seas,other,h*seas,h*other]);names=['Intercept','high','seasonal','other','high_seasonal','high_other']
for col in ['elevation','slope','lon','lat','climate_annual_temp_C','climate_annual_precip_mm']:
 v=s[col].to_numpy();v=np.log1p(v) if col=='climate_annual_precip_mm' else v;x=np.column_stack([x,(v-v.mean())/v.std(ddof=1)]);names.append(col+'_scaled')
x=np.column_stack([x,x[:,8]**2,x[:,9]**2,x[:,8]*x[:,9]]);names+=['lon_squared','lat_squared','lon_lat']
for yr in sorted(s.year.unique())[1:]:x=np.column_stack([x,(s.year==yr).to_numpy(float)]);names.append('year_'+str(yr))
for lc in sorted(s.LC2018_grid_mode.unique())[1:]:x=np.column_stack([x,(s.LC2018_grid_mode==lc).to_numpy(float)]);names.append('landcover_'+str(int(lc)))
rows=[];frags=[];diags=[]
for e in ['EVI','LAI','combined']:
 fit=sm.GLM(s[e].to_numpy(),x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':s.grid_id.to_numpy()},maxiter=100);V=fit.cov_params();b=fit.params;mu={};gr={}
 for ctx in ['permanent','seasonal','other']:
  for hh in [0,1]:
   xp=x.copy();xp[:,1]=hh;xp[:,2]=ctx=='seasonal';xp[:,3]=ctx=='other';xp[:,4]=hh*(ctx=='seasonal');xp[:,5]=hh*(ctx=='other');p=expit(xp@b);mu[ctx,hh]=p.mean();gr[ctx,hh]=(xp*(p*(1-p))[:,None]).mean(axis=0)
  est=mu[ctx,1]-mu[ctx,0];grad=gr[ctx,1]-gr[ctx,0];se=np.sqrt(max(0,grad@V@grad));sub=s[s.context==ctx];frags.append(dict(group=ctx,response=e,risk_low=mu[ctx,0],risk_high=mu[ctx,1],fragility=est,ci_low=est-1.96*se,ci_high=est+1.96*se,n_total_eligible=len(c),n_model=len(s),n_grid_model=s.grid_id.nunique(),n_context_model=len(sub),n_context_grid=sub.grid_id.nunique(),n_low=int((sub.HDSI<q.iloc[0]).sum()),n_high=int((sub.HDSI>q.iloc[1]).sum())))
 for ctx in ['seasonal','other']:
  est=mu[ctx,1]-mu[ctx,0]-mu['permanent',1]+mu['permanent',0];grad=gr[ctx,1]-gr[ctx,0]-gr['permanent',1]+gr['permanent',0];se=np.sqrt(max(0,grad@V@grad));m1=pd.read_csv(ROOT/out/'M1_AFB.csv');base=m1[(m1.comparison==ctx)&(m1.response==e)].AFB.iloc[0];rows.append(dict(model='M2_harmonized',response=e,comparison=ctx,AFB=est,ci_low=est-1.96*se,ci_high=est+1.96*se,M1_AFB=base,change_from_M1=est-base,attenuation_fraction=1-est/base,n_total=len(c),n_model=len(s),n_grid=s.grid_id.nunique(),n_low=int((s.HDSI<q.iloc[0]).sum()),n_high=int((s.HDSI>q.iloc[1]).sum()),ci_method='grid-cluster sandwich delta method'))
 save(pd.DataFrame({'term':names,'coefficient':fit.params,'se_grid_cluster':fit.bse}),out+'M2_coefficients_'+e+'.csv');diags.append(dict(response=e,converged=bool(fit.converged),n_parameters=len(b),rank=int(np.linalg.matrix_rank(x)),condition_number=float(np.linalg.cond(x))));print(e,rows[-2:],flush=True)
save(pd.DataFrame(rows),out+'M2_AFB.csv');save(pd.DataFrame(frags),out+'M2_fragility.csv');save(pd.DataFrame(diags),out+'model_diagnostics.csv')
(ROOT/out/'model_specification.json').write_text(json.dumps({'formula':'high * context + elevation + slope + lon + lat + lon^2 + lat^2 + lon*lat + annual_temperature + log1p(annual_precipitation) + C(2018_grid_modal_landcover) + C(year)','covariate_scaling':'continuous predictors centered and scaled in complete-case low/high fit sample','landcover':'new 2018 MCD12Q1.061 native categories aggregated by spatial mode, not original averaged category codes','climate':'WorldClim v1 1960-1990 independent long-term background; not event-anomaly reference','uncertainty':'grid-cluster sandwich conditional on empirical standardization distribution; cross-grid spatial dependence remains','residual_uncontrolled':['irrigation','groundwater','soil properties'],'distance_inputs':'original manuscript panel; corrected-distance sensitivity separate'},indent=2))
