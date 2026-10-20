from pathlib import Path
import sys,json,datetime
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
import statsmodels.api as sm
from scipy.special import expit
raw,d=load_panel(); out='experiments/E02_confounding/'
q=d.HDSI.quantile([.25,.75]);cov=['elevation','slope','lon','lat','year']
cc=d[cov].notna().all(axis=1); c=d.loc[cc].copy()
for name,s in [('M0',d),('M1',c)]:
 r,a=group_effects(s);r['model']=name;a['model']=name;save(r,out+name+'_fragility.csv');save(a,out+name+'_AFB.csv')
s=c[(c.HDSI<q.iloc[0])|(c.HDSI>q.iloc[1])].copy(); hi=(s.HDSI>q.iloc[1]).to_numpy(float);ss=s.context.eq('seasonal').to_numpy(float);oth=s.context.eq('other').to_numpy(float)
x=np.column_stack([np.ones(len(s)),hi,ss,oth,hi*ss,hi*oth]);names=['Intercept','high','seasonal','other','high_seasonal','high_other']
for col in ['elevation','slope','lon','lat']:
 z=(s[col].to_numpy()-s[col].mean())/s[col].std();x=np.column_stack([x,z]);names.append(col+'_scaled')
x=np.column_stack([x,x[:,8]**2,x[:,9]**2,x[:,8]*x[:,9]]);names+=['lon_squared','lat_squared','lon_lat']
for yr in sorted(s.year.unique())[1:]:x=np.column_stack([x,(s.year==yr).to_numpy(float)]);names.append('year_'+str(yr))
risks=[];afbs=[];diagnostics=[]
for event in ['EVI','LAI','combined']:
 fit=sm.GLM(s[event].to_numpy(),x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':s.grid_id.to_numpy()},maxiter=100)
 V=fit.cov_params();b=fit.params;means={};grads={}
 save(pd.DataFrame({'term':names,'coefficient':b,'se_grid_cluster':fit.bse}),out+'M2_coefficients_'+event+'.csv')
 for ctx in ['permanent','seasonal','other']:
  for h in [0,1]:
   xp=x.copy();xp[:,1]=h;xp[:,2]=ctx=='seasonal';xp[:,3]=ctx=='other';xp[:,4]=h*(ctx=='seasonal');xp[:,5]=h*(ctx=='other')
   p=expit(xp@b);mu=p.mean();g=(xp*(p*(1-p))[:,None]).mean(axis=0);means[ctx,h]=mu;grads[ctx,h]=g
  est=means[ctx,1]-means[ctx,0];g=grads[ctx,1]-grads[ctx,0];se=float(np.sqrt(max(0,g@V@g)))
  sub=s[s.context==ctx]
  risks.append(dict(model='M2_partial',response=event,group=ctx,risk_low=means[ctx,0],risk_high=means[ctx,1],fragility=est,ci_low=est-1.96*se,ci_high=est+1.96*se,n_total_eligible=len(c),n_model=len(s),n_grid=s.grid_id.nunique(),n_context_model=len(sub),n_context_grid=sub.grid_id.nunique(),n_low=int((sub.HDSI<q.iloc[0]).sum()),n_high=int((sub.HDSI>q.iloc[1]).sum()),ci_method='grid-cluster sandwich and delta method'))
 for ctx in ['seasonal','other']:
  est=means[ctx,1]-means[ctx,0]-means['permanent',1]+means['permanent',0];g=grads[ctx,1]-grads[ctx,0]-grads['permanent',1]+grads['permanent',0];se=float(np.sqrt(max(0,g@V@g)))
  base=pd.read_csv(ROOT/out/'M1_AFB.csv');base=base[(base.response==event)&(base.comparison==ctx)].AFB.iloc[0]
  afbs.append(dict(model='M2_partial',response=event,comparison=ctx,AFB=est,ci_low=est-1.96*se,ci_high=est+1.96*se,M1_AFB=base,change_from_M1=est-base,attenuation_fraction=1-est/base,n_total=len(c),n_model=len(s),n_grid=s.grid_id.nunique(),n_low=int((s.HDSI<q.iloc[0]).sum()),n_high=int((s.HDSI>q.iloc[1]).sum()),ci_method='grid-cluster sandwich and delta method'))
 diagnostics.append(dict(response=event,converged=bool(fit.converged),rank=int(np.linalg.matrix_rank(x)),n_columns=x.shape[1],n_model=len(s),condition_number=float(np.linalg.cond(x))))
 print(event,afbs[-2:],flush=True)
save(pd.DataFrame(risks),out+'M2_fragility.csv');save(pd.DataFrame(afbs),out+'M2_AFB.csv');save(pd.DataFrame(diagnostics),out+'model_diagnostics.csv')
# Diagnose covariate support on grid summaries rather than repeated grid-years.
grid=c.groupby(['grid_id','context'],observed=True)[['elevation','slope','lon','lat','LC_Type1']].mean().reset_index();balance=[]
for ctx in ['seasonal','other']:
 for col in ['elevation','slope','lon','lat','LC_Type1']:
  a=grid.loc[grid.context=='permanent',col];b=grid.loc[grid.context==ctx,col];sd=np.sqrt((a.var()+b.var())/2)
  balance.append(dict(comparison=ctx,covariate=col,permanent_mean=a.mean(),comparison_mean=b.mean(),SMD=(b.mean()-a.mean())/sd,n_permanent=len(a),n_comparison=len(b),permanent_p05=a.quantile(.05),permanent_p95=a.quantile(.95),comparison_p05=b.quantile(.05),comparison_p95=b.quantile(.95),caution='averaged categorical codes are invalid landcover control' if col=='LC_Type1' else ''))
save(pd.DataFrame(balance),out+'covariate_balance.csv')
(ROOT/out/'model_specification.json').write_text(json.dumps({'model':'logistic high * context + scaled elevation + scaled slope + lon + lat + lon^2 + lat^2 + lon*lat + C(year)','reference':'permanent, low HDSI, 2018','standardization':'common complete-case low/high sample, observed covariate distribution','uncertainty':'grid clustered sandwich, delta method; does not resolve spatial cross-grid dependence','landcover':'omitted: 41.95% fractional LC_Type1 codes; original categorical landcover cannot be reconstructed','uncontrolled':['irrigation','groundwater','soil properties','independent long-term climate','valid categorical landcover'],'M0_M1_identical':bool(cc.all()),'seed':SEED},indent=2))
