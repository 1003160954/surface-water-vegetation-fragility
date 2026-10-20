from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
import statsmodels.api as sm
from scipy.special import expit
from patsy import dmatrix,build_design_matrices
_,d=load_panel();q=d.HDSI.quantile([.25,.75]);s=d[(d.HDSI<q.iloc[0])|(d.HDSI>q.iloc[1])].copy();h=(s.HDSI>q.iloc[1]).to_numpy(float);z=np.log1p(s.dist_perm_water_km.to_numpy());upper=float(np.log1p(d.dist_perm_water_km.max()))
B=dmatrix('0 + bs(z, df=5, degree=3, include_intercept=False, lower_bound=0, upper_bound=upper)',{'z':z,'upper':upper});info=B.design_info;basis=np.asarray(B);x=np.column_stack([np.ones(len(s)),h,basis,basis*h[:,None]])
for col in ['elevation','slope','lon','lat']:x=np.column_stack([x,(s[col].to_numpy()-s[col].mean())/s[col].std()])
x=np.column_stack([x,x[:,-2]**2,x[:,-1]**2,x[:,-2]*x[:,-1]])
for yr in sorted(s.year.unique())[1:]:x=np.column_stack([x,(s.year==yr).to_numpy(float)])
rows=[];meta=[]
for event in ['EVI','LAI','combined']:
 fit=sm.GLM(s[event].to_numpy(),x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':s.grid_id.to_numpy()},maxiter=100);v=fit.cov_params();b=fit.params
 for dist in np.linspace(0,15,31):
  bb=np.asarray(build_design_matrices([info],{'z':np.array([np.log1p(dist)]),'upper':upper})[0])[0];means=[];grads=[]
  for hh in [0,1]:
   xp=x.copy();xp[:,1]=hh;xp[:,2:7]=bb;xp[:,7:12]=bb*hh;p=expit(xp@b);means.append(p.mean());grads.append((xp*(p*(1-p))[:,None]).mean(axis=0))
  g=grads[1]-grads[0];se=np.sqrt(max(0,g@v@g));e=means[1]-means[0];near=d[abs(d.dist_perm_water_km-dist)<=.5]
  rows.append(dict(response=event,distance_km=dist,risk_low=means[0],risk_high=means[1],fragility=e,ci_low=e-1.96*se,ci_high=e+1.96*se,n_model=len(s),n_grid_model=s.grid_id.nunique(),n_grid_within_0p5km=near.grid_id.nunique(),n_rows_within_0p5km=len(near),ci_method='grid-cluster delta method; pointwise'))
 meta.append(dict(response=event,converged=bool(fit.converged),n_parameters=len(b),n_model=len(s),rank=int(np.linalg.matrix_rank(x))))
 print('spline',event,flush=True)
save(pd.DataFrame(rows),'experiments/E04_robustness/spline_predictions.csv');save(pd.DataFrame(meta),'experiments/E04_robustness/spline_diagnostics.csv')
(ROOT/'experiments/E04_robustness/spline_config.json').write_text(json.dumps({'basis':'cubic B-spline on log1p distance, 5 df, no basis intercept; intercept and high indicator separate','patsy_design':str(info),'formula':'high * bs(log1p(distance),df=5) + elevation + slope + lon + lat + lon^2 + lat^2 + lon*lat + year fixed effects','standardization':'full observed low/high sample common covariate distribution','CI':'95% pointwise grid-cluster sandwich delta method','landcover':'omitted invalid fractional codes','interpretation':'descriptive centroid-context association; no ecological reach or causal mechanism inferred'},indent=2))
