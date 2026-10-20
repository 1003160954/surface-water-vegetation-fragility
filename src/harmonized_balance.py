from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
_,d=load_panel();c=pd.read_csv(ROOT/'experiments/E02_confounding/harmonized_landcover_climate_covariates.csv');g=d.drop_duplicates('grid_id').merge(c,on='grid_id',validate='one_to_one');rows=[]
for ctx in ['seasonal','other']:
 a=g[g.context=='permanent'];b=g[g.context==ctx]
 for col in ['elevation','slope','climate_annual_temp_C','climate_annual_precip_mm']:
  sd=np.sqrt((a[col].var()+b[col].var())/2);rows.append(dict(comparison=ctx,covariate=col,permanent_mean=a[col].mean(),comparison_mean=b[col].mean(),SMD=(b[col].mean()-a[col].mean())/sd,n_permanent=len(a),n_comparison=len(b),permanent_p05=a[col].quantile(.05),permanent_p95=a[col].quantile(.95),comparison_p05=b[col].quantile(.05),comparison_p95=b[col].quantile(.95)))
 for lc in sorted(g.LC2018_grid_mode.unique()):
  p=(a.LC2018_grid_mode==lc).mean();q=(b.LC2018_grid_mode==lc).mean();sd=np.sqrt((p*(1-p)+q*(1-q))/2);rows.append(dict(comparison=ctx,covariate='LC2018_class_'+str(int(lc)),permanent_mean=p,comparison_mean=q,SMD=(q-p)/sd if sd else np.nan,n_permanent=len(a),n_comparison=len(b)))
save(pd.DataFrame(rows),'experiments/E02_confounding/harmonized/covariate_balance.csv')
