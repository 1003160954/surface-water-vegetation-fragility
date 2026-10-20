from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
raw,d=load_panel();out='experiments/E00_original_reproduction/'
# Test whether the pre-existing band labels actually describe any-water rather than permanent-water distance.
lo=d.water_dist_band.astype(str).str.extract(r'^([0-9.]+)-')[0].astype(float);hi=d.water_dist_band.astype(str).str.extract(r'-([0-9.]+)km')[0].astype(float);rows=[]
for col in ['dist_perm_water_km','dist_seasonal_water_km','dist_any_water_km']:
 matches=(d[col]>=lo)&(d[col]<=hi);rows.append(dict(numeric_distance=col,n_matching=int(matches.sum()),n_total=len(d),fraction_matching=matches.mean()))
save(pd.DataFrame(rows),out+'stored_band_numeric_match.csv')
h=pd.read_csv(ROOT/'docs/historical/water_context_HDSI_risk_with_CI.csv');r=pd.read_csv(ROOT/out/'context_fragility.csv');rs=[]
for row in r.itertuples():
 ctx={'permanent':'near_permanent_lt1km','seasonal':'near_seasonal_lt1km','other':'far_or_weak_water'}[row.group];resp='Joint' if row.response=='combined' else row.response
 for lev in ['low','high']:
  old=h[(h.water_context==ctx)&(h.response==resp)&(h.HDSI_group==lev)].iloc[0];rs.append(dict(context=row.group,response=row.response,exposure=lev,recomputed_risk=getattr(row,'risk_'+lev),historical_risk=old.risk,risk_difference=getattr(row,'risk_'+lev)-old.risk,recomputed_n=getattr(row,'n_'+lev),historical_n=old.rows))
save(pd.DataFrame(rs),out+'risk_exact_comparison.csv')
print(pd.DataFrame(rows).to_string(index=False));print('max risk diff',pd.DataFrame(rs).risk_difference.abs().max())
