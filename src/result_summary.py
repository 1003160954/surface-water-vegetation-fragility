import pandas as pd
from pathlib import Path
p=Path(__file__).resolve().parents[1]/'experiments/E04_robustness'
for f,cols in [('HDSI_variant_AFB.csv',['scenario','response','AFB','ci_low','ci_high']),('leave_analysis_year_out_AFB.csv',['excluded_year','response','AFB','ci_low','ci_high']),('event_threshold_AFB.csv',['threshold','response','AFB','ci_low','ci_high'])]:
 a=pd.read_csv(p/f);print(f);print(a[a.comparison=='seasonal'][cols].to_string(index=False))
a=pd.read_csv(p.parent/'E00_original_reproduction/fine_distance_fragility.csv');print(a[a.response=='combined'][['group','n_grid','fragility','ci_low','ci_high']].to_string(index=False))
