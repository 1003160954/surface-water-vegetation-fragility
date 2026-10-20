from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
c=pd.read_parquet(__import__('os').environ.get('AQUATIC_CLIMATOLOGY',str(ROOT/'data/raw/climatology.parquet')));rows=[]
for v in ['EVI','LAI','VPD','LST_Day','SM_rootproxy','Precip','Tair']:
 n=c[v+'_clim_count'];sd=c[v+'_clim_std'];
 for k in sorted(n.unique()):rows.append(dict(variable=v,n_reference=int(k),n_grid_step=int((n==k).sum()),fraction=float((n==k).mean())))
 rows.append(dict(variable=v,n_reference='SUMMARY',n_grid_step=len(c),fraction=np.nan,std_missing=int(sd.isna().sum()),std_zero=int(sd.eq(0).sum()),std_lt_1e_8=int((sd.notna()&(sd<1e-8)).sum()),std_min_positive=float(sd[sd>0].min())))
save(pd.DataFrame(rows),'experiments/E03_independent_baseline/existing_climatology_audit.csv')
(ROOT/'experiments/E03_independent_baseline/status.json').write_text(json.dumps({'status':'blocked_independent_baseline','rows_grid_step':len(c),'unique_grids':c.grid_id.nunique(),'steps':sorted(map(int,c.step16.unique())),'count_max':int(c.EVI_clim_count.max()),'reference_year_metadata_present':False,'conclusion':'Existing reference counts at most six; no evidence of independent 2001-2017 reference. Must not relabel as long-term climatology.','original_product_versions_and_season_windows':'original product versions unverified; season windows recovered from upstream notebook cell28','sha256':hashlib.file_digest(open(__import__('os').environ.get('AQUATIC_CLIMATOLOGY',str(ROOT/'data/raw/climatology.parquet')),'rb'),'sha256').hexdigest()},indent=2))
print(pd.DataFrame(rows).to_string(index=False))
