from pathlib import Path
import argparse,json,hashlib
import pandas as pd,numpy as np
VARS=['EVI','LAI','VPD','LST_Day','SM_rootproxy','Precip','Tair']
def build_reference(history,start=2001,end=2017,min_count=10,min_std=1e-6):
 if not {'grid_id','year','step16',*VARS}.issubset(history.columns):raise ValueError('Missing raw 16-day fields')
 if history.duplicated(['grid_id','year','step16']).any():raise ValueError('Duplicate historical grid-year-step')
 h=history[history.year.between(start,end)].copy()
 if h.empty:raise ValueError('No independent historical observations')
 if end>=2018:raise ValueError('Reference overlaps 2018-2023 analysis')
 # Same conservative raw cleaning as recovered upstream code.
 h.loc[~h.EVI.between(-.2,1.2),'EVI']=np.nan;h.loc[~h.LAI.between(0,12),'LAI']=np.nan
 stats=h.groupby(['grid_id','step16'])[VARS].agg(['mean','std','count']);stats.columns=[v+'_ref_'+s for v,s in stats.columns];stats=stats.reset_index()
 for v in VARS:
  bad=(stats[v+'_ref_count']<min_count)|(stats[v+'_ref_std']<=min_std);stats.loc[bad,v+'_ref_std']=np.nan
 return stats

def apply_reference(current,reference):
 if current.duplicated(['grid_id','year','step16']).any():raise ValueError('Duplicate current grid-year-step')
 if not current.year.between(2018,2023).all():raise ValueError('Unexpected analysis years')
 current=current.copy()
 current.loc[~current.EVI.between(-.2,1.2),'EVI']=np.nan
 current.loc[~current.LAI.between(0,12),'LAI']=np.nan
 x=current.merge(reference,on=['grid_id','step16'],how='left',validate='many_to_one')
 for v in VARS:x[v+'_z_independent']=((x[v]-x[v+'_ref_mean'])/x[v+'_ref_std']).clip(-5,5)
 return x

def main():
 p=argparse.ArgumentParser();p.add_argument('--history',required=True);p.add_argument('--current',required=True);p.add_argument('--provenance',required=True);p.add_argument('--output',required=True);p.add_argument('--min-count',type=int,default=10);a=p.parse_args();prov=json.loads(Path(a.provenance).read_text(encoding='utf8'))
 for key in ['products_and_versions','qa_rules','spatial_transform','season_windows','step_definition','units']:
  if key not in prov or not prov[key]:raise ValueError('Verified provenance required: '+key)
 h=pd.read_parquet(a.history);cur=pd.read_parquet(a.current);ref=build_reference(h,min_count=a.min_count);out=Path(a.output);out.mkdir(parents=True,exist_ok=True);ref.to_parquet(out/'reference_2001_2017.parquet',index=False);apply_reference(cur,ref).to_parquet(out/'analysis_2018_2023_independent.parquet',index=False)
 (out/'reference_config.json').write_text(json.dumps(dict(prov,reference_start=2001,reference_end=2017,min_count=a.min_count,min_std=1e-6,sd_ddof=1,z_clip=[-5,5]),indent=2))
if __name__=='__main__':main()
