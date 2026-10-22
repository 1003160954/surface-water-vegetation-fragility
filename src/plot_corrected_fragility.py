from pathlib import Path
import pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];O=ROOT/'diagnostic_figures';(O/'figures').mkdir(parents=True,exist_ok=True)
a=pd.read_csv(ROOT/'experiments/E01_spatial_support/corrected_1km/context_fragility.csv');fig,axs=plt.subplots(1,3,figsize=(12,5))
for ax,e in zip(axs,['EVI','LAI','combined']):
 t=a[a.response==e].set_index('group').loc[['permanent','seasonal','other']];ax.bar(range(3),t.fragility,color=['#3678a8','#dc9c22','#777777']);ax.errorbar(range(3),t.fragility,yerr=[t.fragility-t.ci_low,t.ci_high-t.fragility],fmt='none',ecolor='black',capsize=4);ax.set_xticks(range(3),['Permanent','Seasonal','Other'],rotation=20);ax.set_title(e);ax.set_ylabel('Risk(high HDSI) minus risk(low HDSI)')
fig.suptitle('Corrected aquatic contexts; paired grid-bootstrap95%intervals');fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(O/'figures'/('Fig5_context_fragility.'+ext),dpi=200,bbox_inches='tight')
plt.close(fig)

import shutil
for ext in ['png','pdf']:
 old=O/'figures'/('Fig5_distance_audit.'+ext)
 if old.exists():shutil.copy2(old,O/'figures'/('SFig18_distance_validation.'+ext))
mapping=O/'figure_source_map.csv'
if mapping.exists():
 t=pd.read_csv(mapping);mask=t.figure.eq('Fig5_distance_audit');t.loc[mask,'figure']='Fig5_context_fragility';t.loc[mask,'sources']='experiments/E01_spatial_support/corrected_1km/context_fragility.csv; source calculation audit retained as SFig18';t.to_csv(mapping,index=False)
