from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

O=Path(__file__).resolve().parent
a=pd.read_csv(O/'independent_reference_AFB.csv')
names={'main_reference_paired':'Recent reference, paired sample','independent_response_fixed_exposure':'Historical response reference','independent_response_and_exposure':'Historical response and exposure references'}
colors=['#68737d','#238b8e','#b75c40']
fig,axes=plt.subplots(2,3,figsize=(10.8,6.2),sharex=True,sharey=True)
for row,minimum in enumerate([10,15]):
 for col,event in enumerate(['EVI','LAI','combined']):
  ax=axes[row,col]
  for k,(scenario,label) in enumerate(names.items()):
   for j,comp in enumerate(['seasonal','other']):
    r=a[(a.min_reference_years==minimum)&(a.response==event)&(a.scenario==scenario)&(a.comparison==comp)].iloc[0]
    y=j+(k-1)*.16
    ax.errorbar(r.AFB,y,xerr=[[r.AFB-r.ci_low],[r.ci_high-r.AFB]],fmt='o',markersize=4,capsize=3,color=colors[k],label=label if j==0 else None)
  ax.axvline(0,color='#555',lw=.8,ls='--');ax.set_yticks([0,1],['Seasonal − permanent','Other − permanent']);ax.set_ylim(-.5,1.5);ax.grid(axis='x',alpha=.18)
  ax.set_title(f'{chr(65+row*3+col)}  '+('Combined' if event=='combined' else event)+f' | ≥{minimum} reference years',fontsize=10)
  ax.spines[['top','right']].set_visible(False)
  if row==1:ax.set_xlabel('Fragility contrast (AFB)')
handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=1,frameon=False,fontsize=10)
fig.tight_layout(rect=[0,.15,1,1]);fig.savefig(O/'Fig_S18.png',dpi=300,bbox_inches='tight');fig.savefig(O/'Fig_S18.pdf',bbox_inches='tight');plt.close(fig)
