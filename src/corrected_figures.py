from pathlib import Path
import sys,pandas as pd,numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False});colors={'EVI':'#25805c','LAI':'#7857a3','combined':'#b63d49'};out=ROOT/'experiments/E01_spatial_support/corrected_1km';mapping=pd.read_csv(ROOT/'experiments/E05_figure_audit/figure_source_map.csv');new=[]
def finish(fig,name,source):
 fig.savefig(ROOT/'figures'/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(ROOT/'figures'/(name+'.pdf'),bbox_inches='tight');plt.close(fig);new.append(dict(figure=name,source=source,status='executed corrected physical1km centroid contexts; within-cell national analysis paused'))
a=pd.read_csv(ROOT/'experiments/E00_original_reproduction/AFB.csv');b=pd.read_csv(out/'AFB.csv');c=pd.read_csv(out/'harmonized/M2_AFB.csv');e=pd.read_csv(out/'harmonized/common_landscape_AFB.csv');fig,axs=plt.subplots(1,3,figsize=(14,5))
for ax,response in zip(axs,colors):
 rows=[t[(t.response==response)&(t.comparison=='seasonal')].iloc[0] for t in [a,b,c,e]]
 for i,r in enumerate(rows):ax.errorbar(r.AFB,i,xerr=[[r.AFB-r.ci_low],[r.ci_high-r.AFB]],fmt='o',color=colors[response],capsize=3)
 ax.set_yticks(range(4),['Original contexts / M0','Corrected1km / M0','Corrected1km / M2','Corrected / landscape']);ax.invert_yaxis();ax.set_ylim(3.5,-.5);ax.axvline(0,color='grey',lw=.8);ax.set_title(response);ax.set_xlabel('AFB: seasonal minus permanent')
fig.suptitle('Physical-distance correction and comparable landscapes\nOriginal and corrected contexts contain different grids; all contrasts carry95% intervals');fig.tight_layout();finish(fig,'SFig9_corrected_1km_permanence','Original E00/AFB.csv; corrected1km/AFB.csv, harmonized/M2_AFB.csv and common_landscape_AFB.csv')
b=pd.read_csv(out/'robustness/spatial_block_AFB.csv');fig,axs=plt.subplots(1,2,figsize=(12,6))
for ax,context in zip(axs,['seasonal','other']):
 t=b[b.comparison==context]
 for i,r in enumerate(t.itertuples()):ax.errorbar(r.AFB,i,xerr=[[r.AFB-r.ci_low],[r.ci_high-r.AFB]],fmt='o',color=colors[r.response],capsize=3)
 ax.set_yticks(range(len(t)),t.response+' / '+t.block_degrees.astype(str)+' deg');ax.axvline(0,color='grey',lw=.8);ax.set_xlabel('AFB: comparison minus permanent');ax.set_title(context)
fig.suptitle('Corrected physical1km contexts: spatial-block uncertainty');fig.tight_layout();finish(fig,'SFig10_corrected_1km_spatial_blocks','corrected1km/robustness/spatial_block_AFB.csv;300 draws seed20261020; origin0,0')
_,d=load_panel();g=d.drop_duplicates('grid_id')[['grid_id','lon','lat']].merge(pd.read_csv(ROOT/'experiments/E01_spatial_support/national_corrected_1km_contexts.csv'),on='grid_id',validate='one_to_one');g['context']=np.where(g.near_perm_1km_corrected==1,'permanent',np.where(g.near_seasonal_1km_corrected==1,'seasonal','other'));g.to_csv(out/'corrected_context_map.csv',index=False);cc=['#3678a8','#dc9c22','#777777'];categories=['permanent','seasonal','other'];fig,ax=plt.subplots(figsize=(8,7))
for cat,color in reversed(list(zip(categories,cc))):
 t=g[g.context==cat];ax.scatter(t.lon,t.lat,color=color,s=2,label=cat+' (G='+str(len(t))+')')
ax.legend(loc='lower left',fontsize=9,markerscale=3);ax.set_xlabel('Longitude');ax.set_ylabel('Latitude');ax.set_title('Corrected1km centroid aquatic contexts\nJRC occurrence>=80%, seasonality>=10 months;30m UTM');fig.text(.5,.01,'Source: EC JRC/Google; response grids remain0.1 degrees',ha='center',fontsize=9);fig.tight_layout(rect=[0,.03,1,1]);finish(fig,'SFig11_corrected_1km_context_map','national_corrected_1km_contexts.csv; main analysis grid coordinates; corrected_context_map.csv')
save(pd.concat([mapping[~mapping.figure.isin([n['figure'] for n in new])],pd.DataFrame(new)],ignore_index=True),'experiments/E05_figure_audit/figure_source_map.csv')
print('Corrected figures complete',len(new))
