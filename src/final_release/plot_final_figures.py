from pathlib import Path
import os,shutil
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch
ROOT=Path(__file__).resolve().parents[2];os.chdir(ROOT)
out=Path('final_results');F=out/'figures';f=F;data=out/'data';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Times New Roman','font.size':11,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(F/(name+'.'+ext),dpi=250,bbox_inches='tight')
 plt.close(fig)
E=data/'maps';P=data/'context'
g=pd.read_csv(E/'grid_map_data.csv');ctx=pd.read_csv(P/'corrected_context_map.csv')
colors=['#eeeeee','#1747a6','#25854e','#cc4444','#67b6d3','#a5b930','#d29baf','#e17d25','#c9ad24','#555555','#aaaaaa']
names=['0 Non-water / zero-filled','1 Permanent','2 New permanent','3 Lost permanent','4 Seasonal','5 New seasonal','6 Lost seasonal','7 Seasonal to permanent','8 Permanent to seasonal','9 Ephemeral permanent','10 Ephemeral seasonal']
fig,axs=plt.subplots(1,3,figsize=(14,5.4))
for ax,col,title in zip(axs,['water_occurrence','water_seasonality','water_transition'],['a  Water occurrence','b  Water seasonality','c  Water transition']):
 if col=='water_transition':ax.scatter(g.lon,g.lat,c=g[col],s=1,cmap=ListedColormap(colors),norm=BoundaryNorm(np.arange(-.5,11.5),11))
 else:
  s=ax.scatter(g.lon,g.lat,c=g[col],s=1,cmap='viridis');fig.colorbar(s,ax=ax,shrink=.7,label='%' if col=='water_occurrence' else 'Months')
 ax.set(xlabel='Longitude (°E)',ylabel='Latitude (°N)',title=title)
fig.legend(handles=[Patch(facecolor=c,label=n) for c,n in zip(colors,names)],loc='lower center',ncol=4,fontsize=10);fig.tight_layout(rect=(0,.2,1,1));save(fig,'Fig2')
g=g.merge(ctx[['grid_id','context']],on='grid_id',suffixes=('_old','')) if 'context' not in g else g.drop(columns='context').merge(ctx[['grid_id','context']],on='grid_id')
cc={'permanent':'#3678a8','seasonal':'#dc9c22','other':'#777777'}
fig,axs=plt.subplots(1,2,figsize=(11,5))
for c,col in cc.items():
 z=g[g.context==c];axs[0].scatter(z.lon,z.lat,s=2,c=col,label=f'{c.title()} (n = {len(z):,})')
axs[0].legend(fontsize=9,loc='lower left');axs[0].set_title('a  Aquatic context')
sc=axs[1].scatter(g.lon,g.lat,c=g.combined,s=2,cmap='YlOrRd',vmin=0,vmax=1);fig.colorbar(sc,ax=axs[1],label='Fraction of valid years');axs[1].set_title('b  Combined-event frequency')
for ax in axs:ax.set(xlabel='Longitude (°E)',ylabel='Latitude (°N)')
fig.tight_layout();save(fig,'Fig3')
r=pd.read_csv(E/'exposure_risks.csv');m=pd.read_csv(E/'high_HDSI_map.csv');fig,axs=plt.subplots(2,2,figsize=(10,7))
sc=axs[0,0].scatter(m.lon,m.lat,c=m.high_HDSI_fraction,s=2,cmap='YlOrRd');fig.colorbar(sc,ax=axs[0,0],label='Fraction of valid years');axs[0,0].set(title='a  High-HDSI exposure',xlabel='Longitude (°E)',ylabel='Latitude (°N)')
for ax,ev,col,lab in zip(list(axs.flat)[1:],['EVI','LAI','combined'],['#25805c','#7857a3','#b63d49'],['b','c','d']):
 z=r[r.response==ev].set_index('group').loc[['low','mid','high']];ax.bar(range(3),z.risk,color=col);ax.errorbar(range(3),z.risk,yerr=[z.risk-z.ci_low,z.ci_high-z.risk],fmt='none',ecolor='black',capsize=3);ax.set(xticks=range(3),xticklabels=['Low','Middle','High'],ylabel='Event probability',xlabel='HDSI group',title=lab+'  '+('Combined events' if ev=='combined' else ev))
fig.tight_layout();save(fig,'Fig4')
r=pd.read_csv(P/'context_fragility.csv');fig,axs=plt.subplots(2,3,figsize=(12,7))
# High-exposure intervals are taken from the separately computed risk table.
for j,ev in enumerate(['EVI','LAI','combined']):
 z=r[r.response==ev].set_index('group').loc[list(cc)]
 for i,c in enumerate(cc):
  axs[0,j].bar(i,z.loc[c,'risk_high'],color=cc[c]);axs[1,j].bar(i,z.loc[c,'fragility'],color=cc[c])
 axs[1,j].errorbar(range(3),z.fragility,yerr=[z.fragility-z.ci_low,z.ci_high-z.fragility],fmt='none',ecolor='black',capsize=3)
 for row in range(2):axs[row,j].set(xticks=range(3),xticklabels=['Permanent','Seasonal','Other'],title='abcdef'[row*3+j]+'  '+('Combined events' if ev=='combined' else ev),ylabel='High-HDSI event probability' if row==0 else 'Fragility')

high=pd.read_csv(data/'Fig6_high_risks.csv')
for j,ev in enumerate(['EVI','LAI','combined']):
 for i,c in enumerate(cc):
  z=high[(high.response==ev)&(high.context==c)].iloc[0]
  axs[0,j].errorbar(i,z.risk,yerr=[[z.risk-z.ci_low],[z.ci_high-z.risk]],fmt='none',ecolor='black',capsize=3)
fig.tight_layout();save(fig,'Fig6')

means=pd.read_csv(data/'diagnostics/corrected_component_and_seasonal_means.csv')
colors={'permanent':'#2878a5','seasonal':'#dc9f24','other':'#777777'}
f,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for j,(ctx,color) in enumerate(colors.items()):
 sub=means[(means.metric!='LAI_seasonal_change')&means.context.eq(ctx)];x=np.arange(5)+(j-1)*.22;axes[0,0].errorbar(x,sub['mean'],yerr=[sub['mean']-sub.ci_low,sub.ci_high-sub['mean']],fmt='o',color=color,label=ctx,capsize=3)
 sub=means[means.metric.eq('LAI_seasonal_change')&means.context.eq(ctx)];axes[0,1].errorbar(np.arange(2)+(j-1)*.15,sub['mean'],yerr=[sub['mean']-sub.ci_low,sub.ci_high-sub['mean']],fmt='o-',color=color,label=ctx,capsize=3)
axes[0,0].set(xticks=range(5),xticklabels=['VPD','LST','SM deficit','P deficit','Tair'],ylabel='Mean standardized component',title='a  Components under high HDSI');axes[0,0].legend(fontsize=8)
axes[0,1].set(xticks=[0,1],xticklabels=['Low HDSI','High HDSI'],ylabel='Late-summer minus spring LAI z',title='b  Seasonal LAI change');axes[0,1].axhline(0,color='grey',lw=.7)
risks=pd.read_csv(data/'context/context_fragility.csv')
for ax,event,label in zip(axes[1],['LAI','combined'],['c  LAI event risks','d  Combined-event risks']):
 for j,(ctx,color) in enumerate(colors.items()):
  z=risks.query('group==@ctx and response==@event').iloc[0];ax.plot([z.risk_low,z.risk_high],[j,j],'o-',color=color,label=ctx)
 ax.set(yticks=range(3),yticklabels=list(colors),xlabel='Negative-anomaly event probability',title=label)
for ext in ['png','pdf']:f.savefig(F/f'Fig7.{ext}',dpi=200)
plt.close(f)

from pathlib import Path
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('final_results/figures');out.mkdir(exist_ok=True)
d=pd.read_csv('final_results/data/Figure5_data.csv');events=['EVI','LAI','combined'];counts=d[d.response.eq('EVI')].sort_values('bin_index').n_grid.astype(int).tolist();assert sum(counts)==17041
for e in events:assert d[d.response.eq(e)].sort_values('bin_index').n_grid.astype(int).tolist()==counts
d.to_csv(out/'Figure5_data.csv',index=False)
plt.rcParams.update({'font.family':'Times New Roman','font.size':12,'axes.labelsize':13,'axes.titlesize':15,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
f,axs=plt.subplots(3,1,figsize=(8,9),sharey=True)
intervals=['[0, 1]','(1, 5]','(5, 10]','(10, 20]','(20, 40]','(40, 61.44]']
for ax,e,title,color in zip(axs,events,['(a) EVI','(b) LAI','(c) Combined events'],['#267895','#568A64','#B67B37']):
 s=d[d.response.eq(e)].sort_values('bin_index');x=np.arange(6);ax.errorbar(x,s.fragility,yerr=[s.fragility-s.ci_low,s.ci_high-s.fragility],fmt='o-',lw=1.6,markersize=6,color=color,capsize=4);ax.set_xticks(x,[f'{lab}\nn = {n:,}' for lab,n in zip(intervals,counts)],fontsize=10);ax.set_title(title,loc='left');ax.set_xlabel('Centroid distance interval (km)',labelpad=10);ax.set_ylim(0,.52);ax.set_yticks(np.arange(0,.51,.1));ax.set_xlim(-.45,5.6);ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)
axs[0].set_ylabel('Vegetation fragility\n(high-HDSI risk − low-HDSI risk)')
f.subplots_adjust(left=.14,right=.98,top=.94,bottom=.1,hspace=.62)
f.text(.53,.98,'Distance-analysis subset: 17,041 unique grid cells · 2018–2023',ha='center',fontsize=15)
for ext in ['png','pdf']:f.savefig(out/f'Fig5.{ext}',dpi=400,facecolor='white')
plt.close(f)

from pathlib import Path
import pandas as pd,numpy as np,shutil
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
src=Path('final_results/data/diagnostics');out=Path('final_results/figures');out.mkdir(exist_ok=True)
a=pd.read_csv(src/'corrected_regional_AFB.csv');d=a[a.comparison.eq('other')].copy();order=['Northeast','North_China','East_China','South_China','Southwest','Northwest_Arid','Tibetan_Plateau','Other'];names=['Northeast China','North China','East China','South China','Southwest China','Northwest arid region','Tibetan Plateau','Other transition regions'];assert len(d)==24 and d.supported.all() and d.valid_draws.eq(300).all();d.to_csv(out/'Figure8_data.csv',index=False)
shutil.copy2(src/'corrected_regional_support.csv',out/'Figure8_sample_support.csv');shutil.copy2(src/'diagnostic_config.json',out/'analysis_config.json')
plt.rcParams.update({'font.family':'Times New Roman','font.size':12,'axes.titlesize':15,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
f,axs=plt.subplots(1,3,figsize=(13,6),sharey=True,sharex=True)
for ax,event,title,color in zip(axs,['EVI','LAI','combined'],['(a) EVI','(b) LAI','(c) Combined events'],['#267895','#568A64','#B67B37']):
 s=d[d.response.eq(event)].set_index('region').loc[order];ax.errorbar(s.AFB,np.arange(8),xerr=[s.AFB-s.ci_low,s.ci_high-s.AFB],fmt='o',color=color,markersize=6,capsize=4,elinewidth=1.3);ax.axvline(0,color='#666666',ls='--',lw=1);ax.set_yticks(np.arange(8),names);ax.set_title(title,loc='left');ax.set_xlabel('AFB (other − permanent)');ax.set_xlim(-.15,.30);ax.set_xticks([-.1,0,.1,.2,.3]);ax.grid(axis='x',alpha=.18);ax.set_axisbelow(True)
axs[0].invert_yaxis();f.subplots_adjust(left=.19,right=.985,bottom=.16,top=.87,wspace=.14);f.text(.59,.95,'Regional vegetation fragility contrasts',ha='center',fontsize=17);f.text(.59,.035,'Nationwide aquatic contexts; 95% grid-cluster bootstrap intervals.',ha='center',fontsize=11)
for ext in ['png','pdf']:f.savefig(out/f'Fig8.{ext}',dpi=400,facecolor='white')
plt.close(f)

out=Path('final_results');F=out/'figures';fig=F;data=out/'data';base=data/'distance'
b=pd.read_csv(base/'distance_band_fragility.csv');b=b[b.metric.eq('dist_perm_corrected_km')];colors=['#267895','#568A64','#B67B37'];events=['EVI','LAI','combined']
f,axs=plt.subplots(2,3,figsize=(13,7.5),sharey=True,layout='constrained')
for row,scheme in enumerate(['fine','grid_quintiles']):
 for ax,e,color in zip(axs[row],events,colors):
  s=b[(b.scheme==scheme)&(b.response==e)].sort_values('bin_index');x=np.arange(len(s));ax.errorbar(x,s.fragility,yerr=[s.fragility-s.ci_low,s.ci_high-s.fragility],fmt='o-',color=color,capsize=3,ms=4);labels=s.group.tolist() if row==0 else ['Q1','Q2','Q3','Q4','Q5'];ax.set_xticks(x,labels,rotation=65 if row==0 else 0,ha='right' if row==0 else 'center');ax.set_title(e if e!='combined' else 'Combined events');ax.set_xlabel('Centroid distance class (km)' if row==0 else 'Grid distance quintile');ax.grid(axis='y',alpha=.2)
axs[0,0].set_ylabel('Fragility — fine classes');axs[1,0].set_ylabel('Fragility — quintiles')
for ext in ['png','pdf']:f.savefig(fig/f'Fig_S3.{ext}',dpi=300)
plt.close(f)
c=pd.read_csv(base/'adjusted_distance_contrasts.csv');f,axs=plt.subplots(1,3,figsize=(12,4.3),sharey=True,layout='constrained');kinds=['grid','block_1deg','block_2deg','block_5deg']
for ax,e,color in zip(axs,events,colors):
 s=c[(c.response==e)&(c.near_km==5)&(c.far_km==20)].set_index('uncertainty').loc[kinds];y=np.arange(4);ax.errorbar(s.fragility_difference,y,xerr=[s.fragility_difference-s.ci_low,s.ci_high-s.fragility_difference],fmt='o',color=color,capsize=4);ax.axvline(0,color='grey',lw=1);ax.set_yticks(y,['Grid','1° blocks','2° blocks','5° blocks']);ax.set_title(e if e!='combined' else 'Combined events');ax.set_xlabel('Adjusted fragility at 20 km minus 5 km');ax.invert_yaxis();ax.grid(axis='x',alpha=.2)
for ext in ['png','pdf']:f.savefig(fig/f'Fig_S5.{ext}',dpi=300)
plt.close(f)

root=Path('.');out=Path('final_results');fig=out/'figures';data=out/'data';full=data/'adjusted';non=data/'noncropland'
u=pd.read_csv(full/'M1_AFB.csv');a=pd.read_csv(full/'M2_AFB.csv');n=pd.read_csv(non/'M2_AFB.csv');land=pd.read_csv(full/'common_landscape_AFB.csv');bal=pd.read_csv(data/'covariate_balance.csv');events=['EVI','LAI','combined'];labels={'EVI':'EVI','LAI':'LAI','combined':'Combined events'}
assert np.allclose(a.M1_AFB,[u[(u.response==r.response)&(u.comparison==r.comparison)].AFB.iloc[0] for r in a.itertuples()])
assert (a.n_model==118498).all() and (a.n_low==59249).all() and (a.n_high==59249).all()
fr=pd.read_csv(full/'M2_fragility.csv')
for r in a.itertuples():
 vals=fr[fr.response.eq(r.response)].set_index('group').fragility
 assert abs(vals[r.comparison]-vals['permanent']-r.AFB)<1e-12
plt.rcParams.update({'font.family':'Times New Roman','font.size':12,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def savefig(f,name):
 for ext in ['png','pdf']:f.savefig(fig/(name+'.'+ext),dpi=300,bbox_inches='tight')
 plt.close(f)
def contrasts(frames,name):
 f,axs=plt.subplots(1,2,figsize=(11,4.8),sharex=True)
 for ax,ctx in zip(axs,['seasonal','other']):
  for j,(frame,title,color,marker) in enumerate(frames):
   for i,e in enumerate(events):
    r=frame[(frame.response==e)&(frame.comparison==ctx)].iloc[0];yy=i+(j-.5)*.2;ax.errorbar(r.AFB,yy,xerr=[[r.AFB-r.ci_low],[r.ci_high-r.AFB]],fmt=marker,color=color,capsize=3,label=title if i==0 else None)
  ax.set_yticks(range(3),[labels[e] for e in events]);ax.invert_yaxis();ax.axvline(0,color='grey',lw=.8);ax.set_title(('Seasonal' if ctx=='seasonal' else 'Other')+' minus permanent');ax.set_xlabel('Fragility contrast (AFB)');ax.grid(axis='x',alpha=.2)
 axs[0].legend(fontsize=10);f.tight_layout();savefig(f,name)
contrasts([(u,'Unadjusted','#777777','o'),(a,'Adjusted','#267895','s')],'Fig_S6_adjusted_contrasts')
f,axs=plt.subplots(1,3,figsize=(11,3.6),sharey=True)
for ax,e in zip(axs,events):
 for j,(frame,title,color) in enumerate([(a,'National adjusted','#267895'),(land,'Common strata','#B67B37')]):
  r=frame[(frame.response==e)&(frame.comparison=='seasonal')].iloc[0];ax.errorbar(r.AFB,j,xerr=[[r.AFB-r.ci_low],[r.ci_high-r.AFB]],fmt='o',color=color,capsize=4)
 ax.set_yticks([0,1],['National adjusted','Common strata']);ax.set_title(labels[e]);ax.axvline(0,color='grey',lw=.8);ax.set_xlabel('Seasonal minus permanent AFB');ax.grid(axis='x',alpha=.2)
f.tight_layout();savefig(f,'Fig_S7_common_landscapes')
f,ax=plt.subplots(figsize=(8,4.4))
for stage,title,color,marker in [('all_corrected_grids','All eligible context grids','#777777','o'),('retained_common_strata_weighted','Retained strata with common weights','#267895','s')]:
 s=bal[bal.stage==stage];ax.plot(abs(s.SMD),range(4),marker=marker,color=color,label=title)
ax.set_yticks(range(4),['Elevation','Slope','Annual temperature','Annual precipitation']);ax.invert_yaxis();ax.axvline(.1,color='grey',ls=':',lw=.8);ax.set_xlabel('Absolute standardized mean difference');ax.legend(fontsize=10);f.tight_layout();savefig(f,'Fig_S8_covariate_balance')
contrasts([(a,'Full sample adjusted','#267895','o'),(n,'Noncropland adjusted','#B67B37','s')],'Fig_S9_noncropland')

F=Path('final_results/figures')
g=pd.read_csv('final_results/data/native_distance_checks.csv');fig,ax=plt.subplots(figsize=(8,3.8));x=np.arange(3)
ax.bar(x-.17,g.native_pixel_geodesic_min_km,width=.34,label='Native geodesic',color='#777777');ax.bar(x+.17,g.UTM_fast_distance_km,width=.34,label='30 m UTM',color='#267895');ax.set(xticks=x,xticklabels=g.grid_id.str.replace('g_','',regex=False),ylabel='Permanent-water distance (km)',xlabel='Centroid longitude_latitude');ax.legend();fig.tight_layout();save(fig,'Fig_S2')
base=Path('final_results/data/robustness')
for file,key,name in [('HDSI_variant_AFB.csv','scenario','Fig_S10'),('event_threshold_AFB.csv','threshold','Fig_S11'),('leave_analysis_year_out_AFB.csv','excluded_year','Fig_S12'),('spatial_block_AFB.csv','block_degrees','Fig_S13')]:
 d=pd.read_csv(base/file);levels=list(d[key].unique());fig,axs=plt.subplots(1,3,figsize=(12,max(3.8,len(levels)*.4+1.4)),sharey=True)
 for ax,ev in zip(axs,['EVI','LAI','combined']):
  for ctx,off,col in [('seasonal',-.13,'#c18a30'),('other',.13,'#267895')]:
   z=d[(d.response==ev)&(d.comparison==ctx)].set_index(key).loc[levels];yy=np.arange(len(levels))+off;ax.errorbar(z.AFB,yy,xerr=[z.AFB-z.ci_low,z.ci_high-z.AFB],fmt='o',color=col,capsize=3,label=ctx.title()+' − permanent')
  ax.axvline(0,color='grey',lw=.8);ax.set(yticks=range(len(levels)),yticklabels=[str(v).replace('omit_','Omit ').replace('only_','Only ') for v in levels],title='Combined events' if ev=='combined' else ev,xlabel='AFB');ax.invert_yaxis()
 axs[0].legend(fontsize=8,loc='best');fig.tight_layout();save(fig,name)

d=pd.read_csv('final_results/data/distance/adjusted_spline_predictions.csv');fig,axs=plt.subplots(1,3,figsize=(12,4),sharey=True)
for ax,ev in zip(axs,['EVI','LAI','combined']):
 z=d[(d.response==ev)&(d.uncertainty=='grid')].sort_values('distance_km');ax.plot(z.distance_km,z.fragility,color='#267895');ax.fill_between(z.distance_km,z.ci_low,z.ci_high,color='#267895',alpha=.2,label='Grid-cluster 95% interval')
 z=d[(d.response==ev)&(d.uncertainty=='block_2deg')].sort_values('distance_km');ax.plot(z.distance_km,z.ci_low,'--',color='#b67b37',label='2°-block 95% interval');ax.plot(z.distance_km,z.ci_high,'--',color='#b67b37');ax.set(title='Combined events' if ev=='combined' else ev,xlabel='Centroid distance (km)');ax.axhline(0,color='grey',lw=.5)
axs[0].set_ylabel('Adjusted fragility');axs[0].legend(fontsize=8);fig.tight_layout();save(fig,'Fig_S4')

# Coverage figure uses the published sample coordinates and composition table.
g=pd.read_csv(data/'distance/sample_coverage.csv')
rep=pd.read_csv(data/'distance/sample_representativeness.csv')
f,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
axs[0].scatter(g.lon,g.lat,s=2,c=np.where(g.available,'#217b9c','#d8d8d8'),rasterized=True)
axs[0].set(xlabel='Longitude',ylabel='Latitude',title='Continuous-distance coverage (blue)');axs[0].set_aspect('equal')
axs[1].barh(rep.variable,rep.SMD,color='#217b9c');axs[1].axvline(0,color='grey')
axs[1].set(xlabel='Standardized mean difference\nAvailable minus unavailable grids',title='Sample composition')
save(f,'Fig_S1')
print('Regenerated 20 quantitative figures; Fig1 is supplied as editable PowerPoint.',flush=True)
