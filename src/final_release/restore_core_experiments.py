from pathlib import Path
import sys,ast,json
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
b=Path.cwd();project=b;sys.path.insert(0,str(project/'src'))
from audit_pipeline import load_panel,group_effects
out=b/'outputs/manuscript_revision_v4';out.mkdir(parents=True,exist_ok=True);t=out/'supporting_data';t.mkdir(exist_ok=True);fig=out/'figures';fig.mkdir(exist_ok=True)
_,d=load_panel();c=pd.read_csv(project/'experiments/E01_spatial_support/national_corrected_1km_contexts.csv');d=d.merge(c,on='grid_id',validate='many_to_one');d['context']=np.where(d.near_perm_1km_corrected==1,'permanent',np.where(d.near_seasonal_1km_corrected==1,'seasonal','other'));q=d.HDSI.quantile([.25,.75]);d['exposure']=np.where(d.HDSI<q.iloc[0],'low',np.where(d.HDSI>q.iloc[1],'high','middle'));d['LAI_seasonal_change']=d.late_summer_LAI_z-d.spring_LAI_z
# Bootstrap whole grids, retaining all available years and fixed national HDSI cutpoints.
grids=pd.Index(d.grid_id.unique());ix=grids.get_indexer(d.grid_id);rng=np.random.default_rng(20261020);W=np.array([np.bincount(rng.integers(len(grids),size=len(grids)),minlength=len(grids)) for _ in range(300)],dtype=np.float32)
rows=[];boots={}
for metric,col,sign in [('VPD','summer_VPD_z',1),('surface_heat','summer_LST_Day_z',1),('soil_moisture_deficit','summer_SM_rootproxy_z',-1),('precipitation_deficit','summer_Precip_z',-1),('air_temperature','summer_Tair_z',1),('LAI_seasonal_change','LAI_seasonal_change',1)]:
 for ctx in ['permanent','seasonal','other']:
  for lev in (['low','high'] if metric=='LAI_seasonal_change' else ['high']):
   m=(d.context.eq(ctx)&d.exposure.eq(lev)&d[col].notna()).to_numpy();vals=d[col].fillna(0).to_numpy()*sign;count=np.bincount(ix,weights=m,minlength=len(grids));sums=np.bincount(ix,weights=m*vals,minlength=len(grids));den=W@count;v=(W@sums)/np.where(den>0,den,np.nan);boots[metric,ctx,lev]=v
   rows.append(dict(metric=metric,context=ctx,exposure=lev,n_rows=int(m.sum()),n_grid=int((count>0).sum()),mean=sums.sum()/count.sum(),ci_low=np.nanquantile(v,.025),ci_high=np.nanquantile(v,.975),valid_draws=int(np.isfinite(v).sum())))
means=pd.DataFrame(rows);means.to_csv(t/'corrected_component_and_seasonal_means.csv',index=False)
contr=[]
for (metric,ctx,lev),v in boots.items():
 if ctx=='permanent':continue
 base=means.query('metric==@metric and context=="permanent" and exposure==@lev').iloc[0];a=means.query('metric==@metric and context==@ctx and exposure==@lev').iloc[0];vv=v-boots[metric,'permanent',lev]
 contr.append(dict(metric=metric,comparison=ctx,exposure=lev,difference=a['mean']-base['mean'],ci_low=np.nanquantile(vv,.025),ci_high=np.nanquantile(vv,.975)))
pd.DataFrame(contr).to_csv(t/'corrected_diagnostic_contrasts.csv',index=False)
del W
src=(project/'docs/original_sources/cell_59.py').read_text(encoding='utf8');node=next(n for n in ast.parse(src).body if isinstance(n,ast.FunctionDef) and n.name=='assign_macro_region');exec(compile(ast.Module(body=[node],type_ignores=[]),'<original-region-rule>','exec'));d['macro_region']=d.apply(assign_macro_region,axis=1)
d.groupby(['macro_region','context','exposure']).agg(n_rows=('grid_id','size'),n_grid=('grid_id','nunique')).reset_index().to_csv(t/'corrected_regional_support.csv',index=False)
regionrows=[]
for reg,s in d.groupby('macro_region'):
 r,a=group_effects(s,B=300,q_fixed=q)
 for ctx in ['seasonal','other']:
  counts=[int((s.context.eq(cc)&s.exposure.eq(lev)).sum()) for cc in ['permanent',ctx] for lev in ['low','high']];totals=[int(s.context.eq(cc).sum()) for cc in ['permanent',ctx]];supported=min(counts)>=30 and min(totals)>=100
  for response in ['EVI','LAI','combined']:
   z=a.loc[a.comparison.eq(ctx)&a.response.eq(response)] if len(a) else a
   row=dict(region=reg,comparison=ctx,response=response,supported=supported,n_permanent_low=counts[0],n_permanent_high=counts[1],n_comparison_low=counts[2],n_comparison_high=counts[3])
   if len(z):row.update(z.iloc[0][['AFB','ci_low','ci_high','valid_draws']].to_dict())
   regionrows.append(row)
regions=pd.DataFrame(regionrows);regions.to_csv(t/'corrected_regional_AFB.csv',index=False)
colors={'permanent':'#2878a5','seasonal':'#dc9f24','other':'#777777'}
f,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for j,(ctx,color) in enumerate(colors.items()):
 sub=means[(means.metric!='LAI_seasonal_change')&means.context.eq(ctx)];x=np.arange(5)+(j-1)*.22;axes[0,0].errorbar(x,sub['mean'],yerr=[sub['mean']-sub.ci_low,sub.ci_high-sub['mean']],fmt='o',color=color,label=ctx,capsize=3)
 sub=means[means.metric.eq('LAI_seasonal_change')&means.context.eq(ctx)];axes[0,1].errorbar(np.arange(2)+(j-1)*.15,sub['mean'],yerr=[sub['mean']-sub.ci_low,sub.ci_high-sub['mean']],fmt='o-',color=color,label=ctx,capsize=3)
axes[0,0].set(xticks=range(5),xticklabels=['VPD','LST','SM deficit','P deficit','Tair'],ylabel='Mean standardized component',title='a  Components under high HDSI');axes[0,0].legend(fontsize=8)
axes[0,1].set(xticks=[0,1],xticklabels=['Low HDSI','High HDSI'],ylabel='Late-summer minus spring LAI z',title='b  Seasonal LAI change');axes[0,1].axhline(0,color='grey',lw=.7)
risks=pd.read_csv(project/'experiments/E01_spatial_support/corrected_1km/context_fragility.csv')
for ax,event,label in zip(axes[1],['LAI','combined'],['c  LAI event risks','d  Combined-event risks']):
 for j,(ctx,color) in enumerate(colors.items()):
  z=risks.query('group==@ctx and response==@event').iloc[0];ax.plot([z.risk_low,z.risk_high],[j,j],'o-',color=color,label=ctx)
 ax.set(yticks=range(3),yticklabels=list(colors),xlabel='Negative-anomaly event probability',title=label)
for ext in ['png','pdf']:f.savefig(fig/f'Fig7_exposure_and_seasonal_diagnostics.{ext}',dpi=200)
plt.close(f)
f,axes=plt.subplots(1,3,figsize=(13,5),layout='constrained');regs=sorted(regions.region.unique())
for ax,event in zip(axes,['EVI','LAI','combined']):
 for j,reg in enumerate(regs):
  z=regions.query('region==@reg and response==@event and comparison=="other"').iloc[0]
  if np.isfinite(z.get('AFB',np.nan)):
   ax.errorbar(z.AFB,j,xerr=[[z.AFB-z.ci_low],[z.ci_high-z.AFB]],fmt='o',mfc='#2878a5' if z.supported else 'white',mec='#2878a5',ecolor='#2878a5',capsize=3)
  else:ax.text(0,j,'NE',fontsize=8)
 ax.axvline(0,color='grey',lw=.8);ax.set(yticks=range(len(regs)),yticklabels=[x+(' *' if not regions.query('region==@x and response==@event and comparison=="other"').iloc[0].supported else '') for x in regs],title=event,xlabel='Other minus permanent AFB')
f.suptitle('Regional contrasts; * limited sample support; NE not estimable')
for ext in ['png','pdf']:f.savefig(fig/f'Fig8_corrected_regional_comparison.{ext}',dpi=200)
plt.close(f)
(t/'diagnostic_config.json').write_text(json.dumps({'national_q25':q.iloc[0],'national_q75':q.iloc[1],'bootstrap':300,'resampling':'paired grid cluster, all years retained','seed':20261020,'region_rule':'original assign_macro_region; coordinate/elevation rules, not administrative polygons','support_rule':'both context totals >=100 grid-years; all four low/high cells >=30; no effect-sign selection','causal_mediation':False,'heavy_GEE_tasks':'remain paused'},indent=2))
print(means.query('metric=="LAI_seasonal_change"').to_string(index=False));print(regions.query('comparison=="other"').to_string(index=False));print('Core diagnostics recomputed',flush=True)
