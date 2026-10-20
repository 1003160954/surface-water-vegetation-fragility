from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch,FancyBboxPatch
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
_,d=load_panel();g=d.groupby('grid_id').agg(lon=('lon','first'),lat=('lat','first'),context=('context','first'),combined=('combined','mean'),water_occurrence=('water_occurrence','first'),water_seasonality=('water_seasonality','first'),water_transition=('water_transition','first'),dist_perm_water_km=('dist_perm_water_km','first'),dist_seasonal_water_km=('dist_seasonal_water_km','first'),dist_any_water_km=('dist_any_water_km','first')).reset_index();q=d.HDSI.quantile([.25,.75]);d['exposure']=np.where(d.HDSI<q.iloc[0],'low',np.where(d.HDSI>q.iloc[1],'high','mid'))
colors={'EVI':'#25805c','LAI':'#7857a3','combined':'#b63d49'};contexts=['permanent','seasonal','other'];cc=['#3678a8','#dc9c22','#777777'];maps=[]
def finish(fig,name,source):
 fig.text(.5,-.025,'Audit reconstruction of original-distance inputs; corrected-distance analysis is separate.',ha='center',fontsize=8,color='#7b3a35')
 fig.savefig(ROOT/'figures'/str(name+'.png'),bbox_inches='tight');fig.savefig(ROOT/'figures'/str(name+'.pdf'),bbox_inches='tight');plt.close(fig);maps.append(dict(figure=name,source=source,status='recomputed' if not name.startswith('Fig8') else 'archival redraw only'))
def forest(ax,data,label='response',value='AFB'):
 for i,r in enumerate(data.itertuples()):
  ax.errorbar(getattr(r,value),i,xerr=[[getattr(r,value)-r.ci_low],[r.ci_high-getattr(r,value)]],fmt='o',color=colors.get(r.response,'#333333'),capsize=3)
 ax.set_yticks(range(len(data)),data[label]);ax.axvline(0,color='#777777',lw=.8);ax.set_xlabel('AFB (comparison minus permanent)')
# Figure 1 deterministic flowchart with verified definitions.
fig,ax=plt.subplots(figsize=(8,8));ax.axis('off');labels=['Input panel: 244,306 grid-years / 55,418 grids\nChina subset, 2018-2023; grid spacing 0.1 degrees','Validity filter: valid_model_both_y\n236,995 grid-years / 48,843 grids','Recent seasonal negative anomalies\nEVI < -0.5; LAI < -0.5; combined = either event','Summer HDSI = (VPD + LST - SM - Precip)/4\nLow < Q25; high > Q75; existing seasonal z-scores','Aquatic contexts from numeric centroid distances\nPermanent <=1 km; otherwise seasonal <=1 km; otherwise other','Risk = event proportion; fragility = risk(high) - risk(low)\nAFB = fragility(comparison) - fragility(permanent)','Uncertainty and exploratory diagnostics\nGrid bootstrap; spatial blocks; partial covariate adjustment\nWithin-cell pilot; baseline audit; no causal mechanism inference']
for i,text in enumerate(labels):
 y=.95-i*.135;ax.text(.5,y,text,ha='center',va='center',transform=ax.transAxes,bbox=dict(boxstyle='round,pad=.55',facecolor='#edf3f8',edgecolor='#567a99'))
 if i<len(labels)-1:ax.annotate('',xy=(.5,y-.082),xytext=(.5,y-.045),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color='#567a99'))
finish(fig,'Fig1_verified_workflow','E00/sample_flow.csv; verified numerical definitions')
# Map figure: discrete transition legend, all source categories retained.
fig,axs=plt.subplots(2,3,figsize=(15,8));transcols=['#ffffff','#0000ff','#22b14c','#d1102d','#99d9ea','#b5e61d','#e6a1aa','#ff7f27','#ffc90e','#7f7f7f','#c3c3c3'];transnames=['0 No change','1 Permanent','2 New permanent','3 Lost permanent','4 Seasonal','5 New seasonal','6 Lost seasonal','7 Seasonal to permanent','8 Permanent to seasonal','9 Ephemeral permanent','10 Ephemeral seasonal']
cols=['water_occurrence','water_seasonality','water_transition','dist_perm_water_km','dist_seasonal_water_km','dist_any_water_km']
for ax,col in zip(axs.ravel(),cols):
 if col=='water_transition':
  cm=ListedColormap(transcols);ax.set_facecolor('#eeeeee');sc=ax.scatter(g.lon,g.lat,c=g[col],s=1,cmap=cm,norm=BoundaryNorm(np.arange(-.5,11.5),11));ax.legend(handles=[Patch(facecolor=c,edgecolor='#999999',label=n) for c,n in zip(transcols,transnames)],loc='lower left',fontsize=5,ncol=2)
 else:
  sc=ax.scatter(g.lon,g.lat,c=g[col],s=1,cmap='viridis');fig.colorbar(sc,ax=ax,shrink=.7,label='percent' if col=='water_occurrence' else 'months' if col=='water_seasonality' else 'km')
 ax.set_title(col.replace('_',' '));ax.set_xlabel('Longitude');ax.set_ylabel('Latitude')
fig.suptitle('JRC metrics at sampled analysis-grid centroids; Source: EC JRC/Google');fig.tight_layout();finish(fig,'Fig2_JRC_metrics_discrete_transition','main panel; E00 validity filter; unique grid centroids; JRC GSW v1.4 legend')
save(g,'experiments/E05_figure_audit/grid_map_data.csv')
fig,axs=plt.subplots(1,2,figsize=(11,5));codes=g.context.map(dict(zip(contexts,range(3))));axs[0].scatter(g.lon,g.lat,c=codes,s=1,cmap=ListedColormap(cc),vmin=0,vmax=2);axs[0].legend(handles=[Patch(color=c,label=n) for c,n in zip(cc,contexts)],fontsize=8);sc=axs[1].scatter(g.lon,g.lat,c=g.combined,s=1,cmap='YlOrRd',vmin=0,vmax=1);fig.colorbar(sc,ax=axs[1],label='Combined-event fraction of valid years');axs[0].set_title('Aquatic context');axs[1].set_title('Combined vegetation negative-anomaly events')
for ax in axs:ax.set_xlabel('Longitude');ax.set_ylabel('Latitude')
fig.tight_layout();finish(fig,'Fig3_context_and_combined_events','E05/grid_map_data.csv')
# Grid-cluster bootstrap for group risks (not binomial independent-observation intervals).
def risks(data,group):
 cats=list(data[group].unique());grids=pd.Index(data.grid_id.unique());ix=grids.get_indexer(data.grid_id);N=len(grids);rng=np.random.default_rng(SEED);w=np.array([np.bincount(rng.integers(N,size=N),minlength=N) for _ in range(300)],dtype=np.float32);rows=[]
 for c in cats:
  m=data[group].eq(c).to_numpy();n=np.bincount(ix,weights=m,minlength=N);den=w@n
  for e in colors:
   sums=np.bincount(ix,weights=m*data[e].to_numpy(),minlength=N);b=w@sums/den;lo,hi=np.quantile(b,[.025,.975]);rows.append(dict(group=c,response=e,risk=sums.sum()/n.sum(),ci_low=lo,ci_high=hi,n_total=int(n.sum()),n_grid=int((n>0).sum()),bootstrap=300,seed=SEED))
 return pd.DataFrame(rows)
risk=risks(d,'exposure');save(risk,'experiments/E05_figure_audit/exposure_risks.csv');gh=d.assign(high=d.exposure.eq('high')).groupby('grid_id').high.mean();g['high_HDSI_fraction']=g.grid_id.map(gh);save(g[['grid_id','lon','lat','high_HDSI_fraction']],'experiments/E05_figure_audit/high_HDSI_map.csv')
fig,axs=plt.subplots(1,4,figsize=(17,4));sc=axs[0].scatter(g.lon,g.lat,c=g.high_HDSI_fraction,s=1,cmap='YlOrRd');fig.colorbar(sc,ax=axs[0],label='Fraction of valid years');axs[0].set_title('High-HDSI frequency')
for ax,e in zip(axs[1:],colors):
 t=risk[risk.response==e].set_index('group').loc[['low','mid','high']];ax.bar(range(3),t.risk,color=colors[e]);ax.errorbar(range(3),t.risk,yerr=[t.risk-t.ci_low,t.ci_high-t.risk],fmt='none',ecolor='black',capsize=3);ax.set_xticks(range(3),['low','mid','high']);ax.set_title(e);ax.set_ylabel('Negative-anomaly event risk')
fig.tight_layout();finish(fig,'Fig4_HDSI_and_event_risks','E05/high_HDSI_map.csv; E05/exposure_risks.csv; grid-cluster bootstrap')
band=pd.read_csv(ROOT/'experiments/E00_original_reproduction/fine_distance_fragility.csv');order=sorted(band[band.response=='EVI'].group,key=lambda v: float(v.strip('()[]').split(',')[0]));fig,axs=plt.subplots(1,3,figsize=(16,6))
for ax,e in zip(axs,colors):
 t=band[band.response==e].set_index('group').loc[order];ax.errorbar(range(len(t)),t.fragility,yerr=[t.fragility-t.ci_low,t.ci_high-t.fragility],fmt='o',color=colors[e],capsize=3);labels=[a.replace('(-1e-09','[0')+'\nG='+str(int(n)) for a,n in zip(t.index,t.n_grid)];ax.set_xticks(range(len(t)),labels,rotation=70,ha='right',fontsize=7);ax.set_title(e);ax.set_ylabel('Fragility: risk(high) - risk(low)');ax.set_xlabel('Centroid distance interval (km); G = unique grids');ax.axhline(0,color='grey',lw=.7);ax.text(len(t)-1,.025,'not estimable',rotation=90,ha='center',va='bottom',fontsize=7,transform=ax.get_xaxis_transform())
fig.suptitle('Centroid aquatic contexts; these intervals do not establish within-cell ecological reach');fig.tight_layout();finish(fig,'Fig5_numeric_distance_bands','E00/fine_distance_fragility.csv; numeric permanent distance; 300 grid bootstrap draws')
r=d.copy();r['context_exposure']=r.context+'_'+r.exposure;rr=risks(r,'context_exposure');save(rr,'experiments/E05_figure_audit/context_exposure_risks.csv');frag=pd.read_csv(ROOT/'experiments/E00_original_reproduction/context_fragility.csv');fig,axs=plt.subplots(2,3,figsize=(13,8))
for j,e in enumerate(colors):
 t=rr[(rr.response==e)&rr.group.str.endswith('_high')].copy();t['group']=t.group.str.replace('_high','');t=t.set_index('group').loc[contexts];a=axs[0,j];a.bar(range(3),t.risk,color=cc);a.errorbar(range(3),t.risk,yerr=[t.risk-t.ci_low,t.ci_high-t.risk],fmt='none',ecolor='black',capsize=3);a.set_ylabel('High-HDSI event risk');a.set_title(e)
 t=frag[frag.response==e].set_index('group').loc[contexts];a=axs[1,j];a.bar(range(3),t.fragility,color=cc);a.errorbar(range(3),t.fragility,yerr=[t.fragility-t.ci_low,t.ci_high-t.fragility],fmt='none',ecolor='black',capsize=3);a.set_ylabel('Fragility')
 for a in axs[:,j]:a.set_xticks(range(3),['Permanent','Seasonal','Other']);a.set_xlabel('Aquatic context')
fig.tight_layout();finish(fig,'Fig6_context_risk_and_fragility','E05/context_exposure_risks.csv; E00/context_fragility.csv')
# Figure 7 exploratory summaries; no mechanism claim.
rows=[]
for ctx in contexts:
 s=d[(d.context==ctx)&(d.exposure=='high')]
 for col in ['summer_VPD_z','summer_LST_Day_z','summer_SM_rootproxy_z','summer_Precip_z']:
  rows.append(dict(context=ctx,variable=col,mean=s[col].mean(),n=len(s)))
save(pd.DataFrame(rows),'experiments/E05_figure_audit/exposure_component_means.csv');d['delta_LAI_z']=d.late_summer_LAI_z-d.spring_LAI_z;delta=d.groupby(['context','exposure']).agg(mean_delta=('delta_LAI_z','mean'),n=('delta_LAI_z','count'),n_grid=('grid_id','nunique')).reset_index();save(delta,'experiments/E05_figure_audit/seasonal_LAI_delta.csv')
fig,axs=plt.subplots(1,2,figsize=(11,5));t=pd.DataFrame(rows)
for i,ctx in enumerate(contexts):axs[0].plot(range(4),t[t.context==ctx]['mean'],marker='o',color=cc[i],label=ctx)
axs[0].set_xticks(range(4),['VPD','LST','Soil moisture','Precipitation']);axs[0].set_ylabel('Mean existing seasonal z-score');axs[0].legend();axs[0].set_title('Exposure components under high HDSI')
for i,ctx in enumerate(contexts):
 t=delta[delta.context==ctx].set_index('exposure').reindex(['low','mid','high']);axs[1].plot(range(3),t.mean_delta,marker='o',color=cc[i],label=ctx)
axs[1].set_xticks(range(3),['low','mid','high']);axs[1].set_ylabel('Late-summer LAI z minus spring LAI z');axs[1].set_title('Seasonal LAI change (descriptive means)');fig.suptitle('Exploratory diagnostics; no identified causal pathway');fig.tight_layout();finish(fig,'Fig7_exploratory_diagnostics','E05/exposure_component_means.csv; E05/seasonal_LAI_delta.csv; no inferential intervals')
# Supplementary forest plots.
m0=pd.read_csv(ROOT/'experiments/E02_confounding/M0_AFB.csv');m2=pd.read_csv(ROOT/'experiments/E02_confounding/M2_AFB.csv');fig,axs=plt.subplots(1,2,figsize=(11,4))
for ax,c in zip(axs,['seasonal','other']):
 a=m0[m0.comparison==c].copy();a['label']=a.response+' M0/M1';b=m2[m2.comparison==c].copy();b['label']=b.response+' partial M2';t=pd.concat([a,b],ignore_index=True);forest(ax,t,label='label');ax.set_title(c+' minus permanent')
fig.tight_layout();finish(fig,'SFig1_same_sample_adjustment','E02/M0_AFB.csv; M1 identical; M2_AFB.csv; partial adjustment excludes invalid landcover')
blocks=pd.read_csv(ROOT/'experiments/E04_robustness/spatial_block_AFB.csv');fig,axs=plt.subplots(1,2,figsize=(12,5))
for ax,c in zip(axs,['seasonal','other']):
 t=blocks[blocks.comparison==c].copy();t['label']=t.response+' / '+t.block_degrees.astype(str)+' deg';forest(ax,t,label='label');ax.set_title(c+' minus permanent')
fig.tight_layout();finish(fig,'SFig2_spatial_block_uncertainty','E04/spatial_block_AFB.csv; blocks origin (0,0), 300 draws, seed 20261020')
for filename,scenario,name in [('HDSI_variant_AFB.csv','scenario','SFig3_HDSI_variants'),('event_threshold_AFB.csv','threshold','SFig4_event_thresholds'),('leave_analysis_year_out_AFB.csv','excluded_year','SFig5_leave_analysis_year_out')]:
 table=pd.read_csv(ROOT/'experiments/E04_robustness'/filename);fig,axs=plt.subplots(1,3,figsize=(15,7))
 for ax,e in zip(axs,colors):
  t=table[(table.response==e)&(table.comparison=='seasonal')].copy();t['label']=t[scenario].astype(str);forest(ax,t,label='label');ax.set_title(e+'; seasonal minus permanent')
 fig.tight_layout();finish(fig,name,'E04/'+filename)
sp=pd.read_csv(ROOT/'experiments/E04_robustness/spline_predictions.csv');fig,axs=plt.subplots(1,3,figsize=(14,4))
for ax,e in zip(axs,colors):
 t=sp[sp.response==e];ax.plot(t.distance_km,t.fragility,color=colors[e]);ax.fill_between(t.distance_km,t.ci_low,t.ci_high,color=colors[e],alpha=.2);ax.set_title(e);ax.set_xlabel('Centroid distance (km)');ax.set_ylabel('Adjusted fragility');ax.axhline(0,color='grey',lw=.7)
fig.tight_layout();finish(fig,'SFig6_centroid_distance_spline','E04/spline_predictions.csv; pointwise grid-cluster delta CI; partial adjustment')
save(pd.DataFrame(maps),'experiments/E05_figure_audit/figure_source_map.csv')
print('FIGURES GENERATED',len(maps),flush=True)

