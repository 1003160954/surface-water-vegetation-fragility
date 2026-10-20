from pathlib import Path
import shutil,json
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap,BoundaryNorm
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'revised_figures';F=OUT/'figures';F.mkdir(parents=True,exist_ok=True);D=OUT/'figure_data';D.mkdir(exist_ok=True);T=OUT/'evidence_tables';T.mkdir(exist_ok=True)
P=ROOT/'experiments/E01_spatial_support/corrected_1km';E=ROOT/'experiments/E05_figure_audit'
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'font.family':'DejaVu Sans'})
colors={'EVI':'#25805c','LAI':'#7857a3','combined':'#b63d49'};contextcols={'permanent':'#3678a8','seasonal':'#dc9c22','other':'#777777'};sources=[]
def read(p,name):
 d=pd.read_csv(p);d.to_csv(D/(name+'.csv'),index=False);return d
def finish(fig,name,source):
 fig.savefig(F/(name+'.png'),dpi=200,bbox_inches='tight');fig.savefig(F/(name+'.pdf'),bbox_inches='tight');plt.close(fig);sources.append({'figure':name,'sources':source,'scope':'completed corrected-context revision; no national distance reach inferred'})

fig,ax=plt.subplots(figsize=(8,8));ax.axis('off')
labels=[
 'China panel 2018–2023 at 0.1°\n244,306 rows / 55,418 grids before validity filtering',
 'Valid recent-anomaly sample\n236,995 rows / 48,843 grids; combined event = EVI OR LAI',
 'Corrected water-context measurement\n30 m UTM; centroid distance ≤1 km; JRC 80% / 10 months',
 'Risk and fragility\nSummer HDSI tails; risk(high) minus risk(low)',
 'Context contrasts and comparability\nGrid intervals; common-covariate M2; supported landscape strata',
 'Spatial uncertainty and sensitivity\n1°, 2°, 5° blocks; 10 HDSI indices; 4 thresholds; 6 year omissions',
 'Unfinished analyses remain explicit\nIndependent baseline; full distances; national within-cell exposure']
for i,label in enumerate(labels):
 y=.95-i*.145;ax.text(.5,y,label,ha='center',va='center',transform=ax.transAxes,fontsize=11,bbox={'boxstyle':'round,pad=.6','facecolor':'#edf3f8','edgecolor':'#567a99'})
 if i<6:ax.annotate('',xy=(.5,y-.095),xytext=(.5,y-.05),xycoords='axes fraction',arrowprops={'arrowstyle':'->','color':'#567a99'})
finish(fig,'Fig1_workflow','Executed sample flow, corrected method configuration, analysis status')

g=read(E/'grid_map_data.csv','JRC_grid_attributes')
fig,axs=plt.subplots(1,3,figsize=(14,5.3))
transcols=['#ffffff','#0000ff','#22b14c','#d1102d','#99d9ea','#b5e61d','#e6a1aa','#ff7f27','#ffc90e','#7f7f7f','#c3c3c3']
transnames=['0 No change','1 Permanent','2 New permanent','3 Lost permanent','4 Seasonal','5 New seasonal','6 Lost seasonal','7 Seasonal to permanent','8 Permanent to seasonal','9 Ephemeral permanent','10 Ephemeral seasonal']
for ax,col,label in zip(axs,['water_occurrence','water_seasonality','water_transition'],['a  Occurrence','b  Seasonality','c  Transition']):
 if col=='water_transition':
  ax.set_facecolor('#eeeeee');ax.scatter(g.lon,g.lat,c=g[col],s=1,cmap=ListedColormap(transcols),norm=BoundaryNorm(np.arange(-.5,11.5),11))
 else:
  sc=ax.scatter(g.lon,g.lat,c=g[col],s=1,cmap='viridis');fig.colorbar(sc,ax=ax,shrink=.72,label='%' if col=='water_occurrence' else 'months')
 ax.set_title(label);ax.set_xlabel('Longitude');ax.set_ylabel('Latitude')
fig.legend(handles=[Patch(facecolor=c,edgecolor='#777777',label=n) for c,n in zip(transcols,transnames)],loc='lower center',ncol=4,fontsize=12,bbox_to_anchor=(.5,-.015));fig.tight_layout(rect=(0,.18,1,1));finish(fig,'Fig2_JRC_attributes','figure_data/JRC_grid_attributes.csv; JRC GSW v1.4 discrete transition labels')

ctx=read(P/'corrected_context_map.csv','corrected_context_map');g=g.merge(ctx[['grid_id','context']].rename(columns={'context':'corrected_context'}),on='grid_id',validate='one_to_one')
fig,axs=plt.subplots(1,2,figsize=(11,5.4))
for c in ['other','seasonal','permanent']:
 t=g[g.corrected_context==c];axs[0].scatter(t.lon,t.lat,s=2,color=contextcols[c],label=c+' (G='+str(len(t))+')')
axs[0].legend(fontsize=9,loc='lower left');axs[0].set_title('a  Corrected 1 km centroid contexts')
sc=axs[1].scatter(g.lon,g.lat,c=g.combined,s=2,cmap='YlOrRd',vmin=0,vmax=1);fig.colorbar(sc,ax=axs[1],label='Fraction of valid years');axs[1].set_title('b  Combined negative-anomaly frequency')
for ax in axs:ax.set_xlabel('Longitude');ax.set_ylabel('Latitude')
fig.tight_layout();finish(fig,'Fig3_corrected_contexts','figure_data/corrected_context_map.csv and JRC_grid_attributes.csv; counts are unique grids')

r=read(E/'exposure_risks.csv','exposure_risks');m=read(E/'high_HDSI_map.csv','high_HDSI_map')
fig,axs=plt.subplots(2,2,figsize=(10,8));sc=axs[0,0].scatter(m.lon,m.lat,c=m.high_HDSI_fraction,s=2,cmap='YlOrRd');fig.colorbar(sc,ax=axs[0,0],label='Fraction of valid years');axs[0,0].set_title('a  High-HDSI exposure');axs[0,0].set_xlabel('Longitude');axs[0,0].set_ylabel('Latitude')
for ax,(e,color),label in zip(list(axs.flat)[1:],colors.items(),['b','c','d']):
 t=r[r.response==e].set_index('group').loc[['low','mid','high']];ax.bar(range(3),t.risk,color=color);ax.errorbar(range(3),t.risk,yerr=[t.risk-t.ci_low,t.ci_high-t.risk],fmt='none',ecolor='black',capsize=3);ax.set_xticks(range(3),['Low','Middle','High']);ax.set_title(label+'  '+e);ax.set_xlabel('HDSI group');ax.set_ylabel('Negative-anomaly risk')
fig.tight_layout();finish(fig,'Fig4_exposure_risks','figure_data/exposure_risks.csv and high_HDSI_map.csv; 300 paired grid bootstrap draws')

geo=read(ROOT/'experiments/E01_spatial_support/independent_native_distance_check_original_threshold.csv','distance_native_UTM_audit');cross=read(P/'context_reclassification_grid.csv','context_reclassification_grid')
fig,axs=plt.subplots(1,2,figsize=(11,5));xx=np.arange(len(geo));w=.25
for i,(col,label,color) in enumerate([('original_table_km','Original table','#888888'),('native_pixel_geodesic_min_km','Native geodesic','#d18f22'),('UTM_fast_distance_km','30 m UTM','#3678a8')]):axs[0].bar(xx+(i-1)*w,geo[col],width=w,label=label,color=color)
axs[0].set_xticks(xx,[v.replace('g','') for v in geo.grid_id],rotation=20,ha='right',fontsize=9);axs[0].set_ylabel('Distance to permanent water (km)');axs[0].set_title('a  Independent physical checks');axs[0].legend(fontsize=9)
mat=cross.set_index('original_context').reindex(index=['permanent','seasonal','other'],columns=['permanent','seasonal','other']);axs[1].imshow(mat.to_numpy(),cmap='Blues')
for i in range(3):
 for j in range(3):axs[1].text(j,i,f'{int(mat.iloc[i,j]):,}',ha='center',va='center',color='white' if mat.iloc[i,j]>15000 else 'black')
axs[1].set_xticks(range(3),['Permanent','Seasonal','Other']);axs[1].set_yticks(range(3),['Permanent','Seasonal','Other']);axs[1].set_xlabel('Corrected context');axs[1].set_ylabel('Original context');axs[1].set_title('b  Reclassification of 48,843 grids');fig.tight_layout();finish(fig,'Fig5_distance_audit','figure_data/distance_native_UTM_audit.csv and context_reclassification_grid.csv')

m0=read(P/'AFB.csv','corrected_M0_AFB');m2=read(P/'harmonized/M2_AFB.csv','corrected_M2_AFB');land=read(P/'harmonized/common_landscape_AFB.csv','corrected_landscape_AFB');blocks=read(P/'robustness/spatial_block_AFB.csv','corrected_spatial_block_AFB')
fig,axs=plt.subplots(1,2,figsize=(11,5.5))
for ax,c in zip(axs,['seasonal','other']):
 for j,e in enumerate(colors):
  for k,(d,marker,label) in enumerate([(m0,'o','M0'),(m2,'s','M2')]):
   rr=d[(d.response==e)&(d.comparison==c)].iloc[0];y=j*2+k;ax.errorbar(rr.AFB,y,xerr=[[rr.AFB-rr.ci_low],[rr.ci_high-rr.AFB]],fmt=marker,color=colors[e],capsize=3)
 ax.set_yticks(range(6),[e+' / '+label for e in colors for label in ['M0','M2']]);ax.set_ylim(5.6,-.6);ax.axvline(0,color='grey',lw=.8);ax.set_xlabel('AFB (absolute risk difference)');ax.set_title(c.capitalize()+' minus permanent')
fig.tight_layout();finish(fig,'Fig6_corrected_contrasts','figure_data/corrected_M0_AFB.csv; M2_AFB.csv; M0 grid percentile CI versus M2 grid-cluster delta CI')

fig,ax=plt.subplots(figsize=(8,4.2))
for i,(e,color) in enumerate(colors.items()):
 rr=land[land.response==e].iloc[0];ax.errorbar(rr.AFB,i,xerr=[[rr.AFB-rr.ci_low],[rr.ci_high-rr.AFB]],fmt='o',color=color,capsize=4)
ax.set_yticks(range(3),list(colors));ax.set_ylim(2.5,-.5);ax.axvline(0,color='grey',lw=.8);ax.set_xlabel('Seasonal minus permanent AFB');ax.set_title('46 supported strata / 13,881 grids / 39,524 rows');fig.tight_layout();finish(fig,'Fig7_landscape_contrasts','figure_data/corrected_landscape_AFB.csv; 296 valid of 300 paired grid bootstrap draws; fixed weights')

fig,axs=plt.subplots(1,2,figsize=(11,6))
for ax,c in zip(axs,['seasonal','other']):
 t=blocks[blocks.comparison==c].sort_values(['block_degrees','response'])
 for i,rr in enumerate(t.itertuples()):ax.errorbar(rr.AFB,i,xerr=[[rr.AFB-rr.ci_low],[rr.ci_high-rr.AFB]],fmt='o',color=colors[rr.response],capsize=3)
 ax.set_yticks(range(len(t)),[rr.response+' / '+str(rr.block_degrees)+'°' for rr in t.itertuples()]);ax.set_ylim(len(t)-.5,-.5);ax.axvline(0,color='grey',lw=.8);ax.set_xlabel('AFB (absolute risk difference)');ax.set_title(c.capitalize()+' minus permanent')
fig.tight_layout();finish(fig,'Fig8_spatial_uncertainty','figure_data/corrected_spatial_block_AFB.csv; 300 paired block resamples; origin0,0')
pd.DataFrame(sources).to_csv(OUT/'figure_source_map.csv',index=False)

a=m0.copy();a['method']='M0 paired-grid bootstrap';b=m2.copy();b['method']='M2 grid-cluster sandwich delta'
pd.concat([a,b],ignore_index=True).to_csv(T/'Table_R1_corrected_context_contrasts.csv',index=False)
shutil.copy2(P/'context_fragility.csv',T/'Table_R1_context_risks_and_counts.csv')
land.to_csv(T/'Table_R2_landscape_contrasts.csv',index=False)
shutil.copy2(P/'harmonized/retained_strata_support.csv',T/'Table_R2_stratum_support.csv')
blocks.to_csv(T/'Table_R3_spatial_blocks.csv',index=False)
shutil.copy2(P/'robustness/scenario_summary.csv',T/'Table_R4_scenario_summary.csv')
for n in ['HDSI_variant_AFB','event_threshold_AFB','leave_analysis_year_out_AFB']:
 shutil.copy2(P/'robustness'/(n+'.csv'),T/('Table_R4_'+n+'.csv'))
geo.to_csv(T/'Table_R5_distance_checks.csv',index=False);cross.to_csv(T/'Table_R5_reclassification_grid_counts.csv',index=False)
shutil.copy2(ROOT/'experiments/E00_original_reproduction/original_regional_AFB_rerun_with_omissions.csv',T/'Table_R6_original_regional_audit_ONLY.csv')
pd.DataFrame([
 {'item':'National corrected 1km centroid classification','status':'complete','scope':'48843 grids; 17056-grid available full-distance check: zero binary mismatches'},
 {'item':'Corrected AFB, adjustment, landscape strata, stress/threshold/year sensitivities, blocks','status':'complete','scope':'recorded model and bootstrap outputs; residual confounding remains'},
 {'item':'Full national distance magnitudes and within-cell metrics','status':'incomplete; paused by user','scope':'completed exports preserved; no resubmission'},
 {'item':'Reaggregation of vegetation and climate after grid shift/resolution change','status':'not executed','scope':'metric-only 24-grid pilot is not this test'},
 {'item':'Independent 2001-2017 anomaly reference','status':'not executed','scope':'processor software checks only; real historical inputs unavailable'},
 {'item':'Original vegetation/thermal collection versions and quality masks','status':'not fully verified','scope':'do not infer IDs from generic MODIS product names'},
 {'item':'Corrected national mechanism and regional effects','status':'not established','scope':'old-context outputs are archived, not transferred'},
 {'item':'Irrigation, groundwater, soil and dynamic management adjustment','status':'not included','scope':'background climate and modal land cover are included'},
 {'item':'Public GitHub/data release','status':'not created','scope':'local Git and evidence archive only'}
]).to_csv(T/'Table_R7_completion_and_provenance_gaps.csv',index=False)
print('EIGHT REVISION FIGURES AND EVIDENCE TABLES SAVED')
