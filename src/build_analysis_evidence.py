from pathlib import Path
import sys,json,hashlib,html,shutil
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from audit_pipeline import load_panel
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=ROOT/'analysis_evidence';O.mkdir(exist_ok=True);F=O/'figures';F.mkdir(exist_ok=True)
P=ROOT/'experiments/E01_spatial_support/corrected_1km';N=ROOT/'experiments/E02_confounding/corrected_noncrop'
raw,d=load_panel();flags=pd.read_csv(ROOT/'experiments/E01_spatial_support/national_corrected_1km_contexts.csv');d=d.merge(flags,on='grid_id',validate='many_to_one');d['original_context']=d.context;d['context']=np.where(d.near_perm_1km_corrected==1,'permanent',np.where(d.near_seasonal_1km_corrected==1,'seasonal','other'))
q=d.HDSI.quantile([.25,.75]);d['exposure']=np.where(d.HDSI<q.iloc[0],'low',np.where(d.HDSI>q.iloc[1],'high','middle'))
rows=[]
for ctx in ['permanent','seasonal','other']:
 for response in ['EVI','LAI','combined']:
  s=d[d.context==ctx];lo=s[s.exposure=='low'];hi=s[s.exposure=='high'];a=int(lo[response].sum());b=int(hi[response].sum())
  rows.append(dict(context=ctx,response=response,n_context_rows=len(s),n_context_grids=s.grid_id.nunique(),low_events=a,low_n=len(lo),high_events=b,high_n=len(hi),risk_low=a/len(lo),risk_high=b/len(hi),fragility=b/len(hi)-a/len(lo)))
counts=pd.DataFrame(rows);counts.to_csv(O/'risk_numerators_and_arithmetic.csv',index=False)
expected=pd.read_csv(P/'context_fragility.csv');test=counts.merge(expected,on=['response'],suffixes=('_fresh','_stored'));test=test[test.context==test.group]
assert abs(test.fragility_fresh-test.fragility_stored).max()<1e-13
afb=[]
for response in ['EVI','LAI','combined']:
 b=counts[(counts.context=='permanent')&(counts.response==response)].fragility.iloc[0]
 for ctx in ['seasonal','other']:
  a=counts[(counts.context==ctx)&(counts.response==response)].fragility.iloc[0];afb.append(dict(response=response,comparison=ctx,comparison_fragility=a,permanent_fragility=b,AFB=a-b))
pd.DataFrame(afb).to_csv(O/'AFB_arithmetic.csv',index=False)
ids=[]
for i,ctx in enumerate(['permanent','seasonal','other']):ids+=d[d.context==ctx].grid_id.drop_duplicates().sample(n=2,random_state=20261020+i).tolist()
cols=['grid_id','year','lon','lat','summer_VPD_z','summer_LST_Day_z','summer_SM_rootproxy_z','summer_Precip_z','HDSI','late_summer_EVI_z','late_summer_LAI_z','EVI','LAI','combined','original_context','context','near_perm_1km_corrected','near_seasonal_1km_corrected','exposure']
examples=d[d.grid_id.isin(ids)][cols].copy();examples['HDSI_recomputed']=(examples.summer_VPD_z+examples.summer_LST_Day_z-examples.summer_SM_rootproxy_z-examples.summer_Precip_z)/4;examples['HDSI_Q25']=q.iloc[0];examples['HDSI_Q75']=q.iloc[1];examples.sort_values(['grid_id','year']).to_csv(O/'worked_grid_year_examples.csv',index=False)
dictionary=pd.read_csv(ROOT/'experiments/E00_original_reproduction/data_dictionary.csv')
def source_for(f):
 if f.startswith(('spring_','summer_','late_summer_','autumn_')):return 'Supplied seasonal panel; recovered seasonal recipe; exact original collection IDs/QA partly unresolved'
 if 'water' in f:return 'Original JRC v1.4 outputs; original distances retained only for audit; use corrected1km flags for inference'
 if f=='HDSI':return '(summer_VPD_z + summer_LST_Day_z - summer_SM_rootproxy_z - summer_Precip_z) / 4'
 return 'Supplied analysis export; field-level dtype/missingness recorded'
dictionary['lineage_note']=dictionary.field.map(source_for);dictionary.to_csv(O/'all_panel_fields_dictionary.csv',index=False)

# Balance is a descriptive check on observed covariates, not proof of exchangeability.
cv=pd.read_csv(ROOT/'experiments/E02_confounding/harmonized_landcover_climate_covariates.csv');d=d.merge(cv,on='grid_id',validate='many_to_one');grid=d.drop_duplicates('grid_id').copy();s=d[d.context.isin(['permanent','seasonal'])&(d.exposure!='middle')].copy();config=json.loads((P/'harmonized/stratification_config.json').read_text());s['stratum']=s.LC2018_grid_mode.astype(int).astype(str)
for col in ['climate_annual_temp_C','climate_annual_precip_mm','elevation']:s['stratum']+='_'+pd.cut(s[col],config['cuts'][col],labels=False).astype(str)
keep=pd.read_csv(P/'harmonized/retained_strata_support.csv').stratum.unique();s=s[s.stratum.isin(keep)].copy();assert len(s)==39524
g=s.drop_duplicates('grid_id').copy();ws=g.groupby('stratum').size();ws=ws/ws.sum();n_by=g.groupby(['stratum','context']).size();g['weight']=[ws[r.stratum]/n_by[r.stratum,r.context] for r in g.itertuples()]
bal=[]
fields=['elevation','slope','climate_annual_temp_C','climate_annual_precip_mm']
for stage,gg in [('all_corrected_grids',grid[grid.context.isin(['permanent','seasonal'])].assign(weight=1.)),('retained_common_strata_weighted',g)]:
 for col in fields:
  stats=[]
  for ctx in ['permanent','seasonal']:
   x=gg[gg.context==ctx];w=x.weight.to_numpy();v=x[col].to_numpy();w=w/w.sum();mu=float(w@v);var=float(w@((v-mu)**2));stats.append((mu,var))
  sd=np.sqrt((stats[0][1]+stats[1][1])/2);bal.append(dict(stage=stage,covariate=col,permanent_mean=stats[0][0],seasonal_mean=stats[1][0],SMD=(stats[1][0]-stats[0][0])/sd,n_grids=gg.grid_id.nunique()))
balance=pd.DataFrame(bal);balance.to_csv(O/'corrected_covariate_balance.csv',index=False)
(O/'covariate_balance_scope.json').write_text(json.dumps({'comparison':'seasonal minus permanent','before':'all corrected permanent/seasonal grids','after':'unique grids with low/high observations in46 retained common strata, pooled-grid stratum weights shared across contexts','SMD':'difference in weighted means / sqrt(mean of two weighted descriptive variances)','warning':'before/after target populations differ; coarse strata do not guarantee all covariates balanced; no irrigation/groundwater/soil balance can be claimed'},indent=2))
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False});colors={'EVI':'#25805c','LAI':'#7857a3','combined':'#b63d49'};maps=[]
def finish(fig,name,source):
 fig.savefig(F/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(F/(name+'.pdf'),bbox_inches='tight');plt.close(fig);maps.append(dict(figure=name,source=source,scope='corrected1km contexts; executed local analyses'))
full=pd.read_csv(P/'harmonized/M2_AFB.csv');non=pd.read_csv(N/'M2_AFB.csv');fig,axs=plt.subplots(1,2,figsize=(12,5.4))
for ax,ctx in zip(axs,['seasonal','other']):
 for i,(response,color) in enumerate(colors.items()):
  for j,(a,label,marker) in enumerate([(full,'Full sample','o'),(non,'Exclude crop classes12/14','s')]):
   r=a[(a.response==response)&(a.comparison==ctx)].iloc[0];y=i*2+j;ax.errorbar(r.AFB,y,xerr=[[r.AFB-r.ci_low],[r.ci_high-r.AFB]],fmt=marker,color=color,capsize=3)
 ax.set_yticks(range(6),[r+' / '+l for r in colors for l in ['Full','Noncrop']]);ax.set_ylim(5.5,-.5);ax.axvline(0,color='grey',lw=.8);ax.set_title(ctx+' minus permanent');ax.set_xlabel('Adjusted AFB; grid-cluster delta95%CI')
fig.suptitle('Agricultural-composition sensitivity; not an irrigation-free experiment');fig.tight_layout();finish(fig,'SFig12_noncrop_sensitivity','E02/corrected_noncrop/M2_AFB.csv; full corrected1km harmonized/M2_AFB.csv')
fig,ax=plt.subplots(figsize=(9,4.5))
for stage,marker,label in [('all_corrected_grids','o','All corrected grids'),('retained_common_strata_weighted','s','Selected common strata / common weights')]:
 t=balance[balance.stage==stage];ax.plot(abs(t.SMD),range(4),marker=marker,label=label)
ax.set_yticks(range(4),['Elevation','Slope','Annual temperature','Annual precipitation']);ax.set_xlabel('Absolute descriptive standardized mean difference');ax.axvline(.1,color='grey',ls=':',lw=.8);ax.legend(fontsize=9);ax.set_title('Observed-covariate balance; target populations differ');fig.tight_layout();finish(fig,'SFig13_corrected_balance','analysis_evidence/corrected_covariate_balance.csv; .1 line is a descriptive reference, not a validity threshold')
for fn,scenario,name in [('HDSI_variant_AFB.csv','scenario','SFig14_corrected_HDSI'),('event_threshold_AFB.csv','threshold','SFig15_corrected_thresholds'),('leave_analysis_year_out_AFB.csv','excluded_year','SFig16_corrected_year_omissions')]:
 a=pd.read_csv(P/'robustness'/fn);fig,axs=plt.subplots(1,3,figsize=(14,6))
 for ax,(response,color) in zip(axs,colors.items()):
  t=a[(a.response==response)&(a.comparison=='seasonal')]
  for i,r in enumerate(t.itertuples()):ax.errorbar(r.AFB,i,xerr=[[r.AFB-r.ci_low],[r.ci_high-r.AFB]],fmt='o',color=color,capsize=3)
  ax.set_yticks(range(len(t)),t[scenario].astype(str));ax.set_ylim(len(t)-.5,-.5);ax.axvline(0,color='grey',lw=.8);ax.set_title(response);ax.set_xlabel('Seasonal minus permanent AFB')
 fig.suptitle('Corrected contexts; grid-bootstrap95%intervals');fig.tight_layout();finish(fig,name,'corrected1km/robustness/'+fn)
flow=pd.read_csv(ROOT/'experiments/E00_original_reproduction/sample_flow.csv');fig,ax=plt.subplots(figsize=(8,4));ax.barh(range(4),flow.n_rows,color='#3678a8');ax.set_yticks(range(4),flow.stage);ax.invert_yaxis();ax.set_xlabel('Grid-year records')
for i,r in enumerate(flow.itertuples()):ax.text(r.n_rows*.97,i,f'{r.n_rows:,} rows / {r.n_grids:,} grids',ha='right',va='center',color='white')
fig.tight_layout();finish(fig,'SFig17_sample_flow','E00/sample_flow.csv; no new filtering invented')
pd.DataFrame(maps).to_csv(O/'figure_source_map.csv',index=False)
protection=[
 dict(factor='Different mean vegetation levels',implemented='grid x seasonal-position recent anomaly standardization',addresses='location/season mean and scale differences',does_not_address='different underlying sensitivities or unstable short references'),
 dict(factor='Land-cover composition',implemented='categorical2018 grid modal land cover in M2',addresses='observed broad cover differences',does_not_address='within-class management or dynamic crop changes'),
 dict(factor='Terrain and background climate',implemented='elevation,slope,WorldClim MAT/MAP,quadratic geography,year indicators',addresses='measured landscape and broad location differences',does_not_address='unobserved irrigation,groundwater,soil; guarantee of causal exchangeability'),
 dict(factor='Comparable landscapes',implemented='46 supported LC x MAT/MAP/elevation-tertile strata; common weights',addresses='restricts to shared coarse settings',does_not_address='all within-stratum differences; different target population'),
 dict(factor='Agricultural landscape contribution',implemented='exclude modal LCclasses12and14; fixed national HDSI cutpoints; refit M0/M1/M2',addresses='dependence on dominant cropland and mosaic classes',does_not_address='identifying irrigation-free or groundwater-independent land'),
 dict(factor='Repeated years and neighbouring cells',implemented='paired grid resampling,cluster sandwich,1/2/5degree blocks',addresses='specified temporal-cluster and spatial resampling dependence',does_not_address='unmeasured confounding or established true spatial correlation range')]
pd.DataFrame(protection).to_csv(O/'confounding_controls_and_limits.csv',index=False)
stage_rows=[
 ('01','Original data products','docs/data_sources.csv; analysis_evidence/all_panel_fields_dictionary.csv','Known JRC/land-cover/climate/pilot products; original vegetation/thermal IDs partly unresolved','gee/; docs/original_sources/upstream_seasonal_source_only.ipynb'),
 ('02','Supplied raw analysis export','logs/E00_manifest.json; analysis_evidence/raw_input_manifest.json','244306 grid-years,55418grids; no duplicate grid-year keys','src/audit_pipeline.py'),
 ('03','Validity and seasonal anomaly recipe','experiments/E00_original_reproduction/sample_flow.csv','236995 rows,48843grids; at most6 annual observations per seasonal position','docs/original_sources/upstream_seasonal_source_only.ipynb'),
 ('04','Stress and event assignments','analysis_evidence/worked_grid_year_examples.csv','HDSI four-component mean; Q25/Q75; EVI/LAI<-.5; OR event','src/audit_pipeline.py'),
 ('05','Corrected aquatic contexts','experiments/E01_spatial_support/national_corrected_1km_contexts.csv','2549permanent,15262seasonal,31032other; 17056-grid binary check zero mismatches','gee/corrected_near_contexts.py'),
 ('06','Event risk and fragility','analysis_evidence/risk_numerators_and_arithmetic.csv; analysis_evidence/AFB_arithmetic.csv','Integer event counts / denominators -> high-low risk -> signed context contrast','src/audit_pipeline.py; src/analyze_corrected_1km.py'),
 ('07','Measured-confounder adjustment','experiments/E01_spatial_support/corrected_1km/harmonized/','M0/M1 same eligible sample; M2 common-covariate predictions,grid-cluster delta95%CI','src/harmonized_confounding.py; src/analyze_corrected_1km.py'),
 ('08','Landscape and agricultural sensitivity','analysis_evidence/corrected_covariate_balance.csv; experiments/E02_confounding/corrected_noncrop/','46strata; excludes crop/mosaic in additional analysis; no full irrigation control claimed','src/common_landscape.py; src/noncropland_sensitivity.py'),
 ('09','Robustness and spatial uncertainty','experiments/E01_spatial_support/corrected_1km/robustness/','10 indices,4 thresholds,6 analysis-year omissions,1/2/5degree blocks; all scenarios retained','src/robustness.py; src/analyze_corrected_1km.py'),
 ('10','Quantitative visualizations','diagnostic_figures/figure_source_map.csv; analysis_evidence/figure_source_map.csv','Main8 figures plus corrected support plots; raw tables retained; no AI image generation','src/plot_diagnostics.py; src/build_analysis_evidence.py')]
steps=pd.DataFrame(stage_rows,columns=['step','operation','evidence','executed_result_or_boundary','code']);steps.to_csv(O/'computation_steps.csv',index=False)
input_path=__import__('audit_pipeline').DATA
(O/'raw_input_manifest.json').write_text(json.dumps({'file_name':input_path.name,'bytes':input_path.stat().st_size,'sha256':hashlib.file_digest(input_path.open('rb'),'sha256').hexdigest(),'level':'supplied exported analysis panel, not original satellite observations','raw_rows':len(raw),'raw_grids':int(raw.grid_id.nunique()),'duplicate_grid_year_keys':int(raw.duplicated(['grid_id','year']).sum()),'sample_selection':'valid_model_both_y','selected_grids_for_worked_examples':ids,'example_selection':'two unique grids per corrected context, pandas random_state20261020,20261021,20261022; illustrative, not inference sample'},indent=2))
htmlparts=['<!doctype html><meta charset="utf-8"><title>Reproducibility evidence guide</title><style>body{max-width:1150px;margin:35px auto;font:16px Arial;line-height:1.65;padding:0 24px}table{border-collapse:collapse;font-size:12px}td,th{border:1px solid #ddd;padding:8px;vertical-align:top}img{max-width:100%}h1,h2{color:#111}</style><h1>Evidence from source tables to quantitative figures</h1><p>This guide documents the calculations linking source tables, event definitions, model estimates and quantitative figures. Source coverage and inference limits are reported alongside the results.</p>',steps.to_html(index=False),'<h2>Every main risk estimate can be recomputed from integer counts</h2>',counts.to_html(index=False,float_format=lambda v:f'{v:.6f}'),'<h2>Confounding reduction and remaining limits</h2>',pd.DataFrame(protection).to_html(index=False),'<h2>Additional agricultural exclusion results</h2>',non[['comparison','response','AFB','ci_low','ci_high','n_total','n_model','n_grid']].to_html(index=False,float_format=lambda v:f'{v:.6f}')]
for m in maps:htmlparts+=['<h2>'+m['figure']+'</h2><img src="figures/'+m['figure']+'.png"><p>'+html.escape(m['source'])+'</p>']
(O/'evidence_guide.html').write_text(''.join(htmlparts),encoding='utf8')
print('ANALYSIS EVIDENCE SAVED',len(examples),'worked grid-years; balance\n',balance.to_string(index=False),flush=True)
