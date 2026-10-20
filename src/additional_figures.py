from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_pipeline import *
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
colors={'EVI':'#25805c','LAI':'#7857a3','combined':'#b63d49'};mapping=pd.read_csv(ROOT/'experiments/E05_figure_audit/figure_source_map.csv');new=[]
def finish(fig,name,source):
 fig.savefig(ROOT/'figures'/str(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(ROOT/'figures'/str(name+'.pdf'),bbox_inches='tight');plt.close(fig);new.append(dict(figure=name,source=source,status='recomputed; original distance audit context'))
r=pd.read_csv(ROOT/'experiments/E05_figure_audit/regional_AFB_grid_cluster.csv');f=pd.read_csv(ROOT/'experiments/E00_original_reproduction/original_regional_AFB_rerun_with_omissions.csv');order=['Tibetan_Plateau','Southwest','North_China','Northeast','Northwest_Arid','Other','East_China','South_China'];fig,axs=plt.subplots(1,3,figsize=(16,6))
for ax,e in zip(axs,colors):
 t=r[(r.response==e)&(r.comparison=='other')].set_index('region').reindex(order);labels=[]
 for i,region in enumerate(order):
  row=t.loc[region];counts=f[(f.region==region)&(f.response==e)].iloc[0];label=region.replace('_',' ')+'\nrows near/other: '+str(counts.n_near)+'/'+str(counts.n_other)
  if row.insufficient_other_vs_permanent:
   label='\u25cb '+label;ax.text(.03,i,'insufficient support',transform=ax.get_yaxis_transform(),va='center',fontsize=7)
  else:ax.errorbar(row.AFB,i,xerr=[[row.AFB-row.ci_low],[row.ci_high-row.AFB]],fmt='o',color=colors[e],capsize=3)
  labels.append(label)
 ax.set_yticks(range(len(order)),labels,fontsize=7);ax.invert_yaxis();ax.set_ylim(len(order)-.5,-.5);ax.axvline(0,color='grey',lw=.7);ax.set_title(e);ax.set_xlabel('AFB: other minus permanent')
fig.suptitle('All recovered macroregions shown; open circle in label = inadequate contrast support\n300 grid-cluster draws; original distance inputs are under audit');fig.tight_layout();finish(fig,'Fig8_all_regions_with_support','E05/regional_AFB_grid_cluster.csv; E00/original_regional_AFB_rerun_with_omissions.csv; recovered rectangular macroregions')
# Corrected same-raster pilot diagnostic; physical distances are not silently treated as equivalent.
a=pd.read_csv(ROOT/'experiments/E01_spatial_support/pilot_centroid_comparison.csv');orig=pd.read_csv(ROOT/'experiments/E01_spatial_support/original_vs_recomputed_centroid_audit.csv');fig,axs=plt.subplots(1,3,figsize=(14,4));axs[0].scatter(orig.dist_perm_water_km,orig.recomputed_centroid_km,c='#3678a8');lim=max(orig.recomputed_centroid_km.max(),orig.dist_perm_water_km.max());axs[0].plot([0,lim],[0,lim],ls='--',color='grey');axs[0].set_xlabel('Original table centroid distance (km label)');axs[0].set_ylabel('Recomputed 30m UTM distance (km)');axs[0].set_title('Original-distance audit (80/9 pilot)')
for scenario,t in a.groupby('scenario'):axs[1].scatter(t.recomputed_centroid_km,t.mean_vegetation_distance_km,label=scenario,s=15)
axs[1].plot([0,24],[0,24],ls='--',color='grey');axs[1].set_xlabel('Same-raster centroid distance (km)');axs[1].set_ylabel('Within-cell vegetation mean distance (km)');axs[1].legend(fontsize=7);axs[1].set_title('Within-cell support: 24 grids')
b=pd.read_csv(ROOT/'experiments/E01_spatial_support/pilot_alignment_metric_comparison.csv')
for scenario,t in b.groupby('scenario'):axs[2].scatter(t.water_fraction_original,t.water_fraction_change,label=scenario,s=15)
axs[2].axhline(0,color='grey',ls='--');axs[2].set_xlabel('Original-cell permanent-water fraction');axs[2].set_ylabel('Fraction change after NE0.05 degree shift');axs[2].set_title('Metric-only alignment sensitivity');fig.tight_layout();finish(fig,'SFig7_within_cell_pilot','E01/pilot_centroid_comparison.csv; original_vs_recomputed_centroid_audit.csv; pilot_alignment_metric_comparison.csv; no shifted response reaggregation')
a=pd.read_csv(ROOT/'experiments/E02_confounding/harmonized/M0_AFB.csv');b=pd.read_csv(ROOT/'experiments/E02_confounding/harmonized/M2_AFB.csv');c=pd.read_csv(ROOT/'experiments/E02_confounding/harmonized/common_landscape_AFB.csv');fig,axs=plt.subplots(1,3,figsize=(14,4))
for ax,e in zip(axs,colors):
 rows=[a[(a.response==e)&(a.comparison=='seasonal')].iloc[0],b[(b.response==e)&(b.comparison=='seasonal')].iloc[0],c[c.response==e].iloc[0]]
 for i,r in enumerate(rows):ax.errorbar(r.AFB,i,xerr=[[r.AFB-r.ci_low],[r.ci_high-r.AFB]],fmt='o',color=colors[e],capsize=3)
 ax.set_yticks(range(3),['M0 / M1','Harmonized M2','Common landscape']);ax.axvline(0,color='grey',lw=.7);ax.set_title(e);ax.set_xlabel('AFB: seasonal minus permanent')
fig.suptitle('Original distance grouping; reanalysis after distance correction still required');fig.tight_layout();finish(fig,'SFig8_harmonized_and_landscape','E02/harmonized/M0_AFB.csv; M2_AFB.csv; common_landscape_AFB.csv')
save(pd.concat([mapping,pd.DataFrame(new)],ignore_index=True),'experiments/E05_figure_audit/figure_source_map.csv')
print('additional figures complete')
