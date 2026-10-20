from pathlib import Path
import ee,pandas as pd,numpy as np,json,time
R=Path(__file__).resolve().parents[1];out=R/'experiments/E01_spatial_support';ee.Initialize(project=__import__('os').environ.get('EE_PROJECT','ee-abao9474'));ee.data.setDeadline(55000);d=pd.read_parquet(__import__('os').environ.get('AQUATIC_PANEL',str(R/'data/raw/china_aquatic_buffering_panel_v2.parquet')));g=d[d.valid_model_both_y].drop_duplicates('grid_id')[['grid_id','lon','lat']].copy();g['zone']=np.floor((g.lon+180)/6).astype(int)+1;gsw=ee.Image('JRC/GSW1_4/GlobalSurfaceWater');occ=gsw.select('occurrence').unmask(0);sea=gsw.select('seasonality').unmask(0);rows=[]
for zone,s in g.groupby('zone'):
 proj=ee.Projection('EPSG:'+str(32600+int(zone))).atScale(30);bands=[]
 for name,mask in [('near_perm_1km_corrected',occ.gte(80).And(sea.gte(10))),('near_seasonal_1km_corrected',sea.gte(1).And(sea.lt(10)))]:
  # A64-pixel (1.92km) search is sufficient for a1km boolean threshold. No far distances are reported.
  dist=mask.reproject(proj).fastDistanceTransform(64,'pixels','squared_euclidean').sqrt().multiply(30).reproject(proj);bands.append(dist.lte(1000).rename(name))
 image=ee.Image.cat(bands)
 for start in range(0,len(s),2000):
  sub=s.iloc[start:start+2000];fc=ee.FeatureCollection([ee.Feature(ee.Geometry.Point([float(r.lon),float(r.lat)]),{'grid_id':r.grid_id,'zone':int(zone)}) for r in sub.itertuples()]);v=image.sampleRegions(collection=fc,properties=['grid_id','zone'],scale=30,projection=proj,tileScale=1,geometries=False).getInfo();rows.extend([f['properties'] for f in v['features']]);pd.DataFrame(rows).to_csv(out/'national_corrected_1km_contexts_received.csv',index=False);print('received',len(rows),'of',len(g),'zone',zone,flush=True)
a=pd.DataFrame(rows);assert len(a)==len(g) and a.grid_id.nunique()==len(g);a.to_csv(out/'national_corrected_1km_contexts.csv',index=False);(out/'corrected_1km_method.json').write_text(json.dumps({'resolution_m':30,'projection':'UTM per6-degree zone, northern hemisphere','classification':'native JRC80%/10months; seasonal1<=months<10, nearest-neighbour metric projection','search_pixels':64,'search_m':1920,'reported_quantity':'boolean distance<=1000m only; no magnitude beyond search radius inferred','n_grid':len(a),'scope':'centroid context correction; within-cell distance separate'},indent=2));print('COMPLETE',len(a),flush=True)
