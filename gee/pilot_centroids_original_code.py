from pathlib import Path
import ee,pandas as pd,time
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'experiments/E01_spatial_support';ee.Initialize(project=__import__('os').environ.get('EE_PROJECT','ee-abao9474'));ee.data.setDeadline(55000)
g=pd.read_csv(out/'pilot_selected_grids.csv');fc=ee.FeatureCollection([ee.Feature(ee.Geometry.Point([r.lon,r.lat]),{'grid_id':r.grid_id}) for r in g.itertuples()]);gsw=ee.Image('JRC/GSW1_4/GlobalSurfaceWater');occ=gsw.select('occurrence').unmask(0);sea=gsw.select('seasonality').unmask(0);proj=ee.Projection('EPSG:32650').atScale(30);rows=[]
for scenario,ot,st in [('original_code_main',80,10)]:
 p=occ.gte(ot).And(sea.gte(st)).reproject(proj);dist=p.fastDistanceTransform(1024,'pixels','squared_euclidean').sqrt().multiply(.03).rename('recomputed_centroid_km').reproject(proj);res=dist.sampleRegions(collection=fc,properties=['grid_id'],scale=30,projection=proj,geometries=False).getInfo()
 for f in res['features']:rows.append(dict(scenario=scenario,**f['properties']))
 print(scenario,len(res['features']),flush=True)
pd.DataFrame(rows).to_csv(out/'pilot_original_code_centroid_distances.csv',index=False)
