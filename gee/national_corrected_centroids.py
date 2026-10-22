from pathlib import Path

if (Path(__file__).resolve().parents[1]/'experiments/E01_spatial_support/heavy_pause.flag').exists():
 raise SystemExit('Heavy jobs paused by user; explicit authorization is required before removing the pause flag.')
import ee,pandas as pd,numpy as np,json,datetime
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'experiments/E01_spatial_support';ee.Initialize(project=__import__('os').environ.get('EE_PROJECT','ee-abao9474'));g=pd.read_parquet(__import__('os').environ.get('AQUATIC_PANEL',str(Path(__file__).resolve().parents[1]/'data/raw/china_aquatic_buffering_panel_v2.parquet')));g=g[g.valid_model_both_y].drop_duplicates('grid_id')[['grid_id','lon','lat']].copy();g['zone']=np.floor((g.lon+180)/6).astype(int)+1
manifest=OUT/'national_corrected_centroid_tasks.json'
if manifest.exists():raise RuntimeError('Task manifest exists; do not submit duplicates')
gsw=ee.Image('JRC/GSW1_4/GlobalSurfaceWater');occ=gsw.select('occurrence').unmask(0);sea=gsw.select('seasonality').unmask(0);rows=[]
for zone,s in g.groupby('zone'):
 proj=ee.Projection('EPSG:'+str(32600+int(zone))).atScale(30);bands=[]
 for name,mask in [('dist_perm_corrected_km',occ.gte(80).And(sea.gte(10))),('dist_seasonal_corrected_km',sea.gte(1).And(sea.lt(10)))]:
  native=mask.reproject(proj);dist=native.fastDistanceTransform(2048,'pixels','squared_euclidean').sqrt().multiply(.03).rename(name).reproject(proj);bands.append(dist);bands.append(dist.lt(61.44).rename(name+'_in_range'))
 image=ee.Image.cat(bands);fc=ee.FeatureCollection([ee.Feature(ee.Geometry.Point([float(r.lon),float(r.lat)]),{'grid_id':r.grid_id,'zone':int(zone)}) for r in s.itertuples()]);samples=image.sampleRegions(collection=fc,properties=['grid_id','zone'],scale=30,projection=proj,tileScale=4,geometries=False);desc='E01_corrected_centroids_20261020_UTM'+str(zone);task=ee.batch.Export.table.toDrive(collection=samples,description=desc,folder='aquatic_spatial_metrics',fileNamePrefix=desc,fileFormat='CSV');task.start();rows.append(dict(zone=int(zone),n_grid=len(s),task_id=task.id,description=desc,crs='EPSG:'+str(32600+int(zone))));manifest.write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'full-sample correction of original distance algorithm; 30m UTM per 6-degree zone; no within-cell support yet','occurrence':80,'seasonality':10,'seasonal_definition':'1 <= seasonality <10','search_radius_km':61.44,'range_flag':'dist <61.44; farther distances are censored, not valid exact distances','folder':'aquatic_spatial_metrics','tasks':rows},indent=2));print(rows[-1],flush=True)
