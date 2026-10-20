from pathlib import Path
import sys,json,time
import ee,pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'experiments/E01_spatial_support'
ee.Initialize(project=__import__('os').environ.get('EE_PROJECT','ee-abao9474'));ee.data.setDeadline(55000)
d=pd.read_parquet(__import__('os').environ.get('AQUATIC_PANEL',str(Path(__file__).resolve().parents[1]/'data/raw/china_aquatic_buffering_panel_v2.parquet')));g=d[d.valid_model_both_y].drop_duplicates('grid_id');g=g[g.lon.between(115.5,117.3)&g.lat.between(28.3,30.0)].sample(n=24,random_state=20261020).sort_values('grid_id');g[['grid_id','lon','lat','dist_perm_water_km']].to_csv(OUT/'pilot_selected_grids.csv',index=False)
# Classify on the source JRC grid, then nearest-neighbour reproject binary masks to metric 30 m UTM.
gsw=ee.Image('JRC/GSW1_4/GlobalSurfaceWater');occ=gsw.select('occurrence').unmask(0);sea=gsw.select('seasonality').unmask(0);wc=ee.ImageCollection('ESA/WorldCover/v200').first().select('Map');veg=wc.remap([10,20,30,40,90,95,100],[1]*7,0).eq(1).And(occ.eq(0));proj=ee.Projection('EPSG:32650').atScale(30);area=ee.Image.pixelArea()
rows=[]
for scenario,ot,st in [('main',80,9),('moderate',70,6),('strict',90,12)]:
 perm=occ.gte(ot).And(sea.gte(st)).reproject(proj)
 # Search radius 30.72 km; output masks missing neighbours. No distances are silently zero-filled.
 dist=perm.fastDistanceTransform(1024,'pixels','squared_euclidean').sqrt().multiply(30).rename('dist_m').reproject(proj)
 valid=dist.lt(30720); va=area.updateMask(veg).rename('vegetation_area'); measured=va.updateMask(valid).rename('distance_supported_area')
 bands=ee.Image.cat([area.rename('cell_area'),area.multiply(perm).rename('permanent_area'),va,measured,dist.multiply(measured).rename('distance_area_sum')])
 for align,dx in [('original',0),('shift_NE_0p05',.05)]:
  feats=[]
  for r in g.itertuples():
   x,y=float(r.lon)+dx,float(r.lat)+dx;geom=ee.Geometry.Rectangle([x-.05,y-.05,x+.05,y+.05],proj='EPSG:4326',geodesic=False);feats.append(ee.Feature(geom,{'grid_id':r.grid_id,'scenario':scenario,'alignment':align,'center_lon':x,'center_lat':y}))
  fc=ee.FeatureCollection(feats);t=time.time()
  try:
   response=bands.reduceRegions(collection=fc,reducer=ee.Reducer.sum(),crs='EPSG:32650',scale=30,tileScale=8,maxPixelsPerRegion=3000000).getInfo(); props=[f['properties'] for f in response['features']]
   for p in props:
    p['water_fraction']=p['permanent_area']/p['cell_area'];p['mean_vegetation_distance_km']=p['distance_area_sum']/p['distance_supported_area']/1000 if p.get('distance_supported_area',0)>0 else None;p['vegetation_distance_coverage']=p.get('distance_supported_area',0)/p['vegetation_area'] if p.get('vegetation_area',0)>0 else None;p['seconds']=time.time()-t
   rows.extend(props);pd.DataFrame(rows).to_csv(OUT/'pilot_within_cell_metrics.csv',index=False);print(scenario,align,'done',len(props),flush=True)
  except Exception as e:
   with (OUT/'gee_errors.txt').open('a',encoding='utf8') as f:f.write(scenario+' '+align+': '+str(e)+'\n')
   print(scenario,align,type(e).__name__,str(e)[:300],flush=True)
(OUT/'pilot_method.json').write_text(json.dumps({'status':'pilot only; not national E01 completion','region':'Poyang Lake surroundings; seeded sample of 24 existing analysis grids','grid_size_degrees':.1,'alignments':[0,.05],'classification':'JRC occurrence and seasonality native source grid, nearest neighbour metric reproject 30m','thresholds':{'moderate':[70,6],'main':[80,9],'strict':[90,12]},'threshold_warning':'moderate and strict are declared new sensitivity settings; original settings unverified','vegetation_support':'ESA WorldCover v200 2021 classes 10,20,30,40,90,95,100; exclude any JRC occurrence >0','metric_crs':'EPSG:32650','distance_search_m':30720,'mask_note':'mean distance only valid when vegetation_distance_coverage approx 1','alignment_limit':'shifted metric supports only; vegetation/climate responses not reaggregated so this is not a completed alignment sensitivity of AFB','seed':20261020},indent=2))

