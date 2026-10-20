from pathlib import Path
import json,ee,collections
out=Path(__file__).resolve().parents[1]/'experiments/E01_spatial_support';ee.Initialize(project=__import__('os').environ.get('EE_PROJECT','ee-abao9474'))
for name in ['national_corrected_centroid','national_within_cell']:
 m=json.loads((out/(name+'_tasks.json')).read_text());r=ee.data.getTaskStatus([t['task_id'] for t in m['tasks']]);(out/(name+'_task_status.json')).write_text(json.dumps(r,indent=2));print(name,dict(collections.Counter(x['state'] for x in r)),[{'description':x.get('description'),'error':x.get('error_message')} for x in r if x['state']=='FAILED'])
