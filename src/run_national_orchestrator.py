from pathlib import Path

if (Path(__file__).resolve().parents[1]/'experiments/E01_spatial_support/heavy_pause.flag').exists():
 raise SystemExit('Heavy jobs paused by user; explicit authorization is required before removing the pause flag.')
import subprocess,sys,time,json,datetime
BASE=Path(__file__).resolve().parents[1];root=BASE;out=root/'experiments/E01_spatial_support';finished=set();start=time.time();log=root/'logs/national_analysis_orchestrator.log'
def note(msg):
 with log.open('a',encoding='utf8') as f:f.write(datetime.datetime.now(datetime.timezone.utc).isoformat()+' '+msg+'\n')
 print(msg,flush=True)
while time.time()-start<7200:
 proc=subprocess.run([sys.executable,'-X','utf8',str(BASE/'src'/('collect_national.py'))],capture_output=True,text=True,encoding='utf8');note('collector exit '+str(proc.returncode)+' '+proc.stdout.strip())
 for name,file in [('centroid','national_corrected_centroids.csv'),('within_cell','national_within_cell_metrics.csv')]:
  if name not in finished and (out/file).exists():
   run=subprocess.run([sys.executable,'-X','utf8',str(BASE/'src'/('analyze_corrected.py'))],capture_output=True,text=True,encoding='utf8');(root/'logs'/('corrected_'+name+'_analysis.log')).write_text(run.stdout+'\n'+run.stderr,encoding='utf8');note(name+' analysis exit '+str(run.returncode));
   if run.returncode==0:
    finished.add(name)
    subprocess.run([sys.executable,'-X','utf8',str(BASE/'src'/('build_report.py'))],capture_output=True,text=True,encoding='utf8')
 if len(finished)==2:note('Both complete corrected datasets analyzed');break
 time.sleep(30)
else:note('Two-hour orchestration bound reached; inspect cloud task manifests before resuming')
