from pathlib import Path
import json,ast,hashlib,re
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parents[1]
checks=[]
def check(name,ok):
 checks.append(dict(check=name,pass_=bool(ok)))
 if not ok:raise AssertionError(name)
a=pd.read_csv(ROOT/'reviewer_evidence/risk_numerators_and_arithmetic.csv');b=pd.read_csv(ROOT/'experiments/E01_spatial_support/corrected_1km/context_fragility.csv').rename(columns={'group':'context'});m=a.merge(b,on=['context','response'],suffixes=('_fresh','_saved'))
check('all nine context risks reconstruct from integer event counts',np.allclose(m.risk_low_fresh,m.risk_low_saved,atol=1e-13) and np.allclose(m.risk_high_fresh,m.risk_high_saved,atol=1e-13))
check('all context fragilities reconstruct',np.allclose(m.fragility_fresh,m.fragility_saved,atol=1e-13))
a=pd.read_csv(ROOT/'reviewer_evidence/AFB_arithmetic.csv');b=pd.read_csv(ROOT/'experiments/E01_spatial_support/corrected_1km/AFB.csv');m=a.merge(b,on=['comparison','response'],suffixes=('_fresh','_saved'));check('six AFB arithmetic contrasts match',np.allclose(m.AFB_fresh,m.AFB_saved,atol=1e-13))
e=pd.read_csv(ROOT/'reviewer_evidence/worked_grid_year_examples.csv');check('worked examples have six unique grids',e.grid_id.nunique()==6);check('worked HDSI formula exact',abs(e.HDSI-e.HDSI_recomputed).max()<1e-13);check('worked EVI/LAI events and OR exact',((e.EVI==(e.late_summer_EVI_z<-.5))&(e.LAI==(e.late_summer_LAI_z<-.5))&(e.combined==((e.EVI==1)|(e.LAI==1)))).all())
d=pd.read_csv(ROOT/'experiments/E02_confounding/corrected_noncrop/model_diagnostics.csv');check('noncrop fits converge and full rank',d.converged.all() and (d.n_parameters==d['rank']).all())
c=json.loads((ROOT/'experiments/E02_confounding/corrected_noncrop/exclusion_config.json').read_text());check('noncrop has fixed national thresholds',np.allclose(c['HDSI_cutpoints_fixed_to_full_sample'],[-.2974858281280519,.2882174368150908]) and c['after']=={'n_rows':172463,'n_grids':36238})
check('heavy pause remains in place',(ROOT/'experiments/E01_spatial_support/heavy_pause.flag').exists())
check('no Word documents in public tree',not list(ROOT.rglob('*.docx')))
check('no manuscript working-draft authoring sources',not any((ROOT/'src'/n).exists() for n in ['revise_manuscript.py','build_balanced_revision.py','build_report.py']))
for p in ROOT.rglob('*.py'):ast.parse(p.read_text(encoding='utf-8-sig'))
check('all public Python source files parse',True)
readme=(ROOT/'README.md').read_text(encoding='utf8');broken=[]
for target in re.findall(r'\]\(([^)]+)\)',readme):
 if not target.startswith(('https:','http:','#')) and not (ROOT/target.split('#')[0]).exists():broken.append(target)
check('README relative links resolve',not broken)
figs=pd.read_csv(ROOT/'docs/figure_inventory.csv');check('34 distinct available figures are indexed',len(figs)==34 and all((ROOT/p).exists() for p in figs.path))
rows=[]
for p in ROOT.rglob('*'):
 if p.is_file() and '.git' not in p.parts and '__pycache__' not in p.parts and p.name!='public_manifest.csv':rows.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
pd.DataFrame(rows).to_csv(ROOT/'public_manifest.csv',index=False)
(ROOT/'logs/public_evidence_validation.json').write_text(json.dumps({'checks':checks,'scope':'public code and derived evidence only; no authoring drafts or raw large exports','broken_links':broken},indent=2));print(json.dumps(checks,indent=2))
