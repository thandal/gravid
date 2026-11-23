"""Summarize fixed-model prediction, coverage and exact input provenance."""
from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[3];BASE=ROOT/'results/timing/portal_arrival_times';W=BASE/'validation'
d=json.loads((W/'results.json').read_text());checks={}
for name,label in [('J0437-4715','adopted')]:
 ep=np.loadtxt(W/name/(label+'_O4b_epochs.txt'));p={l.split()[0]:l.split()[1] for l in (W/name/(label+'.par')).read_text().splitlines() if len(l.split())>1};t=ep[:,0]
 phase=np.remainder((t-float(p.get('TASC',p.get('T0'))))/float(p['PB']),1);gaps=np.diff(np.r_[np.sort(phase),np.sort(phase)[0]+1])
 checks[name]={'O4b_epochs':len(t),'largest_time_gap_days':float(np.diff(np.sort(t)).max()),'largest_sampled_orbital_phase_gap_cycles':float(gaps.max()),'O4b_sampled_loss_percent':100*d[name][label]['O4b']['raw_equal_epoch_power_loss'],'O4ab_sampled_loss_percent':100*d[name][label]['O4ab']['raw_equal_epoch_power_loss']}
 for log in (W/name).glob('*.log'):
  s=log.read_text()
  for flag in ['[CLK3]','[CLK4]','[CLK5]','[CLK6]','[CLK8]','[CLK9]']:
   if flag in s:raise ValueError(f'Clock fallback: {log}')
checks['clock_fallback_scan']='passed'
(W/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
files=list(W.rglob('*.par'))+list(W.rglob('*.tim'))+list(W.rglob('*.res'))+list(W.rglob('*.json'))+list((BASE/'clock').glob('*.clk'))+list(Path(__file__).parent.glob('*.py'))+[ROOT/'scripts/timing/orbital_offset_fit/evaluate_timing.C']
(W/'manifest.json').write_text(json.dumps({str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},indent=2)+'\n')
print(json.dumps(checks,indent=2))
