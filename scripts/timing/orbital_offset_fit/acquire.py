"""Download the public J0437 subsets; run from repository root with network access."""
from pathlib import Path
import json
import subprocess
import tarfile
import shutil
import hashlib

B = Path('results/timing/j0437_external_validation')
REV = 'e2806fcb2c238ec38dd80784dfcb1e6c82f54728'
def get(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(['curl', '-sSL', '--fail', '--max-time', '120', url, '-o', str(path)], check=True)

get(f'https://api.github.com/repos/inpta/InPTA.DR2/git/trees/{REV}?recursive=1', B/'inpta/tree.json')
for entry in json.loads((B/'inpta/tree.json').read_text())['tree']:
    path = entry['path']
    if entry['type'] == 'blob' and (path.startswith('J0437-4715/') or path == 'README.md') and not path.endswith('.DS_Store'):
        assert '..' not in Path(path).parts and not Path(path).is_absolute()
        get(f'https://raw.githubusercontent.com/inpta/InPTA.DR2/{REV}/{path}', B/'inpta'/path)
NOISE_REV = 'e2f9cc573ba193e227672ec2fbe18dbf9a17f820'
get(f'https://raw.githubusercontent.com/inpta/InPTA.DR2.NA/{NOISE_REV}/par/J0437-4715.nonesw.par', B/'inpta/J0437-noise.par')
get(f'https://raw.githubusercontent.com/inpta/InPTA.DR2.NA/{NOISE_REV}/README.md', B/'inpta/noise_README.md')
get('https://docs.datacentral.org.au/documents/52/partim.tar.gz', B/'meerkat/partim.tar.gz')
with tarfile.open(B/'meerkat/partim.tar.gz') as archive:
    for filename in ['partim/J0437-4715.par', 'partim/J0437-4715.tim']:
        dest = B/'meerkat'/filename
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(archive.extractfile(filename).read())
for name, files in [('inpta', sorted((B/'inpta/J0437-4715/tims').glob('*.tim'))),
                    ('meerkat', [B/'meerkat/partim/J0437-4715.tim'])]:
    lines = [line for file in files for line in file.read_text().splitlines()
             if len(line.split()) > 4 and line.split()[0] not in ['C', '#']]
    (B/name/'flattened.tim').write_text('FORMAT 1\nMODE 1\n'+'\n'.join(lines)+'\n')
(B/'clock').mkdir(exist_ok=True)
for path in Path('.conda/share/tempo2/clock').glob('*.clk'):
    shutil.copy2(path, B/'clock'/path.name)
for name in ['gmrt2gps.clk','gps2utc.clk','tai2tt_bipm2020.clk','tai2tt_bipm2023.clk']:
    get('https://raw.githubusercontent.com/ipta/pulsar-clock-corrections/main/T2runtime/clock/'+name, B/'clock'/name)
# Clock upstream is mutable; preserve the acquired bytes and their hashes.
manifest = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in B.rglob('*') if p.is_file()}
Path('scripts/timing/orbital_offset_fit/download_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
