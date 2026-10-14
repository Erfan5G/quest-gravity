"""Download verified source archives without adding large inputs to Git."""
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile

BASE = 'https://www.geosciencebc.com/i/project_data/QUESTdata/GBCReport2008-8/'
FILES = ['GBCReport2008-8_Gravity_Data.zip', 'Gravity_Technical_Report.pdf']


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while b := f.read(1024*1024):
            h.update(b)
    return h.hexdigest()


def main():
    root = Path(__file__).resolve().parents[2]
    raw = root/'data/raw'
    raw.mkdir(parents=True, exist_ok=True)
    expected = json.loads((root/'docs/source_checksums.json').read_text())
    for name in FILES:
        p = raw/name
        if not p.exists():
            print('Downloading', name, flush=True)
            tmp = p.with_suffix(p.suffix+'.part')
            with urlopen(BASE+name, timeout=120) as response, tmp.open('wb') as out:
                while b := response.read(1024*1024):
                    out.write(b)
            tmp.replace(p)
        if sha256(p) != expected[name]:
            raise ValueError(f'Source checksum changed: {name}; review before processing')
    target = raw/'release'
    target.mkdir(exist_ok=True)
    with ZipFile(raw/FILES[0]) as z:
        for name in z.namelist():
            destination = (target/name).resolve()
            if not destination.is_relative_to(target.resolve()):
                raise ValueError('Unsafe archive path')
        z.extractall(target)
    print('Verified sources and extracted release.')


if __name__ == '__main__':
    main()
