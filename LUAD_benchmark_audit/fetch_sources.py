"""Download the paper supplement and rebuild the public patient-level matrix.

Run in the unpacked audit folder after installing requirements.txt.
The published matrix snapshot is included. Modern remote releases can change.
"""
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile
from io import BytesIO
import os
import re
import time

def main():
    root = Path(__file__).resolve().parent
    os.chdir(root)
    (root / 'supplement').mkdir(exist_ok=True)
    url = 'https://www.mdpi.com/article/10.3390/ijms25073661/s1'
    for attempt in range(3):
        try:
            with urlopen(url, timeout=65) as response:
                archive = response.read()
            break
        except OSError:
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)
    with ZipFile(BytesIO(archive)) as zf:
        sheets = [n for n in zf.namelist() if re.fullmatch(r'TableS[123] [^/]+\.xlsx', n)]
        if len(sheets) != 3:
            raise ValueError('Expected three supplementary workbooks')
        for name in sheets:
            (root / 'supplement' / name).write_bytes(zf.read(name))
    from cohort import prepare
    from xfetch import run
    from build_matrix import build
    prepare()
    run()
    build()

if __name__ == '__main__':
    main()
