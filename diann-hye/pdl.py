"""Parallel ranged downloader: PRIDE throttles each connection (~0.35 MB/s), so split each file."""
import sys, os, threading, urllib.request, time

CH = 10  # connections per file


def get(url, start, end, path, idx, out):
    part = f'{path}.part{idx}'
    for _ in range(60):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have >= end - start + 1:
            break
        try:
            req = urllib.request.Request(url, headers={'Range': f'bytes={start + have}-{end}'})
            with urllib.request.urlopen(req, timeout=60) as r, open(part, 'ab') as f:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    f.write(b)
        except Exception:
            time.sleep(20)
    out[idx] = os.path.exists(part) and os.path.getsize(part) >= end - start + 1


def fetch(url, size, path):
    if os.path.exists(path) and os.path.getsize(path) == size:
        print('HAVE', path, flush=True)
        return
    step = size // CH
    rng = [(i * step, size - 1 if i == CH - 1 else (i + 1) * step - 1) for i in range(CH)]
    out = {}
    ts = [threading.Thread(target=get, args=(url, a, b, path, i, out)) for i, (a, b) in enumerate(rng)]
    [t.start() for t in ts]
    [t.join() for t in ts]
    if not all(out.values()):
        print('FAIL', path, flush=True)
        return
    with open(path, 'wb') as f:
        for i in range(CH):
            with open(f'{path}.part{i}', 'rb') as p:
                while True:
                    b = p.read(1 << 24)
                    if not b:
                        break
                    f.write(b)
            os.remove(f'{path}.part{i}')
    print(('OK ' if os.path.getsize(path) == size else 'BAD ') + path, os.path.getsize(path), size, flush=True)


items = [l.split() for l in open(sys.argv[1]) if l.strip()]
ths = [threading.Thread(target=fetch, args=(u, int(s), os.path.join(sys.argv[2], os.path.basename(u)))) for u, s in items]
[t.start() for t in ths]
[t.join() for t in ths]
print('RAWDONE', flush=True)
