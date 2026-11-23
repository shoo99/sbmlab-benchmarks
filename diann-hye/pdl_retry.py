"""Parallel ranged downloader: PRIDE throttles each connection (~0.35 MB/s), so split each file. Same as pdl.py with 2000 retries per part.

Fixed 2026-11-18: the earlier version wrote whatever PRIDE sent into the file. PRIDE sometimes answers a range request
with a ~200-byte XML error body ("<Error><Code>ConnectionClosedException ...") in place of data; the file size still
matched, so the damage went unnoticed (2 of 6 Orbitrap DDA files and 5 of 6 SWATH .wiff.scan files of these benchmarks).
This version rejects any response that is not 206 with the requested Content-Range, never writes an error body, and scans
the finished file for one. PRIDE gives no checksums for PXD028735, so also compare two independent downloads (cmp)."""
import sys, os, threading, urllib.request, time

CH = 10  # connections per file


def get(url, start, end, path, idx, out):
    f = None
    part = f'{path}.part{idx}'
    for _ in range(2000):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have >= end - start + 1:
            break
        try:
            req = urllib.request.Request(url, headers={'Range': f'bytes={start + have}-{end}'})
            with urllib.request.urlopen(req, timeout=60) as r:
                # 2026-11-18: PRIDE sometimes answers a range request with an XML error body
                # ("<Error><Code>Connection..."); the old version appended it to the file and the size still matched.
                cr = r.headers.get('Content-Range', '')
                ct = r.headers.get('Content-Type', '')
                if r.status != 206 or not cr.startswith(f'bytes {start + have}-') or 'xml' in ct or 'html' in ct:
                    raise IOError(f'bad response {r.status} {cr!r} {ct!r}')
                f = open(part, 'ab')
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    if b'<Error><Code>' in b:
                        raise IOError('error body inside the stream')
                    f.write(b)
                f.close()
        except Exception:
            try:
                f.close()
            except Exception:
                pass
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
    bad = False
    with open(path, 'rb') as g:   # last line of defence: any HTTP error body that slipped through
        tail = b''
        while True:
            b = g.read(1 << 24)
            if not b:
                break
            if b'<Error><Code>' in tail + b:
                bad = True
                break
            tail = b[-32:]
    ok = os.path.getsize(path) == size and not bad
    print(('OK ' if ok else 'BAD ') + path, os.path.getsize(path), size, 'error-body' if bad else '', flush=True)


items = [l.split() for l in open(sys.argv[1]) if l.strip()]
ths = [threading.Thread(target=fetch, args=(u, int(s), os.path.join(sys.argv[2], os.path.basename(u)))) for u, s in items]
[t.start() for t in ths]
[t.join() for t in ths]
print('RAWDONE', flush=True)
