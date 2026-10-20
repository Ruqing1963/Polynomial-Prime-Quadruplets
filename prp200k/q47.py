#!/usr/bin/env python3
"""
q47.py -- sieve + PFGW driver for probable primes  Q(n) = n^47 - (n-1)^47.

Number theory used
  * Q(n) is odd, Q(n) == 1 (mod 47), and every prime factor p of Q(n) has p == 1 (mod 94).
  * For such p:  p | Q(n)  <=>  n == 1/(1-u) (mod p)  for one of the 46 nontrivial 47th
    roots of unity u mod p.  Roots come in pairs r, 1-r because Q(n) = Q(1-n).
  So the sieve only touches primes 283, 659, 941, ... and removes 46 classes for each.

Candidates are written as n = a0*10^e + t with a 12-digit a0, so PFGW lines stay short.

Typical run (200,000-digit search):
  python q47.py init   --digits 200000 --length 150000 --state q47.npz
  python q47.py sieve  --state q47.npz --pmax 1e12 --prp-minutes 10   # resumable; rerun with larger --pmax to deepen
  python q47.py status --state q47.npz
  python q47.py export --state q47.npz --split 8 --outdir pfgw_run
  python q47.py selftest                                            # brute-force check of the sieve

Requires numpy.  Python >= 3.8.
"""
import argparse
import math
import os
import sys
import time
from array import array
from multiprocessing import Pool

import numpy as np

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)  # n has ~4348 digits, above Python's default 4300 limit

EXP = 47
MOD = 2 * EXP  # every prime factor of Q(n) is 1 mod 94
NROOTS = EXP - 1
EULER_GAMMA = 0.5772156649015329


def Q(n):
    return n**EXP - (n - 1) ** EXP


def iroot(x, k):
    """floor(x ** (1/k)) for integers x >= 1."""
    r = 1 << -(-x.bit_length() // k)
    while True:
        s = ((k - 1) * r + x // r ** (k - 1)) // k
        if s >= r:
            return r
        r = s


def parse_num(s):
    return int(float(s)) if any(c in s for c in "eE.") else int(s)


# ---------------------------------------------------------------- primes p == 1 (mod 94)

def small_primes(limit):
    s = np.ones(limit + 1, dtype=bool)
    s[:2] = False
    for i in range(2, math.isqrt(limit) + 1):
        if s[i]:
            s[i * i :: i] = False
    return np.nonzero(s)[0]


def primes_1_mod_94(lo, hi, sp):
    """Primes p == 1 (mod 94) with lo <= p < hi.  sp must contain all primes <= sqrt(hi)."""
    k_lo = max(1, -(-(lo - 1) // MOD))
    k_hi = -(-(hi - 1) // MOD)  # exclusive
    if k_hi <= k_lo:
        return np.zeros(0, dtype=np.int64)
    p_max = MOD * (k_hi - 1) + 1
    seg = np.ones(k_hi - k_lo, dtype=bool)
    for q in sp.tolist():
        if q * q > p_max:
            break
        if q == 2 or q == EXP:  # 94k+1 is never divisible by 2 or 47
            continue
        k0 = -pow(MOD, -1, q) % q
        start = k_lo + (k0 - k_lo) % q
        if MOD * start + 1 == q:
            start += q
        seg[start - k_lo :: q] = False
    return (np.nonzero(seg)[0].astype(np.int64) + k_lo) * MOD + 1


def q_roots(p):
    """The 46 residues n (mod p) with p | Q(n), for a prime p == 1 (mod 94)."""
    e = (p - 1) // EXP
    g = 2
    while (w := pow(g, e, p)) == 1:
        g += 1
    half = NROOTS // 2
    d, pre = [], []
    u = acc = 1
    for _ in range(half):  # u = w^1 .. w^23; the inverses w^-j give the partner roots 1 - r
        u = u * w % p
        x = (1 - u) % p
        d.append(x)
        acc = acc * x % p
        pre.append(acc)
    inv = pow(acc, -1, p)  # batch inversion of all 23 values (1 - u)
    r = [0] * half
    for i in range(half - 1, 0, -1):
        r[i] = inv * pre[i - 1] % p
        inv = inv * d[i] % p
    r[0] = inv
    return r + [(1 - x) % p for x in r]


# ---------------------------------------------------------------- state

def start_of(st):
    return st["a0"] * 10 ** st["e"] + st["r0"]


def save_state(path, st):
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        np.savez_compressed(
            f,
            alive=np.packbits(st["alive"]),
            length=np.int64(len(st["alive"])),
            meta=np.array(
                [str(st[k]) for k in ("digits", "a0", "e", "step", "r0", "pmax", "wheel")]
            ),
        )
    os.replace(tmp, path)


def load_state(path):
    with np.load(path) as z:
        L = int(z["length"])
        alive = np.unpackbits(z["alive"])[:L].astype(bool)
        meta = [int(x) for x in z["meta"]]
    st = dict(zip(("digits", "a0", "e", "step", "r0", "pmax", "wheel"), meta))
    st["alive"] = alive
    return st


def prime_probability(st):
    """Heuristic P(Q(n) prime | Q(n) has no prime factor <= pmax) = e^gamma * ln(pmax) / ln Q."""
    return math.exp(EULER_GAMMA) * math.log(st["pmax"]) / (st["digits"] * math.log(10))


def report(st):
    S = int(np.count_nonzero(st["alive"]))
    L = len(st["alive"])
    pp = prime_probability(st)
    print(f"  sieved to p <= {st['pmax']:.3e}: {S} of {L} candidates left ({S / L:.4%})")
    print(f"  P(PRP) per candidate ~ 1/{1 / pp:,.0f};  expected PRPs in this set ~ {S * pp:.3f}")
    return S


# ---------------------------------------------------------------- sieve

_W = {}


def _init_worker(start, step, length, sp):
    _W.update(start=start, step=step, length=length, sp=sp)


def _sieve_chunk(bounds):
    lo, hi = bounds
    start, step, L = _W["start"], _W["step"], _W["length"]
    small, hits = [], array("q")
    for p in primes_1_mod_94(lo, hi, _W["sp"]).tolist():
        if step % p == 0:  # wheel prime: residue class already chosen to avoid it
            continue
        a = start % p
        minv = pow(step, -1, p) if step != 1 else 1
        offs = [(x - a) * minv % p for x in q_roots(p)]
        if p < L:
            small.append((p, offs))
        else:
            hits.extend(o for o in offs if o < L)
    return hi, small, hits.tobytes()


def run_sieve(st, pmax, workers=1, chunk=200_000_000, prp_minutes=None, path=None, save_every=300):
    alive = st["alive"]
    L = len(alive)
    lo = st["pmax"] + 1
    if pmax < lo:
        return
    chunks = [(a, min(a + chunk, pmax + 1)) for a in range(lo, pmax + 1, chunk)]
    initargs = (start_of(st), st["step"], L, small_primes(math.isqrt(pmax) + 1))
    pool = None
    if workers > 1:
        pool = Pool(workers, _init_worker, initargs)
        results = pool.imap(_sieve_chunk, chunks)
    else:
        _init_worker(*initargs)
        results = map(_sieve_chunk, chunks)

    t0 = t_last = t_saved = time.time()
    left_last = int(np.count_nonzero(alive))
    try:
        for hi, small, hits in results:
            for p, offs in small:
                for o in offs:
                    alive[o::p] = False
            if hits:
                alive[np.frombuffer(hits, dtype=np.int64)] = False
            st["pmax"] = hi - 1
            now = time.time()
            if path and (now - t_last > 10 or hi > pmax):
                left = int(np.count_nonzero(alive))
                line = f"p <= {st['pmax']:.4e}  left {left}  removed {left_last - left}  [{now - t0:.0f}s]"
                if prp_minutes and left_last > left:
                    core_s = (now - t_last) * workers
                    gain = (left_last - left) * prp_minutes * 60 / core_s
                    line += f"  PRP-time saved / sieve-time spent ~ {gain:.1f}x"
                print(line, flush=True)
                t_last, left_last = now, left
                if now - t_saved > save_every:
                    save_state(path, st)
                    t_saved = now
    except KeyboardInterrupt:
        print("interrupted; saving state")
    finally:
        if pool:
            pool.terminate()
        if path:
            save_state(path, st)


# ---------------------------------------------------------------- commands

def cmd_init(args):
    D, L = args.digits, args.length
    T = 10 ** (D - 1)
    n = iroot(T // EXP, EXP - 1)
    while Q(n) < T:
        n += 1
    nd = len(str(n))
    e = max(nd - args.lead_digits, 0)
    a0 = -(-n // 10**e)
    base = a0 * 10**e

    step, r0, wheel = 1, 0, []
    if args.wheel:
        sp = small_primes(10**5)
        wheel = primes_1_mod_94(2, 10**8, sp)[: args.wheel].tolist()
        step = math.prod(wheel)
        r0 = -base % step  # n == 0 (mod p) gives Q(n) == 1 (mod p): never divisible
    n_last = base + r0 + step * (L - 1)
    if n_last > iroot(10**D // EXP, EXP - 1):
        sys.exit("range too long: Q(n) would exceed the requested digit count")

    st = dict(digits=D, a0=a0, e=e, step=step, r0=r0, pmax=1, wheel=len(wheel),
              alive=np.ones(L, dtype=bool))
    save_state(args.state, st)
    print(f"n has {nd} digits; n = {a0}*10^{e} + t")
    if wheel:
        print(f"wheel primes {wheel}: t = {r0} + {step}*k  (n == 0 mod each wheel prime)")
    print(f"{L} candidates, Q(n) has {D} digits. State -> {args.state}")


def cmd_sieve(args):
    st = load_state(args.state)
    pmax = parse_num(args.pmax)
    if pmax <= st["pmax"]:
        print(f"already sieved to {st['pmax']:.3e}")
    else:
        run_sieve(st, pmax, args.workers, parse_num(args.chunk), args.prp_minutes, args.state)
    report(st)


def cmd_status(args):
    st = load_state(args.state)
    print(f"Q(n) digits {st['digits']}, n = {st['a0']}*10^{st['e']} + t, wheel primes {st['wheel']}")
    report(st)


SH = """#!/bin/sh
# One PFGW instance per work_* directory (each needs its own pfgw.ini / pfgw.log).
# PRPs found are appended to work_*/pfgw.log.
PFGW=${PFGW:-pfgw64}
cd "$(dirname "$0")"
for d in work_*; do
  (cd "$d" && nohup nice "$PFGW" -f0 -l"pfgw_out.txt" candidates.txt > /dev/null 2>&1 &)
done
echo "started; check work_*/pfgw.log for PRP results"
"""

BAT = """@echo off
rem One PFGW instance per work_* directory (each needs its own pfgw.ini / pfgw.log).
rem PRPs found are appended to work_*\\pfgw.log.  Edit PFGW if pfgw64.exe is not on PATH.
set PFGW=pfgw64.exe
cd /d "%~dp0"
for /d %%D in (work_*) do (
  start "pfgw %%D" /low /d "%%D" %PFGW% -f0 -l"pfgw_out.txt" candidates.txt
)
echo started; check work_*\\pfgw.log for PRP results
"""


def cmd_export(args):
    st = load_state(args.state)
    S = report(st)
    ks = np.nonzero(st["alive"])[0]
    if args.limit:
        ks = ks[: args.limit]
    head = f"{st['a0']}*10^{st['e']}"
    os.makedirs(args.outdir, exist_ok=True)
    for i, part in enumerate(np.array_split(ks, args.split)):
        d = os.path.join(args.outdir, f"work_{i:02d}")
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "candidates.txt"), "w", newline="\n") as f:
            if args.format == "abc":
                f.write(f"ABC ({head}+$a)^{EXP}-({head}+$a-1)^{EXP}\n")
            for k in part.tolist():
                t = st["r0"] + st["step"] * k
                if args.format == "abc":
                    f.write(f"{t}\n")
                else:
                    f.write(f"({head}+{t})^{EXP}-({head}+{t}-1)^{EXP}\n")
    for name, text in (("run_pfgw.sh", SH), ("run_pfgw.bat", BAT)):
        with open(os.path.join(args.outdir, name), "w", newline="\r\n" if name.endswith("bat") else "\n") as f:
            f.write(text)
    print(f"wrote {len(ks)} of {S} candidates into {args.split} work dirs under {args.outdir}/")


def cmd_selftest(args):
    import random

    rng = random.Random(1)
    sp = small_primes(10**4)
    ps = primes_1_mod_94(2, 10**6, sp).tolist()
    for p in rng.sample(ps, 50):
        r = q_roots(p)
        assert len(set(r)) == NROOTS and all(Q(x) % p == 0 for x in r), p
    print("roots: ok (46 distinct zeros of Q mod p for 50 random p)")

    D, L, pmax = 1000, 4000, 30000
    allp = math.prod(small_primes(pmax).tolist())  # ALL primes, also checks the 1 mod 94 theorem
    for w in (0, 3):
        a = argparse.Namespace(digits=D, length=L, lead_digits=12, wheel=w, state="_selftest.npz")
        cmd_init(a)
        st = load_state(a.state)
        os.remove(a.state)
        run_sieve(st, pmax, workers=1, chunk=7000)
        n0 = start_of(st)
        brute = np.array([math.gcd(Q(n0 + st["step"] * k), allp) == 1 for k in range(L)])
        assert (brute == st["alive"]).all(), "sieve disagrees with brute force"
        print(f"sieve: ok (wheel={w}, {int(brute.sum())} survivors match gcd brute force)")
        try:
            import gmpy2
        except ImportError:
            continue
        found = sum(gmpy2.is_prime(Q(n0 + st["step"] * k), 2) for k in np.nonzero(st["alive"])[0].tolist())
        exp_ = int(st["alive"].sum()) * prime_probability(st)
        print(f"  {found} PRPs among survivors at {D} digits; heuristic expected {exp_:.1f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("init", help="choose the n range and create a state file")
    a.add_argument("--digits", type=int, default=200000, help="decimal digits of Q(n)")
    a.add_argument("--length", type=int, default=150_000,
                   help="number of n values to sieve (150k gives ~2.8 expected PRPs at 200k digits)")
    a.add_argument("--lead-digits", type=int, default=12, help="digits of a0 in n = a0*10^e + t")
    a.add_argument("--wheel", type=int, default=0,
                   help="fix n == 0 mod the first W allowed primes (same effect as sieving them)")
    a.add_argument("--state", default="q47.npz")
    a.set_defaults(func=cmd_init)

    a = sub.add_parser("sieve", help="sieve (or deepen the sieve) up to --pmax")
    a.add_argument("--state", default="q47.npz")
    a.add_argument("--pmax", default="1e9")
    a.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    a.add_argument("--chunk", default="2e8", help="p-range per work unit")
    a.add_argument("--prp-minutes", type=float, default=None,
                   help="measured PFGW minutes per test on one core; prints whether sieving still pays")
    a.set_defaults(func=cmd_sieve)

    a = sub.add_parser("status")
    a.add_argument("--state", default="q47.npz")
    a.set_defaults(func=cmd_status)

    a = sub.add_parser("export", help="write PFGW input files and run scripts")
    a.add_argument("--state", default="q47.npz")
    a.add_argument("--outdir", default="pfgw_run")
    a.add_argument("--split", type=int, default=1, help="number of work dirs (one PFGW per core)")
    a.add_argument("--format", choices=("abc", "expr"), default="abc")
    a.add_argument("--limit", type=int, default=0, help="export only the first N survivors")
    a.set_defaults(func=cmd_export)

    a = sub.add_parser("selftest")
    a.set_defaults(func=cmd_selftest)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
