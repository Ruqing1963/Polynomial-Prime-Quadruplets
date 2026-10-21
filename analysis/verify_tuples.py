import json, gmpy2, time
Q = lambda n: gmpy2.mpz(n) ** 47 - gmpy2.mpz(n - 1) ** 47
S = json.load(open("quad_starts_from_zip.json"))
paper15 = [23159557,117309848,136584738,218787064,411784485,423600750,523331634,640399031,987980498,1163461515,1370439187,1643105964,1691581855,1975860550,1996430175]
st = set(S) | set(paper15)
q5 = sorted(n for n in st if n + 1 in st and n < 8e11)
t = time.time()
bad5 = [n for n in q5 if not all(gmpy2.is_bpsw_prp(Q(n + j)) for j in range(5))]
w20 = sorted(n for n in st if n > 1e19)
bad20 = [n for n in w20 if not all(gmpy2.is_bpsw_prp(Q(n + j)) for j in range(4))]
base = 4797791281515011781318
w1000 = [base + o for o in (720056643, 7642492197, 9120900759, 12970778218, 15262442996, 16351036649,
                            17840013419, 28374845686, 34066697501, 36297085036, 37001546293)]
bad1000 = [n for n in w1000 if not all(gmpy2.is_bpsw_prp(Q(n + j)) for j in range(4))]
# is any 1000-digit quadruplet part of a quintuplet?  check Q(n-1), Q(n+4)
ext = [(n, bool(gmpy2.is_bpsw_prp(Q(n - 1))), bool(gmpy2.is_bpsw_prp(Q(n + 4)))) for n in w1000]
digits = [gmpy2.num_digits(Q(n), 10) for n in (w20[0], w1000[0])]
print(json.dumps({"quintuplets": len(q5), "quintuplet_failures": bad5, "1e20_quads": len(w20), "1e20_failures": bad20,
                  "1000d_quads": len(w1000), "1000d_failures": bad1000, "1000d_neighbors_prime": [e for e in ext if e[1] or e[2]],
                  "digits(1e20, 1000d)": digits, "time_s": round(time.time() - t)}, indent=1))
json.dump({"quintuplets": q5, "window_1e20": w20, "window_1000d": w1000}, open("verified_tuples.json", "w"), indent=1)
