import gmpy2, time
n = gmpy2.mpz(586135293385)*gmpy2.mpz(10)**4336 + 1647
Q = n**47 - (n-1)**47
print("bits Q", Q.bit_length(), "bits Q-1", (Q-1).bit_length(), flush=True)
t=time.time(); print("BPSW (strong base-2 + strong Lucas-Selfridge):", gmpy2.is_bpsw_prp(Q), "%.0fs"%(time.time()-t), flush=True)
