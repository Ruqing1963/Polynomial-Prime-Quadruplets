import gmpy2, time, hashlib, sys
sys.set_int_max_str_digits(0)
n = gmpy2.mpz(586135293385)*gmpy2.mpz(10)**4336 + 1647
Q = n**47 - (n-1)**47
print("sha256(decimal)", hashlib.sha256(str(Q).encode()).hexdigest(), flush=True)
t=time.time()
print("GMP strong PRP base 11:", gmpy2.is_strong_prp(Q, 11), "%.0fs"%(time.time()-t), flush=True)
