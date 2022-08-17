# Coulomb phase = arg Gamma(l+1-i/k)

from scipy.special import gamma
import cmath

k=1.0

for l in range(0,4):
    z=l+1-(1j)/k
    print(("{0:.5f} {1:.5f}").format(z.real, z.imag))
    v=gamma(z)
    t=cmath.phase(v)
    print(("l={0:1d}, arg G(l+1-i/k)={1:8.5f}").format(l, t))
    
