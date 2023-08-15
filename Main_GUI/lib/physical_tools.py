# Tools for physical calculations

import math
import numpy as np

def calcOperatorCoeff(Y_coeff, polarization_i, theta, phi):
    # conversion of angle from degree to radian
    theta=math.radians(theta)
    phi=math.radians(phi)
    if polarization_i==0:
        # linear
        Y_coeff[2]=-math.sqrt(2.0*math.pi/3.0)*math.sin(theta)*(math.cos(phi)-math.sin(phi)*1j)
        Y_coeff[1]= math.sqrt(4.0*math.pi/3.0)*math.cos(theta)
        Y_coeff[0]= math.sqrt(2.0*math.pi/3.0)*math.sin(theta)*(math.cos(phi)+math.sin(phi)*1j)
    elif polarization_i==1:
        # right circular
        Y_coeff[2]=math.sqrt(2.0*math.pi/3.0)*(1+math.cos(theta))*(math.cos(phi)-math.sin(phi)*1j)
        Y_coeff[1]=math.sqrt(4.0*math.pi/3.0)*math.sin(theta)
        Y_coeff[0]=math.sqrt(2.0*math.pi/3.0)*(1-math.cos(theta))*(math.cos(phi)+math.sin(phi)*1j)
    elif polarization_i==2:
        # left circular
        Y_coeff[2]= math.sqrt(2.0*math.pi/3.0)*(-1+math.cos(theta))*(math.cos(phi)-math.sin(phi)*1j)
        Y_coeff[1]= math.sqrt(4.0*math.pi/3.0)*math.sin(theta)
        Y_coeff[0]=-math.sqrt(2.0*math.pi/3.0)*( 1-math.cos(theta))*(math.cos(phi)+math.sin(phi)*1j)
    else:
        print("Error: invalid polarization")
        return

    print(("Coefficients for Y_{{1,+1}}: {0:7.3f} + {1:7.3f}i").format(Y_coeff[2].real, Y_coeff[2].imag))
    print(("Coefficients for Y_{{1,+0}}: {0:7.3f} + {1:7.3f}i").format(Y_coeff[1].real, Y_coeff[1].imag))
    print(("Coefficients for Y_{{1,-1}}: {0:7.3f} + {1:7.3f}i").format(Y_coeff[0].real, Y_coeff[0].imag))

# see OpenMX/source/AngularF.c for the order of the spherical harmonics                            
def convertLCAO_p(px, py, pz, LCAO):
    # px   1/2 sqrt(3/pi) sin(T)cos(P)
    # py   1/2 sqrt(3/pi) sin(T)sin(P)
    # pz   1/2 sqrt(3/pi) cos(T)

    # pm1  1/2 sqrt(3/2pi) sin(T)(cos(P)-i*sin(P)) = (px-i*py)/sqrt(2)
    # pp0  1/2 sqrt(3/pi)  cos(T)                  = pz
    # pp1 -1/2 sqrt(3/2pi) sin(T)(cos(P)+i*sin(P)) = -(px+i*py)/sqrt(2)

    # therefore,
    # phi = LCAO(px) *px  + LCAO(py) *py
    # phi = LCAO(pm1)*pm1 + LCAO(pp1)*pp1
    #   LCAO(pm1)- LCAO(pp1)=sqrt(2)LCAO(px)
    # -iLCAO(pm1)-iLCAO(pp1)=sqrt(2)LCAO(py)
    # LCAO(pm1)= (LCAO(px)+iLCAO(py))/sqrt(2)
    # LCAO(pp1)=-(LCAO(px)-iLCAO(py))/sqrt(2)
    
    LCAO[0]=(px+1j*py)/math.sqrt(2)
    LCAO[1]=pz
    LCAO[2]=-(px-1j*py)/math.sqrt(2)

def convertLCAO_d(d3z2r2, dx2y2, dxy, dxz, dyz, LCAO):
    # d3z2r2 1/4 sqrt(5/pi)   (3cos^2(T)-1)
    # dx2y2  1/4 sqrt(15/pi)  sin^2(T)cos(2P)
    # xy     1/4 sqrt(15/pi)  sin^2(T)sin(2P)
    # xz     1/2 sqrt(15/pi)  sin(T)cos(T)cos(P)
    # yz     1/2 sqrt(15/pi)  sin(T)cos(T)sin(P)

    # dm2    1/4 sqrt(15/2pi) sin^2(T)(cos(2P)-i*sin(2P))   = (dx2y2-i*xy)/sqrt(2)
    # dm1    1/2 sqrt(15/2pi) sin(T)cos(T)(cos(P)-i*sin(P)) = (xz-i*yz)/sqrt(2)
    # dp0    1/4 sqrt(5/pi)   (3cos^2(T)-1)                 = d3z2r2
    # dp1   -1/2 sqrt(15/2pi) sin(T)cos(T)(cos(P)+i*sin(P)) = -(xz+i*yz)/sqrt(2)
    # dp2    1/4 sqrt(15/2pi) sin^2(T)(cos(2P)+i*sin(2P))   = (dx2y2+i*xy)/sqrt(2)
    LCAO[0]=(dx2y2+1j*dxy)/math.sqrt(2)
    LCAO[1]=(dxz+1j*dyz)/math.sqrt(2)
    LCAO[2]=d3z2r2
    LCAO[3]=-(dxz-1j*dyz)/math.sqrt(2)
    LCAO[4]=(dx2y2-1j*dxy)/math.sqrt(2)

def convertLCAO_f(f5z23r2, f5xy2xr2, f5yz2yr2, fzx2zy2, fxyz, fx33xy2, f3yx2y3, LCAO):
    # f5z23r2  1/4 sqrt(7/pi)   (5cos^2(T)-3)cos(T)
    # f5xy2xr2 1/8 sqrt(42/pi)  (5cos^2(T)-1)sin(T)cos(P)
    # f5yz2yr2 1/8 sqrt(42/pi)  (5cos^2(T)-1)sin(T)sin(P)
    # fzx2zy2  1/4 sqrt(105/pi) sin^2(T)cos(T)cos(2P)
    # fxyz     1/4 sqrt(105/pi) sin^2(T)cos(T)sin(2P)
    # fx33xy2  1/8 sqrt(70/pi)  sin^3(T)cos(3P)
    # f3yx2y3  1/8 sqrt(70/pi)  sin^3(T)sin(3P)

    # fm3  1/8 sqrt(35/pi)   sin^3(T)(cos(3P)-i*sin(3P))          = (fx33xy2-i*f3yx2y3)/sqrt(2)
    # fm2  1/4 sqrt(105/2pi) sin^2(T)cos(T)(cos(2P)-i*sin(2P))    = (fzx2zy2-i*fxyz)/sqrt(2)
    # fm1  1/8 sqrt(21/pi)   (5cos^2(T)-1)sin(T)(cos(P)-i*sin(P)) = (f5xy2xr2-i*f5yz2yr2)/sqrt(2)
    # fp0  1/4 sqrt(7/pi)    (5cos^2(T)-3)cos(T)                  = f5z23r2
    # fp1 -1/8 sqrt(21/pi)   (5cos^2(T)-1)sin(T)(cos(P)+i*sin(P)) = -(f5xy2xr2+i*f5yz2yr2)/sqrt(2)
    # fp2  1/4 sqrt(105/2pi) sin^2(T)cos(T)(cos(2P)+i*sin(2P))    = (fzx2zy2+i*fxyz)/sqrt(2)
    # fm3 -1/8 sqrt(35/pi)   sin^3(T)(cos(3P)+i*sin(3P))          = -(fx33xy2+i*f3yx2y3)/sqrt(2)
    LCAO[0]=(fx33xy2+1j*f3yx2y3)/math.sqrt(2)
    LCAO[1]=(fzx2zy2+1j*fxyz)/math.sqrt(2)
    LCAO[2]=(f5xy2xr2+1j*f5yz2yr2)/math.sqrt(2)
    LCAO[3]=f5z23r2
    LCAO[4]=-(f5xy2xr2-1j*f5yz2yr2)/math.sqrt(2)
    LCAO[5]=(fzx2zy2-1j*fxyz)/math.sqrt(2)
    LCAO[6]=-(fx33xy2-1j*f3yx2y3)/math.sqrt(2)

def spBessel(l, x):
    # spherical Bessel function j_l(x)
    # j_l(x)=(-1)^l x^l (1/x d/dx)^l sin(x)/x
    if l==0:
        # sin(x)/x
        return math.sin(x)/x
    elif l==1:
        # {sin(x)-x*cos(x)}/x^2
        return (math.sin(x)-x*math.cos(x))/x**2
    elif l==2:
        # {(3-x^2)sin(x)-3x*cos(x)}/x^3
        return ((3-x**2)*math.sin(x)-3*x*math.cos(x))/x**3
    elif l==3:
        # {(15-6x^2)sin(x)+(x^3-15x)cos(x)}/x^4
        return ((15-6*x**2)*math.sin(x)+(x**3-15*x)*math.cos(x))/x**4
    elif l==4:
        # {(x^4-45x^2+105)sin(x)+(10x^3-105x)cos(x)}/x^5
        return ((x**4-45*x**2+105)*math.sin(x)+(10*x**3-105*x)*math.cos(x))/x**5
    else:
        return 0

# no longer used
def radialIntegral(wfn1, wfn2, r, dr):
    # int wfn1*wfn2*r dr
    ret=0.0
    for i in range (0, r.shape[0]-1):
        ret+=wfn1[i]*wfn2[i]*r[i]*(r[i+1]-r[i])
    return ret

def Gaunt(lp, mp, l, m):
    if abs(mp)>lp or abs(m)>l:
        return 0

    if mp==m:
        if lp==l+1:
            return math.sqrt((3.0/(4.0*math.pi))*(((l+1.0)**2-m**2)/((2.0*l+3.0)*(2.0*l+1.0))))
        elif lp==l-1:
            return math.sqrt((3.0/(4.0*math.pi))*((l**2-m**2)/((2.0*l-1.0)*(2.0*l+1.0))))
        else:
            return 0
    elif mp==m+1:
        if lp==l+1:
            return math.sqrt((3.0/(4.0*math.pi))*(((l+m+2.0)*(l+m+1.0))/(2.0*(2.0*l+3.0)*(2.0*l+1.0))))
        elif lp==l-1:
            return -math.sqrt((3.0/(4.0*math.pi))*(((l-m)*(l-m-1.0))/(2.0*(2.0*l-1.0)*(2.0*l+1.0))))
        else:
            return 0
    elif mp==m-1:
        if lp==l+1:
            return math.sqrt((3.0/(4.0*math.pi))*(((l-m+2.0)*(l-m+1.0))/(2.0*(2.0*l+3.0)*(2.0*l+1.0))))
        elif lp==l-1:
            return -math.sqrt((3.0/(4.0*math.pi))*(((l+m)*(l+m-1.0))/(2.0*(2.0*l-1.0)*(2.0*l+1.0))))
        else:
            return 0
    else:
        return 0
    
def sphericalHarmonics(r):
    # r=[x, y, z]
    # x=r*sin(theta)*cos(phi)
    # y=r*sin(theta)*sin(phi)
    # sqrt(x^2+y^2)=r*sin(theta)
    # z=r*cos(theta)
    Ylm=np.zeros((5,11), dtype=complex)
    
    r_length=math.sqrt(np.inner(r, r))
    cosT=r[2]/r_length
    sinT=math.sqrt(r[0]**2+r[1]**2)/r_length
    cosP=1
    sinP=0
    if abs(sinT)>1e-5:
        cosP=r[0]/(r_length*sinT)
        sinP=r[1]/(r_length*sinT)

    cos2P=cosP**2-sinP**2
    sin2P=2*sinP*cosP

    cos3P=4*cosP**3-3*cosP
    sin3P=3*sinP-4*sinP**3

    cos4P=8*cosP**4-8*cosP**2+1
    sin4P=4*sinP*cosP*(2*cosP**2-1)

    # s(0): 1/2 1/sqrt(pi)
    Ylm[0][0]=1.0/(2.0*math.sqrt(math.pi))
    
    # p(-1): 1/2 sqrt(3/2pi) sin(T)(cos(P)-i*sin(P))
    Ylm[1][0]=(1.0/2.0)*math.sqrt(3.0/(2.0*math.pi))*sinT*(cosP-1j*sinP)
    # p(0): 1/2 sqrt(3/pi) cos(T)
    Ylm[1][1]=(1.0/2.0)*math.sqrt(3.0/math.pi)*cosT
    # p(1): -1/2 sqrt(3/2pi) sin(T)(cos(P)+i*sin(P))
    Ylm[1][2]=-(1.0/2.0)*math.sqrt(3.0/(2.0*math.pi))*sinT*(cosP+1j*sinP)

    # d(-2): 1/4 sqrt(15/2pi) sin^2(T)(cos(2P)-i*sin(2P))
    Ylm[2][0]=(1.0/4.0)*math.sqrt(15.0/(2.0*math.pi))*sinT**2*(cos2P-1j*sin2P)
    # d(-1): 1/2 sqrt(15/2pi) sin(T)cos(T)(cos(P)-i*sin(P))
    Ylm[2][1]=(1.0/2.0)*math.sqrt(15.0/(2.0*math.pi))*sinT*cosT*(cosP-1j*sinP)
    # d(0): 1/4 sqrt(5/pi) (3cos^2(T)-1)
    Ylm[2][2]=(1.0/4.0)*math.sqrt(5.0/math.pi)*(3*cosT**2-1.0)
    # d(1): -1/2 sqrt(15/2pi) sin(T)cos(T)(cos(P)+i*sin(P))
    Ylm[2][3]=-(1.0/2.0)*math.sqrt(15.0/(2.0*math.pi))*sinT*cosT*(cosP+1j*sinP)
    # d(2): 1/4 sqrt(15/2pi) sin^2(T)(cos(2P)+i*sin(2P))
    Ylm[2][4]=(1.0/4.0)*math.sqrt(15.0/(2.0*math.pi))*sinT**2*(cos2P+1j*sin2P)
    
    # f(-3): 1/8 sqrt(35/pi) sin^3(T)(cos(3P)-i*sin(3P))
    Ylm[3][0]=(1.0/8.0)*math.sqrt(35.0/math.pi)*sinT**3*(cos3P-1j*sin3P)
    # f(-2): 1/4 sqrt(105/2pi) sin^2(T)cos(T)(cos(2P)-i*sin(2P))
    Ylm[3][1]=(1.0/4.0)*math.sqrt(105.0/(2.0*math.pi))*sinT**2*cosT*(cos2P-1j*sin2P)
    # f(-1): 1/8 sqrt(21/pi) (5cos^2(T)-1)sin(T)(cos(P)-i*sin(P))
    Ylm[3][2]=(1.0/8.0)*math.sqrt(21.0/math.pi)*(5.0*cosT**2-1.0)*sinT*(cosP-1j*sinP)
    # f(0): 1/4 sqrt(7/pi) (5cos^2(T)-3)cos(T)
    Ylm[3][3]=(1.0/4.0)*math.sqrt(7.0/math.pi)*(5.0*cosT**2-3.0)*cosT
    # f(1): -1/8 sqrt(21/pi) (5cos^2(T)-1)sin(T)(cos(P)+i*sin(P))
    Ylm[3][4]=-(1.0/8.0)*math.sqrt(21.0/math.pi)*(5.0*cosT**2-1.0)*sinT*(cosP+1j*sinP)
    # f(2): 1/4 sqrt(105/2pi) sin^2(T)cos(T)(cos(2P)+i*sin(2P))
    Ylm[3][5]=(1.0/4.0)*math.sqrt(105.0/(2.0*math.pi))*sinT**2*cosT*(cos2P+1j*sin2P)
    # f(3): -1/8 sqrt(35/pi) sin^3(T)(cos(3P)+i*sin(3P))
    Ylm[3][6]=-(1.0/8.0)*math.sqrt(35.0/math.pi)*sinT**3*(cos3P+1j*sin3P)
    
    # g(-4): 3/16 sqrt(35/2pi) sin^4(T)(cos(4P)-i*sin(4P))
    Ylm[4][0]=(3.0/16.0)*math.sqrt(35.0/(2.0*math.pi))*sinT**4*(cos4P-1j*sin4P)
    # g(-3): 3/8 sqrt(35/pi) sin^3(T)cos(T)(cos(3P)-i*sin(3P))
    Ylm[4][1]=(3.0/8.0)*math.sqrt(35.0/math.pi)*sinT**3*cosT*(cos3P-1j*sin3P)
    # g(-2): 3/8 sqrt(5/2pi) (7cos^2(T)-1)sin^2(T)(cos(2P)-i*sin(2P))
    Ylm[4][2]=(3.0/8.0)*math.sqrt(5.0/(2.0*math.pi))*(7.0*cosT**2-1.0)*(cos2P-1j*sin2P)
    # g(-1): 3/8 sqrt(5/pi) (7cos^2(T)-3)sin(T)cos(T)(cos(P)-i*sin(P))
    Ylm[4][3]=(3.0/8.0)*math.sqrt(5.0/math.pi)*(7.0*cosT**2-3.0)*sinT*cosT*(cosP-1j*sinP)
    # g(0): 3/16 sqrt(1/pi) (35cos^4(T)-30cos^2(T)+3)
    Ylm[4][4]=(3.0/16.0)*math.sqrt(1.0/math.pi)*(35.0*cosT**4-30.0*cosT**2+3.0)
    # g(1): -3/8 sqrt(5/pi) (7cos^2(T)-3)sin(T)cos(T)(cos(P)+i*sin(P))
    Ylm[4][5]=-(3.0/8.0)*math.sqrt(5.0/math.pi)*(7.0*cosT**2-3.0)*sinT*cosT*(cosP+1j*sinP)
    # g(2): 3/8 sqrt(5/2pi) (7cos^2(T)-1)sin^2(T)(cos(2P)+i*sin(2P))
    Ylm[4][6]=(3.0/8.0)*math.sqrt(5.0/(2.0*math.pi))*(7.0*cosT**2-1.0)*(cos2P+1j*sin2P)
    # g(3): -3/8 sqrt(35/pi) sin^3(T)cos(T)(cos(3P)+i*sin(3P))
    Ylm[4][7]=-(3.0/8.0)*math.sqrt(35.0/math.pi)*sinT**3*cosT*(cos3P+1j*sin3P)
    # g(4): 3/16 sqrt(35/2pi) sin^4(T)(cos(4P)+i*sin(4P))
    Ylm[4][8]=(3.0/16.0)*math.sqrt(35.0/(2.0*math.pi))*sinT**4*(cos4P+1j*sin4P)

    return Ylm
