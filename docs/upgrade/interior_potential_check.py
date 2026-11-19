# Independent check of the formulas of docs/upgrade/interior_potential.md (not the engine's code).
#   python3 docs/upgrade/interior_potential_check.py      (numpy, scipy, mpmath; about a minute)
# alpha + 12C, l = 1, KTF Woods-Saxon; scipy DOP853 + mpmath Coulomb functions.
import numpy as np, mpmath as mp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
HBARC=197.3269804; ALPHA=1/137.035999084; U=931.49410242
m1,m2=4.00260325413,12.0; mu=m1*m2/(m1+m2)*U; Z=2*6; e2=ALPHA*HBARC
V0,R,A,RC,L=45.891,3.55,0.60,3.55,1
K=2*mu/HBARC**2
def V(r):
    vc = Z*e2/r if r>=RC else Z*e2/(2*RC)*(3-(r/RC)**2)
    return vc - V0/(1+np.exp((r-R)/A)) + L*(L+1)/(K*r*r)
def rhs(r,y,E): return [y[1], K*(V(r)-E)*y[0]]
def f_reg(E,a):                       # a u'(a)/u(a), regular solution outward from r0
    r0=1e-3; s=solve_ivp(rhs,(r0,a),[r0**(L+1),(L+1)*r0**L],args=(E,),method='DOP853',rtol=1e-12,atol=1e-300)
    return a*s.y[1,-1]/s.y[0,-1]
def coul(E,r):
    k=np.sqrt(K*E); eta=Z*e2*mu/(HBARC**2*k); rho=k*r
    F=float(mp.coulombf(L,eta,rho)); G=float(mp.coulombg(L,eta,rho))
    dF=float(mp.diff(lambda x: mp.coulombf(L,eta,x),rho)); dG=float(mp.diff(lambda x: mp.coulombg(L,eta,x),rho))
    return k,F,G,dF,dG
def L_out(E,a,rm=15.0):               # a H'(a)/H(a), H = G+iF at rm, inward; and phi = arg H(a)
    k,F,G,dF,dG=coul(E,rm)
    y0=[G,k*dG,F,k*dF]
    def r4(r,y): c=K*(V(r)-E); return [y[1],c*y[0],y[3],c*y[2]]
    s=solve_ivp(r4,(rm,a),y0,method='DOP853',rtol=1e-12,atol=1e-300)
    H=s.y[0,-1]+1j*s.y[2,-1]; dH=s.y[1,-1]+1j*s.y[3,-1]
    return a*dH/H, np.angle(H)
def delta_direct(E,rm=15.0):          # nuclear phase shift of the potential, matched to Coulomb at rm
    f=f_reg(E,rm); k,F,G,dF,dG=coul(E,rm); rho=k*rm
    return np.arctan2(f*F-rho*dF, rho*dG-f*G)
print("1. no levels: phase from U_pot = exp(-2i phi) (f - L*)/(f - L) against the direct potential phase (deg)")
for a in (5.35,4.0,6.5):
    worst=0
    for E in (0.5,1.0,2.0,2.4,2.6,3.5,5.0):
        f=f_reg(E,a); Lc,phi=L_out(E,a); Upot=np.exp(-2j*phi)*(f-np.conj(Lc))/(f-Lc)
        d=np.angle(Upot)/2; dd=delta_direct(E)
        diff=(d-dd+np.pi/2)%np.pi-np.pi/2; worst=max(worst,abs(diff))
    print(f"   a = {a} fm: worst |difference| {np.degrees(worst):.2e} deg; |U_pot|-1 = {abs(abs(Upot)-1):.1e}")
a=5.35; b=0.0
E0=brentq(lambda E: f_reg(E,a)-b,1.0,3.0,xtol=1e-10); h=1e-4
g2=-1/((f_reg(E0+h,a)-f_reg(E0-h,a))/(2*h)); wl=HBARC**2/(mu*a*a)
print(f"2. R_pot pole (a = 5.35 fm, b = 0): E = {E0:.4f} MeV, gamma^2 = {g2:.4f} MeV, theta^2 = gamma^2/(3 hbar^2/(2 mu a^2)) = {g2/(1.5*wl):.3f}  [KTF: 1.926 MeV, 0.94]")
print("3. one level + R_pot: U from R = R_pot + g^2/(El - E) against the level-matrix form with primed functions")
El,gam2=3.2,0.15
for E in (1.0,1.926,2.5,3.1,4.0):       # the level matrix's boundary condition is taken equal to b here
    f=f_reg(E,a); Lc,phi=L_out(E,a)
    # direct: R = R_pot + the pole term, both with the boundary constant b
    Rr=1/(f-b)+gam2/(El-E); L0=Lc-b
    Ud=np.exp(-2j*phi)*(1-np.conj(L0)*Rr)/(1-L0*Rr)
    Lp=b+(Lc-b)*(f-b)/(f-Lc); Pp=Lp.imag; Sp=Lp.real
    phip=phi+np.angle((f-Lc)/(f-b))
    Ainv=(El-E)-gam2*(Sp-b+1j*Pp)
    Um=np.exp(-2j*phip)*(1+2j*Pp*gam2/Ainv)
    print(f"   E = {E:5.3f}: |U_direct - U_primed| = {abs(Ud-Um):.1e}   P = {Lc.imag:.4f}  P' = {Pp:.4f}  S' = {Sp:+.4f}")
# df/dE identity
E=2.0; f1=(f_reg(E+h,a)-f_reg(E-h,a))/(2*h)
r0=1e-3; s=solve_ivp(rhs,(r0,a),[r0**(L+1),(L+1)*r0**L],args=(E,),method='DOP853',rtol=1e-12,atol=1e-300,dense_output=True)
rr=np.linspace(r0,a,20001); uu=s.sol(rr)[0]; integ=np.trapezoid(uu*uu,rr)
print(f"4. df/dE = -(2 mu a/hbar^2) int_0^a u^2 dr / u(a)^2: finite difference {f1:.6f}, integral {-K*a*integ/uu[-1]**2:.6f} MeV^-1")
