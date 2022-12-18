#dr.gharib85@gmail.com
# correcltioan functon is essntial component for HEOM  construction  
# even ith ikeda nice work for non exponetial form of (relatted item ) correlation function 
# still restrictrd to matsubara or pade decomposition (there are alternative or even improvemnt in particular lo temperatu -> baycentric ,,fano )
# from my point of view   fiting correlation function is best handy implemtation
# but howefer analytical forn of correlatin function is needded see bofin -heom -ex-1d
#here  i start to obtain correlation function numerically then fitting
#code here is just litt modification of parts obtaind (qutip+pyrho)
import numpy as np
from scipy import integrate
from numpy import pi,tan
from mpmath import mp





class specttraldensity():
    def __init__(self,slist):
        sd=slist[0]
        if sd=="deby":
            self.l=slist[1]
            self.wc=slist[2]
            self.temp=slist[3]#add tem
            self.J=self.deby
            self.omega_inf=20*self.wc
        elif sd=="ohmic":
            self.a=slist[1]
            self.wc=slist[2]
            self.temp=slist[3]
            if len(slist)>4:
                self.s=slist[4]
            else:
                self.s=1
            self.J=self.ohmic
            self.omega_inf=100*self.wc #i dont know how 
            
            
            
        else:
            print("spectral density",sd,"not fount")
            raise SystemExit
    def deby(self,omega):
        w=abs(omega)
        jw=2*self.l*self.wc*w/(w**2+self.wc**2)
        return jw*(omega>=0)-jw*(omega<0)
    def ohmic(self,omega):
        w=abs(omega)
        jw=(self.a)*(w**self.s)*(self.wc**(1-self.s))*np.exp(-w/self.wc)
        return jw*(omega>=0)-jw*(omega<0)

    def re_bath_corr(self,omega):
        def coth(x):
            return 1.0/np.tanh(x)
        beta=1/(self.temp)
        omega += 1e-14
        n_omega=.5*(coth(beta*omega/2)-1.0)
        return self.J(omega)*(n_omega+1)
    def bath_corr(self, t):
        re_Ct = integrate.quad(self.re_bath_corr,
                                -self.omega_inf, self.omega_inf,
                                limit=1000, weight='cos', wvar=t)
        im_Ct = integrate.quad(self.re_bath_corr, 
                                -self.omega_inf, self.omega_inf,
                                limit=1000, weight='sin', wvar=t)
        re_Ct, im_Ct = re_Ct[0], -im_Ct[0]
        return (1.0/np.pi)*(re_Ct + 1j*im_Ct)
    def ohmic_corr_analytical(self,t):
        beta=1/self.temp
        alpha=self.a
        wc=self.wc
        s=self.s
        corr = (1/pi) * alpha * wc**(1 - s) * beta**(-(s + 1)) * mp.gamma(s + 1)
        z1_u = (1 + beta * wc - 1.0j * wc * t) / (beta * wc)
        z2_u = (1 + 1.0j * wc * t) / (beta * wc)
        return np.array([
        complex(corr * (mp.zeta(s + 1, u1) + mp.zeta(s + 1, u2)))
        for u1, u2 in zip(z1_u, z2_u)
    ], dtype=np.complex128)


ohmic=specttraldensity(["ohmic",1,2,3])
t=np.linspace(0,100,1000)
corr_analytical=ohmic.ohmic_corr_analytical(t)

corr_numerical=np.array([ohmic.bath_corr(i) for i in t])
import matplotlib.pyplot as plt
c=corr_analytical-corr_numerical
plt.plot(t,np.real(c))
plt.show()

"""plt.plot(t,np.real(corr_analytical))
plt.plot(t,np.imag(corr_analytical))
plt.show()
plt.plot(t,np.real(corr_numerical))
plt.plot(t,np.imag(corr_numerical))
plt.show()
"""
