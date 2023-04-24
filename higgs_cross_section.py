import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

mass = np.linspace(0, 50E-3, 1001)

alpha = (1)**2/(4*np.pi)

F = 6.64
Z = 74

m_e = 0.510E-3

alpha_h = (m_e**2)*(1.166E-5)*np.sqrt(2)/(4*np.pi)

A = (2*(alpha**2)*alpha_h*Z)/(m_e**2)

def B(m):
    return m**2/m_e**2

xsecs = []

for m in mass:
    b = B(m)
    f = lambda z : F*A*z*((1+2*(1-z)/3*z**2)*b)/(1+b*(1-z)/z**2)**2
    xsecs.append(integrate.quad(f, 0, 1)[0])

plt.scatter(mass*1000, xsecs)
plt.xlabel('Higgs boson mass [MeV]')
plt.ylabel('Cross section [b]')
plt.show()
   
nums = [] 
for x in xsecs:
    nums.append(x*2E12)

plt.scatter(mass*1000, nums)
plt.xlabel('Higgs boson mass [MeV]')
plt.ylabel('Number of events produced on target')
plt.show()


def tau_inverse(m):
    return 0.5*alpha_h*m*(1 - 4*m_e**2/m**2)**(3/2)


def gamma(m):
    return energy/m

def momentum(m):
    return np.sqrt(energy**2 - m**2)
energy = 1.6
taus = []
gammas = []
momentums = []

for m in mass:
    taus.append(tau_inverse(m))
    gammas.append(gamma(m))
    momentums.append(momentum(m))

def beta(gamma, p, m):
    return p/(m*gamma)

betas = []
for i, m in enumerate(mass):
    betas.append(beta(gammas[i], momentums[i], m))

def exponent(m, tau, beta, gamma):
    return -2*0.65*tau/(beta*gamma)

def acc_exp(m, tau, beta, gamma):
    return -4*tau/(beta*gamma)

exponents= []
acc_exps = []

for i, m in enumerate(mass):
    exponents.append(exponent(m, taus[i], betas[i], gammas[i]))
    acc_exps.append(acc_exp(m, taus[i], betas[i], gammas[i]))

survival_prob = []
acceptance = []
for i,m in enumerate(mass):
    survival_prob.append(np.exp(exponents[i]))
    acceptance.append(np.exp(acc_exps[i]))

surviving_num = []
accepted = []
for i, j in enumerate(nums):
    surviving_num.append(j*survival_prob[i])
    accepted.append(j*survival_prob[i]*acceptance[i])

plt.scatter(mass*1000, surviving_num)
plt.xlabel('Higgs boson mass [MeV]')
plt.ylabel('Number of Higgs bosons reaching detector')
plt.show()

plt.scatter(mass*1000, accepted)
plt.xlabel('Higgs boson mass [MeV]')
plt.ylabel('Number of decays')
plt.show()




                                            




    
