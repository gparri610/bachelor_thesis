import numpy as np
import matplotlib.pyplot as plt
import random


#constants and parameters
pb_convert = 3.894*10**8 #pb per Gev**(-2)
stw = 0.222246
Qf1 = 2./3.
Af1 = 1./2.
Vf1 = 1./2. - 4./3.*stw
Qf2 = -1./3.
Af2 = -1./2.
Vf2 = -1./2. - 2./3.*stw
E = 91.1876 #GeV
s = E**2 #GeV^2
Mz = 91.188 #Gev
Tz = 2.4414 #Gev
Vf_m = -1./2. + 2.*stw
Af_m = -1./2.
alpha = 1/(132.507)
Gf = 1.16639*10**(-5) #Gev^(-2)

kappa = np.sqrt(2)*Gf*Mz**2/(4*np.pi*alpha)
chi1 = kappa*s*(s-Mz**2)/((s-Mz**2)**2 + Tz**2*Mz**2)
chi2 = kappa**2*s**2/((s-Mz**2)**2 + Tz**2*Mz**2)

A01 = Qf1**2 -2*Qf1*Vf1*Vf_m*chi1 + (Af_m**2 + Vf_m**2)*(Af1**2 + Vf1**2)*chi2
A11 = -4*Qf1*Af1*Af_m*chi1 + 8*Af_m*Vf_m*Af1*Vf1*chi2

A02 = Qf2**2 -2*Qf2*Vf_m*Vf2*chi1 + (Af_m**2 + Vf_m**2)*(Af2**2 + Vf2**2)*chi2
A12 = -4*Qf2*Af2*Af_m*chi1 + 8*Af_m*Vf_m*Af2*Vf2*chi2


#differential cross section

def diff_sigma1(costh):

    return (2 * np.pi)*alpha**2*(A01*(1+costh**2)+ A11*costh)/(4*s)

# set the seed for random numbers
# random.random() will give us random numbers using the Mersenne Twister algorithm
# see: https://docs.python.org/2/library/random.html

def diff_sigma2(costh):

    return (2 * np.pi)*alpha**2*(A02*(1+costh**2)+ A12*costh)/(4*s)


seed = 12342
random.seed(seed)






# range [-1,1]

delta = 2


N = 1000000

sum_W = 0
sum_W_sq = 0

W_max1 = 0

array1 = []

#Monte Carlo (2.3)

for i in range(N):

    #creiamo valori a caso di cos(theta)
    costh_i = -1 + random.random()*delta

    #definiamo il weigth Wi
    W_i = delta*diff_sigma1(costh_i)
    array1.append(W_i)

    #sommiamo tutti i weight
    sum_W = sum_W + W_i

    #sommiamo i quadrati
    sum_W_sq = sum_W_sq + W_i**2

    if W_i > W_max1:
        W_max1 = W_i



    # Monte Carlo (2.3)

W_max2 = 0
array2 = []

for i in range(N):

        # creiamo valori a caso di cos(theta)
    costh_i = -1 + random.random() * delta

        # definiamo il weigth Wi
    W_i = delta * diff_sigma2(costh_i)
    array2.append(W_i)

        # sommiamo tutti i weight
    sum_W = sum_W + W_i

        # sommiamo i quadrati
    sum_W_sq = sum_W_sq + W_i ** 2

    if W_i > W_max2:
        W_max2 = W_i


#cross section (2.4)

#sigma = sum_W/N

#varianza (2.5)

#var = sum_W_sq/N - (sum_W/N)**2

#standard deviation

#err_MC = np.sqrt(var/N)

#cross section finale

#print('total cross section =', sigma * pb_convert, '+-', err_MC * pb_convert, 'pb')
#print("W_max =", W_max, "costh =", costh_max)


#############################################
#############################################
#############################################

#hit or miss method (2.3)

Nevents = 10000
uwe1 = []


for j in range(N):

    costh_j = -1 + random.random()*delta

    dsj = diff_sigma1(costh_j)

    R = random.random()

    if dsj/W_max1 > R:
        uwe1.append(costh_j)
        if len(uwe1) == Nevents:
            break
        else:
            continue

    else:
        continue

uwe2 = []


for n in range(N):

    costh_j = -1 + random.random()*delta

    dsj = diff_sigma2(costh_j)

    R = random.random()

    if dsj/W_max2 > R:
        uwe2.append(costh_j)
        if len(uwe2) == Nevents:
            break
        else:
            continue

    else:
        continue






#histogram



n_bins = 50









#############################################
#############################################
#############################################

# p = [ E, x , y , z ]
#print("\n")
#print("Initial four-momenta of e+ and e-")
#print("\n")
#print("p1 = ", [E/2., 0.0, 0.0, E/2.])
#print("p2 = ", [E/2., 0.0, 0.0, -E/2.])
#print("\n")
#print("Four-momenta of mu+ and mu-")
#print("\n")
#for n in range(Nevents):

    #phi_n = 2.*np.pi*random.random()
    #sinth_n = np.sqrt(1.-uwe[n]**2.)
    #print("p3 =", [E/2., E/2.*sinth_n*np.cos(phi_n), E/2.*sinth_n*np.sin(phi_n), E/2.*uwe[n]])
    #print("p4 =", [E/2., -E/2.*sinth_n*np.cos(phi_n), -E/2.*sinth_n*np.sin(phi_n), -E/2.*uwe[n]])
    #print("\n")


def norm_cross_section2(costh):

    return(3/8*(1+costh**2 + A12*costh/A02))


def norm_cross_section1(costh):

    return(3/8*(1+costh**2 + A11*costh/A01))





xx = np.linspace(-1, 1, 1000)
plt.figure()
plt.plot(xx, norm_cross_section1(xx), color = "red", label  = "Normalised cross section")
plt.plot(xx, norm_cross_section2(xx), color = "blue", label  = "Normalised cross section")
plt.legend()
plt.ylim(0, 1)
plt.xlim(-1,1)
#plt.title("Probability distribution of the $\\theta$ angle after the hard scattering")
plt.xlabel("cos($\\theta$)", fontsize = 20)
#plt.ylabel("Probability")
plt.hist(uwe1, n_bins, [-1,1], normed = 1, color = "blue", alpha = 0.5, edgecolor = "blue", linewidth = 3, facecolor="None", label="u, c, t")
plt.hist(uwe2, n_bins, [-1,1], normed = 1, color = "red", alpha = 0.5, ls='dashed', lw=3, facecolor="None", label="d, s, b")
plt.legend(fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()
