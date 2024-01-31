import math as m
import random as r
from vector import Vec4
from particle import Particle, CheckEvent
from qcd import AlphaS, NC, TR, CA, CF
from quarks_momenta import Quarks_momenta_generator
import pyjet
import numpy as np
from pyjet.utils import ptepm2ep
import matplotlib.pyplot as plt
from scipy import optimize

class Shower:



    def __init__(self, alpha, t0):
        self.t0 = t0
        self.alpha = alpha
        self.alphamax = alpha(t0)
        self.kernels = [Pqq([fl, fl, 21]) for fl in [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]]
        self.kernels += [Pgq([21, fl, -fl]) for fl in [1, 2, 3, 4, 5]]
        self.kernels += [Pgg([21, 21, 21])]


    def MakeKinematics(self, z, y, phi, pa_t, ps_t):

        CM = pa_t + ps_t
        pa_tCM = Vec4.Boost(CM, pa_t)
        p = pa_tCM/pa_tCM.P()
        if p[1] == 0.0:
            e1 = Vec4(0., 1., 0., 0.)
        else:
            e1 = Vec4(0., -p[3]/p[1], 0., 1.)
            e1 = e1/e1.P()
        e2 = p.Cross(e1)
        ktmod= np.sqrt( 2.*(pa_t*ps_t)*z*(1-z)*y)
        ktCM = ktmod*(m.cos(phi)*e1 + m.sin(phi)*e2)
        kt = Vec4.BoostBack(CM, ktCM)
        pb = z*pa_t + (1-z)*y*ps_t + kt
        pc = (1-z)*pa_t + y*z*ps_t - kt
        ps = (1-y)*ps_t

        new_momenta = [pb, pc, ps]

        return new_momenta

    def MakeColors(self, flavs, colbc, cols):
        self.c += 1
        if flavs[0] != 21:
            if flavs[0] > 0:
                return [[self.c, 0], [colbc[0], self.c]]
            else:
                return [[0, self.c], [self.c, colbc[1]]]
        else:
            if flavs[1] == 21:
                if colbc[0] == cols[1]:
                    if colbc[1] == cols[0] and r.random() > 0.5:
                        return [[colbc[0], self.c], [self.c, colbc[1]]]
                    return [[self.c, colbc[1]], [colbc[0], self.c]]
                else:
                    return [[colbc[0], self.c], [self.c, colbc[1]]]
            else:
                if flavs[1] > 0:
                    return [[colbc[0], 0], [0, colbc[1]]]
                else:
                    return [[0, colbc[1]], [colbc[0], 0]]


    def GeneratePoint(self, event):

        while self.t > self.t0:

            t = self.t0

            for split in event[2:]:
                for spect in event[2:]:
                    if spect == split: continue
                    if not split.ColorConnected(spect): continue

                    for sf in self.kernels:
                        if sf.flavs[0] != split.pid: continue

                        m2 = (split.mom + spect.mom).M2()
                        if m2 < 4.*self.t0: continue
                        zp = 0.5*(1. + m.sqrt(1.-4.*self.t0/m2))
                        zm = 1.-zp

                        g = self.alphamax/(2.*m.pi)*sf.Integral(zm,zp)

                        tt = self.t*m.pow(r.random(), 1./g)

                        if tt>t:
                            t = tt
                            s = [split, spect, sf, m2, zp]

            self.t = t

            if t > self.t0:
                z = s[2].GenerateZ(1.-s[4], s[4])

                y = t/s[3]/z/(1.-z)

                if y < 1.:

                    w = (1.-y)*self.alpha(t)/self.alphamax
                    w*= s[2].Value(z,y)/s[2].Estimate(z)

                    if w > r.random():

                        phi = 2*m.pi*r.random()
                        moms = self.MakeKinematics(z,y,phi, s[0].mom, s[1].mom)
                        cols = self.MakeColors(s[2].flavs, s[0].col, s[1].col)
                        event.append(Particle(s[2].flavs[2], moms[1], cols[1]))
                        s[0].Set(s[2].flavs[1],moms[0],cols[0])
                        s[1].mom = moms[2]

                        return


    def Run(self,event):
        self.c = 1
        self.t = (event[0].mom+event[1].mom).M2()
        while self.t > self.t0:
            self.GeneratePoint(event)





class Kernel:

    def __init__(self,flavs):
        self.flavs = flavs



class Pqq (Kernel):

    def Value(self,z,y):
        return CF*(2./(1.-z*(1.-y))-(1.+z))

    def Estimate(self,z):
        return CF*2./(1.-z)

    def Integral(self,zm,zp):

        return CF*2.*m.log((1.-zm)/(1.-zp))

    def GenerateZ(self,zm,zp):
        return 1.+(zp-1.)*m.pow((1.-zm)/(1.-zp),r.random())




class Pgg (Kernel):

    def Value(self, z, y):
        return CA/2.*(2./(1.-z*(1.-y))-2.+z*(1.-z))

    def Estimate(self, z):
        return CA*1./(1-z)

    def Integral(self, zm, zp):
        return CA*m.log((1.-zm)/(1.-zp))

    def GenerateZ(self, zm, zp):
        return 1. + (zp-1.)*m.pow((1.-zm)/(1.-zp), r.random())




class Pgq (Kernel):

    def Value(self, z, y):
        return TR/2.*(1-2*z*(1-z))

    def Estimate(self, z):
        return TR/2.

    def Integral(self, zm, zp):
        return TR/2.*(zp-zm)

    def GenerateZ(self, zm,zp):
        rr = r.random()
        return (zp-zm)*rr + zm

class generate_hist():

    def __init__(self, array):
        self.a = array

    def hist_number_jets(self):

        plt.figure()

        bins = np.arange(1, max(self.a) + 1.5) - 0.5

        plt.hist(self.a, bins, color = "blue", alpha = 0.5 ,edgecolor="black", rwidth=1)
        plt.xticks(np.arange(2, max(self.a) + 1, 1), fontsize = 16)
        plt.xlim(1.5, 14.5)
        plt.yticks(fontsize=16)
        plt.xlabel("Jet multiplicity", fontsize= 20)
        plt.ylabel("Events", fontsize= 20)

        #plt.title("Jet multiplicity at $E_{CM}$ = $M_Z$")

        plt.show()

    def hist_Z_mass(self):

        bins = np.linspace(20, 92, 100)

        plt.figure()


        plt.hist(self.a, bins,  color = "blue", alpha = 0.5 )
        plt.xlabel("$M_z$   [GeV]", fontsize= 20)
        plt.ylabel("Events", fontsize= 20)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        #plt.title("Mass of the boson Z")
        plt.show()


    def hist_Z_p_T(self):

        bins = np.linspace(0, 25, 150)

        plt.figure()


        plt.hist(self.a, bins,  color = "blue", alpha = 0.5)
        #plt.title("Recoil of the jets")
        plt.xlabel("$|p_T^Z|$   [GeV]", fontsize= 20)
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize=16)
        plt.ylabel("Events", fontsize= 20)
        plt.show()

    #def fit_function(self, x, A, B):

     #   return A*m.log(B*x)

    def hist_multiplicity_E(self, N):



        MULTI = []

        for i in range(len(self.a)):
            p_e1 = Vec4(self.a[i] /2., self.a[i] / 2., 0., self.a[i] / 2.)
            p_e2 = Vec4(self.a[i] /2., -self.a[i] / 2., 0., -self.a[i] / 2.)

            quarks = Quarks_momenta_generator(self.a[i])

            NJETS = []
            for k in range(N):
                c = r.randint(1, 5)
                p_q1 = quarks.event_generator_p1(c)
                p_q2 = quarks.event_generator_p2(p_q1)
                event = [Particle(-11, p_e1), Particle(11, p_e2), Particle(c, p_q1, col=[1, 0]),
                         Particle(-c, p_q2, col=[0, 1])]
                myshower.Run(event)
                moms = []
                for x in event:
                    moms.append((x.mom.E, x.mom.px, x.mom.py, x.mom.pz))
                vectors = np.array(moms[2:], dtype=([('E', 'f8'), ('px', 'f8'), ('py', 'f8'), ('pz', 'f8')]))
                sequence = pyjet.cluster(vectors, algo = "ee_genkt", R=0.2, p=-1, ep=True)
                jets = sequence.inclusive_jets()
                njets = len(jets)
                NJETS.append(njets)


            MULTI.append(sum(NJETS)/len(NJETS))

        #params, params_covariance = optimize.curve_fit(self.fit_function, self.a, MULTI, p0=[10., 1.])

        plt.figure()
        #x = np.array([10,20,30,40,50,60,70])
        plt.plot(self.a, MULTI, color = "blue", alpha = 0.5, linewidth = 5)
        plt.scatter(self.a, MULTI, color = "b", linewidth = 5)
        plt.ylabel("Average jet multiplicity", fontsize= 20)
        #my_xticks = ['10', '20', '50', '100', '200','500','1000']
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize=16)
        #plt.ylim(3.35, 3.85)
        plt.xlabel("$\sqrt{s}$   [GeV]", fontsize= 20)
        plt.show()



    def hist_multiplicity_r(self, N):

        MULTI = []

        E = 91.1876

        p_e1 = Vec4(E / 2., E / 2., 0., E / 2.)
        p_e2 = Vec4(E / 2., -E / 2., 0., -E / 2.)

        for i in range(len(self.a)):
            NJETS = []
            for j in range(N):
                c = r.randint(1, 5)
                p_q1 = quarks.event_generator_p1(c)
                p_q2 = quarks.event_generator_p2(p_q1)
                event = [Particle(-11, p_e1), Particle(11, p_e2), Particle(c, p_q1, col=[1, 0]),
                         Particle(-c, p_q2, col=[0, 1])]
                myshower.Run(event)
                moms = []
                for x in event:
                    moms.append((x.mom.E, x.mom.px, x.mom.py, x.mom.pz))
                vectors = np.array(moms[2:], dtype=([('E', 'f8'), ('px', 'f8'), ('py', 'f8'), ('pz', 'f8')]))
                sequence = pyjet.cluster(vectors, algo="ee_genkt", R=self.a[i], p=-1, ep=True)
                jets = sequence.inclusive_jets()
                njets = len(jets)
                NJETS.append(njets)

            MULTI.append(sum(NJETS) / len(NJETS))

        plt.figure()
        plt.bar(self.a, MULTI, width = 0.1, color = "blue", alpha = 0.5 , edgecolor=['black']*len(self.a))
        #plt.title("Jet multiplicity vs Radius ")
        plt.ylabel("Average jet multiplicity", fontsize= 20)
        plt.xticks(self.a, fontsize = 16)
        plt.yticks(fontsize=16)
        plt.xlabel("R", fontsize= 20)
        plt.show()


N = 100000
r.seed(1408)


E = 91.1876

p_e1 = Vec4(E/2., E/2., 0., E/2.)
p_e2 = Vec4(E/2., -E/2., 0., -E/2.)

quarks = Quarks_momenta_generator(E)
myshower = Shower(AlphaS(91.1876, 0.118), 1.)  # myshower = Shower(alpha_(mZ),t0)


NJETS = []
NJETS00 = []
PROVA = []
NJETS03 = []
NJETS002 = []
NJETS032 = []
MOMS = []
JETS = []
for n in range(N):
    c = r.randint(1, 5)
    p_q1 = quarks.event_generator_p1(c)
    p_q2 = quarks.event_generator_p2(p_q1)
    event = [Particle(-11, p_e1), Particle(11, p_e2), Particle(c, p_q1, col=[1,0]), Particle(-c, p_q2, col=[0,1])]
    myshower.Run(event)
    moms=[]
    mm = []
    for x in event:
        #print(x)

        moms.append((x.mom.E,x.mom.px,x.mom.py,x.mom.pz))
    vectors = np.array(moms[2:], dtype=([('E','f8'),('px','f8'),('py','f8'),('pz','f8')]))
    #print("")

    MOMS.append(moms[2:])

    #print(vectors)
    sequence = pyjet.cluster(vectors, algo="ee_genkt", R=0.2, p=-1, ep=True)
    jets = sequence.inclusive_jets()
    njets = len(jets)
    NJETS.append(njets)
    jet1 = ptepm2ep(np.array([(jets[0].pt,jets[0].eta,jets[0].phi,jets[0].mass)],
                              dtype=([('pt','f8'),('eta','f8'),('phi','f8'),('mass','f8')])))
    jet2 = ptepm2ep(np.array([(jets[1].pt, jets[1].eta, jets[1].phi, jets[1].mass)],
                             dtype=([('pt', 'f8'), ('eta', 'f8'), ('phi', 'f8'), ('mass', 'f8')])))


    jet12 = []
    for i in range(4):
       jet12.append(jet1[0][i] + jet2[0][i])
    PROVA.append(jet12)
    for i in range(4):
       jet12.append(moms[2][i] + moms[3][i])
    if(njets>2):
      mass = np.sqrt(jet12[0]**2 - jet12[1]**2 - jet12[2]**2 - jet12[3]**2)
      pt = np.sqrt(jet12[1]**2+jet12[2]**2)
      NJETS00.append(mass)
      NJETS03.append(pt)


    jet122 = []     ##### label "2" means the 2 jets event are also considered
    for i in range(4):
        jet122.append(jet1[0][i] + jet2[0][i])
    for i in range(4):
        jet122.append(moms[2][i] + moms[3][i])
    mass2 = np.sqrt(jet122[0] ** 2 - jet122[1] ** 2 - jet122[2] ** 2 - jet122[3] ** 2)
    pt2 = np.sqrt(jet122[1] ** 2 + jet122[2] ** 2)
    NJETS002.append(mass2)
    NJETS032.append(pt2)

    #for n in range(njets):
        #jets_n = ptepm2ep(np.array([(jets[n].pt,jets[n].eta,jets[n].phi,jets[n].mass)],
                             # dtype=([('pt','f8'),('eta','f8'),('phi','f8'),('mass','f8')])))
        #JETS.append(jets_n)


print(len(NJETS00)/len(PROVA) * 100)
#print(len(MOMS2))
#print(len(MOMS))

##### 1 #####
#x1 = generate_hist(NJETS)
#x1.hist_number_jets()
#############

##### 2 #####
#x2 = generate_hist(NJETS00)
#x2.hist_Z_mass()
#############

##### 3 #####
x3 = generate_hist(NJETS03)
x3.hist_Z_p_T()
#############

#############
#x22 = generate_hist(NJETS002)
#x22.hist_Z_mass()
#############

#############
#x32 = generate_hist(NJETS032)
#x32.hist_Z_p_T()
#############

##### 4 #####
#rr = np.linspace(0.1, 1.7, 17)
#x4 = generate_hist(rr)
#x4.hist_multiplicity_r(N)
#############
#print(rr)

##### 6 ######
def T(x, moms):


    phi = x[0]
    theta = x[1]
    nx = m.sin(theta)*m.cos(phi)
    ny = m.sin(theta)*m.sin(phi)
    nz = m.cos(theta)

    s = 0
    d = 0
    for j in range(len(moms)):
        s += abs(nx*moms[j][1] + ny*moms[j][2] + nz*moms[j][3])

        d += abs(m.sqrt(moms[j][1]**2 + moms[j][2]**2 + moms[j][3]**2))


    return 1.-s/d


THRUST = []

for i in range(len(MOMS)):

    moms_i = MOMS[i]


    if len(moms_i)==2:
        continue

    x0 = np.array([m.pi, m.pi/2.])

    bnds=((0,2*m.pi),(0,m.pi))
    maxx = optimize.minimize(T, x0, args= moms_i, bounds=bnds)

    if 1-T([maxx.x[0], maxx.x[1]], moms_i)<1./2.:
        x0 = np.array([m.pi*3./2., m.pi*2./5.])
        maxx = optimize.minimize(T, x0, args= moms_i, bounds=bnds)
        THRUST.append(1-T([maxx.x[0], maxx.x[1]], moms_i))

    else:
        thrust = 1-T([maxx.x[0], maxx.x[1]], moms_i)
        THRUST.append(thrust)

   # if 1-T([maxx.x[0], maxx.x[1]], moms_i) > 0.999:
        #print(moms_i)
        #print(maxx.fun)
      #  phi = maxx.x[0]
      #  theta = maxx.x[1]
      #  nx = m.sin(theta) * m.cos(phi)
       # ny = m.sin(theta) * m.sin(phi)
        #nz = m.cos(theta)

#print(len(THRUST))

nn = np.linspace(.5, 1., 120)

#plt.figure()
plt.hist(THRUST, nn,  color = "blue", alpha = 0.5)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel("T", fontsize=20)
plt.ylabel("Events", fontsize = 20)
#plt.show()

#############
##### 5 #####
#EE = np.array([10., 20., 50., 100., 200., 500., 1000.])
#x5 = generate_hist(EE)
#x5.hist_multiplicity_E(N)
#############

