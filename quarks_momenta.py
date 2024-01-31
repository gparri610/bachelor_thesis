import numpy as np
import random as r
import math as m
from vector import Vec4
import matplotlib.pyplot as plt

class Quarks_momenta_generator:

    def __init__(self,E):
        self.E = E
        self.xsec_u, self.W_max_u = self.Cross_Section(2)
        self.xsec_d, self.W_max_d = self.Cross_Section(1)

    def constants(self, c):
        stw = 0.222246
        if c == 1 or c == 3 or c == 5:
            Qf = -1. / 3.
            Vf = -1. / 2. - 2. / 3. * stw
            Af = -1./2.
        else:
            Qf = 2. / 3.
            Vf = 1. / 2. - 4. / 3. * stw
            Af = 1. / 2.

        Vf_m = -1. / 2. + 2. * stw
        Af_m = -1. / 2.

        s = self.E ** 2
        Mz = 91.1876
        Tz = 2.4414
        alpha = 1. / (132.507)
        Gf = 1.16639 * 10 ** (-5)

        kappa = m.sqrt(2.) * Gf * Mz ** 2. / (4. * m.pi * alpha)
        chi1 = kappa * s * (s - Mz ** 2.) / ((s - Mz ** 2.) ** 2. + Tz ** 2. * Mz ** 2.)
        chi2 = kappa ** 2. * s ** 2. / ((s - Mz ** 2.) ** 2. + Tz ** 2. * Mz ** 2.)

        A0 = Qf ** 2 - 2 * Qf * Vf * Vf_m * chi1 + (Af_m ** 2 + Vf_m ** 2) * (Af ** 2 + Vf ** 2) * chi2
        A1 = -4 * Qf * Af * Af_m * chi1 + 8 * Af_m * Vf_m * Af * Vf * chi2

        return A0, A1, alpha, s


    def diff_cross_section(self, costh, c):

        A0 = self.constants(c)[0]
        A1 = self.constants(c)[1]
        alpha = self.constants(c)[2]
        s = self.constants(c)[3]

        return (2. * m.pi) * alpha ** 2 * (A0 * (1. + costh ** 2) + A1 * costh) / (4. * s)


    def Cross_Section(self, c):

        pb_convert = 3.894 * 10 ** 8
        delta = 2
        sum_W = 0
        sum_W_sq = 0
        array = []
        N = 10000



        for i in range(N):

            costh_i = -1. + r.random() * delta

            W_i = delta * self.diff_cross_section(costh_i, c)
            array.append(W_i)


            sum_W = sum_W + W_i

            sum_W_sq = sum_W_sq + W_i ** 2


        sigma = sum_W / N

        return(sigma*pb_convert), max(array)


    def event_generator_p1(self, c):


        if c==1 or c==3 or c==5:
            W_max = self.W_max_d
        else:
            W_max = self.W_max_u

        delta = 2
        uwe = []


        for j in range(100000):

            costh_j = -1 + r.random() * delta

            dsj = self.diff_cross_section(costh_j, c)

            if dsj / W_max > r.random():
                uwe.append(costh_j)
                if len(uwe) == 1:
                    break
                else:
                    continue

            else:
                continue


        phi = 2. * m.pi * r.random()
        sinth = np.sqrt(1. - uwe[0]**2.)
        p1 = Vec4(self.E / 2., self.E / 2. * sinth * m.cos(phi), self.E / 2. * sinth * m.sin(phi), self.E / 2. * uwe[0])

        return p1

    def event_generator_p2(self, p1):

        p2 = Vec4(p1[0], -p1[1], -p1[2], -p1[3])

        return p2

E=91.18
a = Quarks_momenta_generator(E)
d,_ =a.Cross_Section(1)
u,_ =a.Cross_Section(2)
print(d)
print(u)
print("P(d)=")
print(d/(2.*u+3.*d))
print("P(u)=")
print(u/(2.*u+3.*d))




