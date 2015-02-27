import numpy as np
import math
import sys

class MutualInfo:
    def joint_prob_3d(self, x, y, z):
        m = np.zeros(shape = [len(np.unique(i)) for i in [x, y, z]])
        for i, j, k in [(i,j,k) for i in range(m.shape[0]) for j in range(m.shape[1]) for k in range(m.shape[2])]:
            m[i,j,k] = ((x == i) & (y == j)  & (z == k)).sum() / float(x.shape[0])
        return m

    def joint_prob(self, x, y):
        m = np.zeros(shape = [len(np.unique(i)) for i in [x, y]])
        for i, j in [(i,j) for i in range(m.shape[0]) for j in range(m.shape[1])]:
            m[i,j] = ((x == i) & (y == j)).sum() / float(x.shape[0])
        return m

    def prob(self, x):
        return np.array([(x==i).sum()/float(x.shape[0]) for i in range(len(np.unique(x)))])

    def cond_prob(self, x, y):
        m = np.zeros(shape = [len(np.unique(i)) for i in [x, y]])
        for i, j in [(i,j) for i in range(m.shape[0]) for j in range(m.shape[1])]:
            if (x == i).sum() == 0:
                m[i, j] = 0
            m[i, j] = ((x == i) & (y == j)).sum() / float((x == i).sum())
        return m
    
    def MI(self, x, y):
        pxy = self.joint_prob(x, y)
        px = self.prob(x)
        py = self.prob(y)
        MI_sum = 0
        for i, j in [(i,j) for i in range(len(np.unique(x))) for j in range(len(np.unique(y)))]:
            if pxy[i,j] == 0:
                continue
            if px[i] == 0:
                continue
            if py[j] == 0:
                continue
            MI_sum += pxy[i, j] * math.log(pxy[i, j] / float(px[i] * py[j]))
        return MI_sum 

    def CMI(self, x, y, z):
        mi_sum = 0
        pxyz = self.joint_prob_3d(x, y, z)
        pxz = self.joint_prob(x, z)
        pyz = self.joint_prob(y, z)
        pz = self.prob(z)
        for i, j, k in [(i,j,k) for i in range(len(np.unique(x))) for j in range(len(np.unique(y))) for k in range(len(np.unique(z)))]:
            if pxyz[i,j, k] == 0:
                continue
            if pxz[i,k] == 0:
                continue
            if pyz[j,k] == 0:
                continue
            if pz[k] == 0:
                continue
            #print pxyz[i, j, k]  , pz[k], pxz[i,k], pyz[j, k]
            mi_sum += pxyz[i, j, k] * math.log((pxyz[i, j, k] * pz[k])/(pxz[i,k]*pyz[j, k]))
        return mi_sum

    def entropy(self, x, y):
        pxy = self.joint_prob(x, y)
        cond_pyx = self.cond_prob(y, x) 
        entropy_sum = 0
        for i, j in [(i,j) for i in range(len(np.unique(x))) for j in range(len(np.unique(y)))]:
            if pxy[i, j] == 0:
                continue
            if cond_pyx[j, i] == 0:
                continue
            #print pxy[i, j]
            #print cond_pxy[i, j]
            entropy_sum += -pxy[i, j] * math.log(cond_pyx[j, i])
        return entropy_sum


    def CWMI(self, x, y, z):
        entropy = self.entropy(y, z)
        if entropy == 0:
            return 0
        return self.CMI(x, y, z) / self.entropy(y, z)

if __name__ == "__main__":
    a = np.array([1, 1, 0])
    b = np.array([0, 1, 0])
    c = np.array([0, 2, 1])
    mi = MutualInfo()
    print "cond prob", mi.cond_prob(a, b)
    print "joint prob", mi.joint_prob(a, b)
    print "joint prob 3d", mi.joint_prob_3d(a, b, c)
    print "single prob", mi.prob(c)
    print "MI", mi.MI(a, b)
    print "CMI", mi.CMI(a, b, c)
    print "entropy", mi.entropy(a, b)
    print "CWMI", mi.CWMI(a, b, c)
