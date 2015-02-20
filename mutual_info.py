import numpy as np
import math

class MutualInfo:
    def joint_prob(self, x, y, z = None):
        m = np.zeros(shape = (2,2,2))
        z_is_none = False 
        if z is None:
            z = np.zeros(shape = x.shape)
            z_is_none = True
        for i, j, k in [(i,j,k) for i in range(2) for j in range(2) for k in range(2)]:
            m[i,j,k] = ((x == i) & (y == j)  & (z == k)).sum() / float(x.shape[0])
        if z_is_none:
            return m[:, :, 0]
        else:
            return m

    def prob(self, x):
        return np.array([x[x==0].sum()/float(x.shape[0]), x[x==1].sum()/float(x.shape[0])])

    def cond_prob(self, x, y):
        m = np.zeros(shape = (2,2))
        for i, j in [(i,j) for i in range(2) for j in range(2)]:
            if (x == i).sum() == 0:
                m[i, j] = 0
            m[i, j] = ((x == i) & (y == j)).sum() / float((x == i).sum())
        return m
    
    def MI(self, x, y):
        pxy = self.joint_prob(x, y)
        px = self.prob(x)
        py = self.prob(y)
        MI_sum = 0
        for i, j in [(i,j) for i in range(2) for j in range(2)]:
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
        pxyz = self.joint_prob(x, y, z)
        pxz = self.joint_prob(x, z)
        pyz = self.joint_prob(y, z)
        pz = self.prob(z)
        for i, j, k in [(i,j,k) for i in range(2) for j in range(2) for k in range(2)]:
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
        cond_pxy = self.cond_prob(x, y) 
        entropy_sum = 0
        for i, j in [(i,j) for i in range(2) for j in range(2)]:
            if pxy[i, j] == 0:
                continue
            if cond_pxy[i, j] == 0:
                continue
            print pxy[i, j]
            print cond_pxy[i, j]
            entropy_sum += -pxy[i, j] * math.log(cond_pxy[i, j])
        return entropy_sum


    def CWMI(self, x, y, z):
        entropy = self.entropy(y, z)
        if entropy == 0:
            return 0
        return self.CMI(x, y, z) / self.entropy(y, z)

if __name__ == "__main__":
    a = np.array([1, 1, 0])
    b = np.array([0, 1, 0])
    c = np.array([1, 1, 1])
    mi = MutualInfo()
    print "cond prob", mi.cond_prob(a, b)
    print "joint prob", mi.joint_prob(a, b)
    print "MI", mi.MI(a, b)
    print "CMI", mi.CMI(a, b, c)
    print "entropy", mi.entropy(a, b)
    print "CWMI", mi.CWMI(a, b, c)
