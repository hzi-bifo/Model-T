import unittest
import mutual_info
import pandas as ps
import numpy as np
import math

class InequalityTest(unittest.TestCase):

    def setUp(self):
        self.mi = mutual_info.MutualInfo()
        self.x = np.array([0, 1, 0, 1])
        self.y = np.array([0, 1, 1, 0])
        self.z = np.array([0, 1, 1, 1])

    def testMarginalProb(self):
        self.assertTrue(np.array_equal(self.mi.prob(self.x), np.array([0.5, 0.5])))

    def testJoint2D(self):
        self.assertTrue(np.array_equal(self.mi.joint_prob(self.x, self.y), np.array([[0.25,0.25],[0.25,0.25]])))

    def testJoint3D(self):
        self.assertTrue(np.array_equal(self.mi.joint_prob_3d(self.x, self.y, self.z), np.array([[[0.25, 0],[0, 0.25]], [[0, 0.25],[0, 0.25]]])))

    def testMI(self):
        self.assertEqual(self.mi.MI(self.x, self.y), 0)
    
    def testCMI(self):
        self.assertEqual(self.mi.CMI(self.x, self.y, self.z), 0.25 * math.log(3.0/4.0) + 0.5 * math.log(3.0/2.0))

    def testEntropy(self):
        self.assertEqual(self.mi.entropy(self.y, self.z), -0.5 * math.log(2.0/3.0))

if __name__ == '__main__':
    unittest.main()
