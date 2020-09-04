import unittest
import pyftk
import numpy as np

class MainTest(unittest.TestCase):
    def test_det4(self):
        m = np.random.rand(4, 4)
        self.assertAlmostEqual(pyftk.numeric.det4(m), np.linalg.det(m), 3)
    
    def test_det3(self):
        m = np.random.rand(3, 3)
        self.assertAlmostEqual(pyftk.numeric.det3(m), np.linalg.det(m), 3)

if __name__ == '__main__':
    unittest.main()
