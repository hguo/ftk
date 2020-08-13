import unittest
import pyftk

class MainTest(unittest.TestCase):
    def test_critical_point_tracking_moving_extremum_2d(self):
        data = pyftk.synth.moving_extremum(11, 13, 20, 5, 5, 0.1, 0.2)
        result = pyftk.track_critical_points_2d_scalar(data)
        self.assertEqual(len(result), 1)

if __name__ == '__main__':
    unittest.main()
