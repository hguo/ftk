import unittest
import pyftk

class MainTest(unittest.TestCase):
    def test_critical_point_tracking_moving_extremum_2d(self):
        data = pyftk.synthesizers.moving_extremum(11, 13, 20, 5, 5, 0.1, 0.2)
        result = pyftk.trackers.track_critical_points_2d_scalar(data)
        self.assertEqual(len(result), 1)

    def test_critical_point_tracking_spiral_woven(self):
        data = pyftk.synthesizers.spiral_woven(10, 10, 20)
        result = pyftk.trackers.track_critical_points_2d_scalar(data)
        self.assertEqual(len(result), 20)

if __name__ == '__main__':
    unittest.main()
