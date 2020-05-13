import unittest
import numpy as np
import matplotlib.pyplot as plt
from pynamix import color, data, exposure, io, measure, plotting


class TestMeasure(unittest.TestCase):
    def testHanningWindow(self):
        """Test case A. note that all test method names must begin with 'test.'"""
        w = measure.hanning_window()
        plt.imshow(w)
        plt.colorbar()
        plt.show()

    def testAngularBinning(self):
        """Test case A. note that all test method names must begin with 'test.'"""
        n_maskQ = measure.angular_binning(N=1000)
        plt.subplot(221)
        plt.imshow(n_maskQ[:, :, 0, 0])
        plt.colorbar()
        plt.subplot(222)
        plt.imshow(n_maskQ[:, :, 0, 1])
        plt.colorbar()
        plt.subplot(223)
        plt.imshow(n_maskQ[:, :, 1, 0])
        plt.colorbar()
        plt.subplot(224)
        plt.imshow(n_maskQ[:, :, 1, 1])
        plt.colorbar()
        plt.show()

    def testMainDirection(self):
        m = np.array([[1, 0], [1, 0]]) / np.sqrt(2)
        angle, dzeta = measure.main_direction(m)
        print(np.degrees(angle), dzeta)
        m = np.array([[1, 0], [0, 1]]) / np.sqrt(2)
        angle, dzeta = measure.main_direction(m)
        print(np.degrees(angle), dzeta)
        m = np.array([[0, 1], [1, 0]]) / np.sqrt(2)
        angle, dzeta = measure.main_direction(m)
        print(np.degrees(angle), dzeta)

    def testRadialGrid(self):
        r_grid = measure.radial_grid()
        print(r_grid)

    def testRadialFFT(self):
        data, logfile = io.pendulum()
        data, logfile = exposure.apply_ROI(data, logfile, left=600)
        print(np.amax(data))
        data = exposure.clamp(data, 10000, 50000)
        logfile["length"] = {}
        logfile["length"]["height"] = 240.0  # mm
        plt.imshow(data[1000])
        plt.show()
        print(len(logfile["frames"]))
        plt.imshow(data[10])
        plt.show()
        wavelength, radialspec = measure.radial_FFT(data, logfile, tmin=10, tmax=12)
        plt.semilogx(wavelength, radialspec)
        plt.show()

    def testAverageSizeMap():
        data, logfile = io.pendulum()
        data, logfile = exposure.apply_ROI(data, logfile, left=600)
        print(np.amax(data))
        data = exposure.clamp(data, 10000, 50000)
        logfile["length"] = {}
        logfile["length"]["height"] = 240.0  # mm
        x, y, size = measure.average_size_map(data, logfile, tmin=10, tmax=11)

        nt, nx, ny = data.shape
        X, Y = np.meshgrid(range(nx), range(ny), indexing="ij")
        plt.subplot(121)
        plt.pcolormesh(X, Y, data[10])
        plt.gca().invert_yaxis()
        plt.subplot(122)
        plt.pcolormesh(x, y, size[0])
        plt.gca().invert_yaxis()
        plt.colorbar()
        plt.show()

    def textOrientationMap():
        virino = color.virino()
        data, logfile = io.load_seq("/Volumes/LTS/DynamiX/FRC/Marta/MartaTest1-D0.log")
        orient, dzeta = measure.orientation_map(data, tmin=1000, tmax=1002)
        plt.subplot(131)
        plt.imshow(data[1000])
        plt.subplot(132)
        plt.pcolormesh(orient[0], cmap=virino)
        plt.colorbar()
        plt.subplot(133)
        plt.pcolormesh(size[:, :, 0])
        plt.colorbar()
        plt.show()


class OtherTestCase(unittest.TestCase):
    def setUp(self):
        blah_blah_blah()

    def tearDown(self):
        blah_blah_blah()

    def testBlah(self):
        assert self.blahblah == "blah", "blah isn't blahing blahing correctly"


if __name__ == "__main__":
    unittest.main()  # run all tests
