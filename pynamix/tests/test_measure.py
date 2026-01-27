import unittest
import numpy as np
from pynamix import measure


class TestMeasureModule(unittest.TestCase):
    """Test cases for the measure module"""

    def test_main_direction_horizontal(self):
        """Test main_direction with horizontal orientation"""
        # Perfectly horizontal tensor
        tensor = np.array([[1, 0], [0, 0]], dtype=float)
        angle, dzeta = measure.main_direction(tensor)
        
        # Angle should be 0 or pi (horizontal)
        self.assertTrue(abs(angle) < 0.01 or abs(angle - np.pi) < 0.01)
        
        # dzeta is magnitude
        self.assertGreater(dzeta, 0)

    def test_main_direction_vertical(self):
        """Test main_direction with vertical orientation"""
        # Perfectly vertical tensor
        tensor = np.array([[0, 0], [0, 1]], dtype=float)
        angle, dzeta = measure.main_direction(tensor)
        
        # Angle should be pi/2 (vertical)
        self.assertAlmostEqual(angle, np.pi / 2, places=5)
        
        # dzeta is magnitude
        self.assertGreater(dzeta, 0)

    def test_main_direction_diagonal(self):
        """Test main_direction with diagonal orientation"""
        # 45-degree diagonal tensor
        tensor = np.array([[1, 1], [1, 1]], dtype=float) / np.sqrt(2)
        angle, dzeta = measure.main_direction(tensor)
        
        # Angle should be pi/4 (45 degrees)
        self.assertAlmostEqual(angle, np.pi / 4, places=5)

    def test_main_direction_range(self):
        """Test that main_direction returns angle in [0, pi]"""
        # Test with various tensors
        for _ in range(10):
            tensor = np.random.rand(2, 2)
            angle, dzeta = measure.main_direction(tensor)
            
            # Angle should be in [0, pi]
            self.assertGreaterEqual(angle, 0)
            self.assertLessEqual(angle, np.pi)


class TestHanningWindow(unittest.TestCase):
    """Test Hanning window generation"""

    def test_hanning_window_default(self):
        """Test hanning_window with default patch size"""
        w = measure.hanning_window()
        
        # Default patchw is 32, so window should be 64x64
        self.assertEqual(w.shape, (64, 64))
        
        # Values should be between 0 and 1
        self.assertGreaterEqual(np.min(w), 0)
        self.assertLessEqual(np.max(w), 1)

    def test_hanning_window_custom_size(self):
        """Test hanning_window with custom patch size"""
        patchw = 16
        w = measure.hanning_window(patchw)
        
        # Window should be 2*patchw x 2*patchw
        self.assertEqual(w.shape, (32, 32))

    def test_hanning_window_center_maximum(self):
        """Test that hanning window has maximum near center"""
        w = measure.hanning_window(32)
        
        # Center should have higher values than edges
        center_val = w[32, 32]
        edge_val = w[0, 0]
        
        self.assertGreater(center_val, edge_val)

    def test_hanning_window_radial_properties(self):
        """Test that hanning window has correct radial properties"""
        w = measure.hanning_window(32)
        
        # Check that center has high value (near 1)
        self.assertGreater(w[32, 32], 0.99)
        
        # Check that values decrease with distance from center
        # Points closer to center should have higher values
        self.assertGreater(w[32, 32], w[25, 32])
        self.assertGreater(w[25, 32], w[20, 32])
        self.assertGreater(w[20, 32], w[10, 32])
        
        # Check that corners (far from center) are zero or very small
        self.assertLess(w[1, 1], 0.01)
        self.assertLess(w[1, 63], 0.01)
        self.assertLess(w[63, 1], 0.01)
        self.assertLess(w[63, 63], 0.01)
        
        # Check that window is zero or very small outside radius (patchw=32)
        # Points at distance > 32 should be zero or nearly zero
        self.assertLess(w[1, 32], 0.01)  # distance ~31, should be small
        self.assertLess(w[63, 32], 0.01)  # distance ~31, should be small

    def test_hanning_window_zero_outside_radius(self):
        """Test that hanning window is zero outside radius"""
        patchw = 32
        w = measure.hanning_window(patchw)
        
        # Corners should be zero (distance > patchw)
        self.assertEqual(w[0, 0], 0)
        self.assertEqual(w[0, -1], 0)
        self.assertEqual(w[-1, 0], 0)
        self.assertEqual(w[-1, -1], 0)


class TestGrid(unittest.TestCase):
    """Test grid generation for patch analysis"""

    def test_grid_basic_no_ROI(self):
        """Test grid without ROI in logfile"""
        data = np.zeros((10, 100, 80))  # nt, nx, ny
        logfile = {"detector": {}}
        xstep, ystep, patchw = 16, 16, 8
        
        gridx, gridy = measure.grid(data, logfile, xstep, ystep, patchw)
        
        # Grid should start at patchw and end at nx-patchw
        self.assertEqual(gridx[0], patchw)
        self.assertLess(gridx[-1], 100 - patchw)
        
        # Check spacing
        if len(gridx) > 1:
            self.assertEqual(gridx[1] - gridx[0], xstep)

    def test_grid_with_ROI(self):
        """Test grid with ROI in logfile"""
        data = np.zeros((10, 100, 80))
        logfile = {
            "detector": {
                "ROI": {
                    "left": 10,
                    "right": 90,
                    "top": 5,
                    "bottom": 75
                }
            }
        }
        xstep, ystep, patchw = 16, 16, 8
        
        gridx, gridy = measure.grid(data, logfile, xstep, ystep, patchw)
        
        # Grid should respect ROI boundaries (within the ROI region)
        # Check that grids are non-empty before accessing elements
        self.assertGreater(len(gridx), 0, "gridx should not be empty")
        self.assertGreater(len(gridy), 0, "gridy should not be empty")
        
        # gridx should start from left + patchw
        self.assertGreaterEqual(gridx[0], logfile["detector"]["ROI"]["left"] + patchw)
        # gridx should end before right - patchw
        self.assertLessEqual(gridx[-1], logfile["detector"]["ROI"]["right"] - patchw)
        
        # gridy should start from top + patchw
        self.assertGreaterEqual(gridy[0], logfile["detector"]["ROI"]["top"] + patchw)
        # gridy should end before bottom - patchw
        self.assertLessEqual(gridy[-1], logfile["detector"]["ROI"]["bottom"] - patchw)

    def test_grid_centered_mode(self):
        """Test grid with centered mode"""
        data = np.zeros((10, 100, 80))
        logfile = {"detector": {}}
        xstep, ystep, patchw = 16, 16, 8
        
        gridx, gridy = measure.grid(data, logfile, xstep, ystep, patchw, mode="center")
        
        # Grid should be non-empty
        self.assertGreater(len(gridx), 0)
        self.assertGreater(len(gridy), 0)
        
        # Grid should be within valid bounds
        self.assertGreaterEqual(gridx[0], patchw)
        self.assertLessEqual(gridx[-1], 100 - patchw)
        
        # Check that grid is reasonably centered
        # The first point should not be exactly at patchw (should have offset)
        self.assertGreater(gridx[0], patchw)

    def test_grid_full_mode(self):
        """Test grid with full mode (no buffer)"""
        data = np.zeros((10, 100, 80))
        logfile = {"detector": {}}
        xstep, ystep, patchw = 16, 16, 8
        
        gridx, gridy = measure.grid(data, logfile, xstep, ystep, patchw, mode="full")
        
        # Grid should start at 0
        self.assertEqual(gridx[0], 0)
        self.assertEqual(gridy[0], 0)
        
        # Grid should be non-empty
        self.assertGreater(len(gridx), 0)
        self.assertGreater(len(gridy), 0)

    def test_grid_invalid_mode(self):
        """Test grid with invalid mode raises error"""
        data = np.zeros((10, 100, 80))
        logfile = {"detector": {}}
        xstep, ystep, patchw = 16, 16, 8
        
        with self.assertRaises(ValueError):
            measure.grid(data, logfile, xstep, ystep, patchw, mode="invalid")

    def test_grid_returns_1d_arrays(self):
        """Test that grid returns 1D arrays"""
        data = np.zeros((10, 100, 80))
        logfile = {"detector": {}}
        xstep, ystep, patchw = 16, 16, 8
        
        gridx, gridy = measure.grid(data, logfile, xstep, ystep, patchw)
        
        self.assertEqual(len(gridx.shape), 1)
        self.assertEqual(len(gridy.shape), 1)


class TestAngularBinning(unittest.TestCase):
    """Test angular binning for Q coefficients"""

    def test_angular_binning_default(self):
        """Test angular_binning with default parameters"""
        # This will take some time, so use smaller N for testing
        n_maskQ = measure.angular_binning(patchw=8, N=100)
        
        # Shape should be [2*patchw, 2*patchw, 2, 2]
        self.assertEqual(n_maskQ.shape, (16, 16, 2, 2))
        
        # Values should be finite
        self.assertTrue(np.all(np.isfinite(n_maskQ)))

    def test_angular_binning_symmetry(self):
        """Test that Q coefficients have expected symmetry"""
        n_maskQ = measure.angular_binning(patchw=8, N=100)
        
        # Q[i,j,0,1] should equal Q[i,j,1,0] (symmetry of tensor)
        diff = np.abs(n_maskQ[:, :, 0, 1] - n_maskQ[:, :, 1, 0])
        # Allow for numerical errors
        self.assertLess(np.max(diff), 0.1)

    def test_angular_binning_values_range(self):
        """Test that Q coefficients are in reasonable range"""
        n_maskQ = measure.angular_binning(patchw=8, N=100)
        
        # Values should be between -1 and 1 for normalized tensor components
        # (after removing NaNs from division by zero)
        finite_vals = n_maskQ[np.isfinite(n_maskQ)]
        self.assertGreaterEqual(np.min(finite_vals), -2)
        self.assertLessEqual(np.max(finite_vals), 2)


class TestRadialGrid(unittest.TestCase):
    """Test radial grid generation"""

    def test_radial_grid_default(self):
        """Test radial_grid with default parameters"""
        # Use smaller parameters for faster testing
        r_grid, nr_pxr = measure.radial_grid(rnb=50, patchw=8, N=100)
        
        # r_grid should be 1D with rnb elements
        self.assertEqual(len(r_grid), 50)
        
        # nr_pxr should be 3D
        self.assertEqual(nr_pxr.shape, (16, 16, 50))
        
        # r_grid should be increasing
        self.assertTrue(np.all(np.diff(r_grid) > 0))

    def test_radial_grid_range(self):
        """Test that radial grid spans expected range"""
        patchw = 8
        r_grid, nr_pxr = measure.radial_grid(rnb=50, patchw=patchw, N=100)
        
        # Grid should start near 0
        self.assertLess(r_grid[0], 1)
        
        # Grid should end around patchw * 1.5
        self.assertGreater(r_grid[-1], patchw * 1.3)
        self.assertLess(r_grid[-1], patchw * 1.7)

    def test_radial_grid_nr_pxr_normalized(self):
        """Test that nr_pxr values are normalized (sum to ~1)"""
        r_grid, nr_pxr = measure.radial_grid(rnb=50, patchw=8, N=100)
        
        # For any pixel, sum over all radii should be ~1
        # Pick a pixel near center
        pixel_sum = np.sum(nr_pxr[8, 8, :])
        
        # Should be close to 1 (normalized probability)
        self.assertGreater(pixel_sum, 0.5)
        self.assertLess(pixel_sum, 1.5)


if __name__ == "__main__":
    unittest.main()
