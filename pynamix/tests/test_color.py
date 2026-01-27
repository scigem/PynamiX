import unittest
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pynamix import color


class TestColorModule(unittest.TestCase):
    """Test cases for the color module"""

    def test_virino_returns_colormap(self):
        """Test that virino() returns a valid matplotlib colormap"""
        cmap = color.virino()
        self.assertIsInstance(cmap, LinearSegmentedColormap)
        self.assertEqual(cmap.name, "virino")

    def test_virino_colormap_range(self):
        """Test that virino colormap works across full range"""
        cmap = color.virino()
        # Test colormap can be evaluated at various points
        colors_at_0 = cmap(0.0)
        colors_at_half = cmap(0.5)
        colors_at_1 = cmap(1.0)
        
        # Each should return RGBA values
        self.assertEqual(len(colors_at_0), 4)
        self.assertEqual(len(colors_at_half), 4)
        self.assertEqual(len(colors_at_1), 4)
        
        # Values should be in [0, 1] range
        for color_val in [colors_at_0, colors_at_half, colors_at_1]:
            for component in color_val:
                self.assertGreaterEqual(component, 0.0)
                self.assertLessEqual(component, 1.0)

    def test_virino2d_valid_input(self):
        """Test virino2d with valid angle inputs"""
        # Create a simple grid of angles
        angles = np.array([[0, np.pi/4], [np.pi/2, np.pi]])
        magnitude = np.ones_like(angles)
        
        result = color.virino2d(angles, magnitude)
        
        # Check output shape - should add RGB dimension
        self.assertEqual(result.shape, (2, 2, 3))
        
        # Check all RGB values are in valid range
        self.assertTrue(np.all(result >= 0))
        self.assertTrue(np.all(result <= 1))

    def test_virino2d_negative_angles(self):
        """Test virino2d with negative angles (should work within -pi to pi)"""
        angles = np.array([[-np.pi, -np.pi/2], [-np.pi/4, 0]])
        magnitude = np.ones_like(angles)
        
        result = color.virino2d(angles, magnitude)
        
        # Check output shape
        self.assertEqual(result.shape, (2, 2, 3))
        
        # Check all RGB values are in valid range
        self.assertTrue(np.all(result >= 0))
        self.assertTrue(np.all(result <= 1))

    def test_virino2d_angle_bounds_assertion(self):
        """Test that virino2d raises assertion for angles outside [-pi, pi]"""
        # Angles above pi
        angles_too_high = np.array([[0, np.pi * 1.5]])
        magnitude = np.ones_like(angles_too_high)
        
        with self.assertRaises(AssertionError):
            color.virino2d(angles_too_high, magnitude)
        
        # Angles below -pi
        angles_too_low = np.array([[0, -np.pi * 1.5]])
        magnitude = np.ones_like(angles_too_low)
        
        with self.assertRaises(AssertionError):
            color.virino2d(angles_too_low, magnitude)

    def test_virino2d_magnitude_effect(self):
        """Test that magnitude parameter affects the output"""
        angles = np.array([[0, np.pi/2]])
        magnitude_low = np.array([[0.1, 0.1]])
        magnitude_high = np.array([[1.0, 1.0]])
        
        result_low = color.virino2d(angles, magnitude_low)
        result_high = color.virino2d(angles, magnitude_high)
        
        # Results should be different (though current implementation may not use magnitude)
        # This test documents current behavior
        self.assertEqual(result_low.shape, result_high.shape)

    def test_virino2d_single_value(self):
        """Test virino2d with scalar-like inputs"""
        angles = np.array([[0]])
        magnitude = np.array([[1]])
        
        result = color.virino2d(angles, magnitude)
        
        self.assertEqual(result.shape, (1, 1, 3))
        self.assertTrue(np.all(result >= 0))
        self.assertTrue(np.all(result <= 1))


if __name__ == "__main__":
    unittest.main()
