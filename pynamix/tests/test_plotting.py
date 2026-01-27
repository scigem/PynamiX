import unittest
import numpy as np
import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend for testing

import matplotlib.pyplot as plt
from pynamix import plotting


class TestPlottingModule(unittest.TestCase):
    """Test cases for the plotting module"""

    def test_hist_basic_3d_data(self):
        """Test hist with 3D data"""
        # Create simple test data
        data = np.random.randint(0, 65535, size=(10, 20, 20))

        # Should not raise error
        try:
            plotting.hist(data, frame=5, vmin=1000, vmax=50000)
            plt.close("all")
        except Exception as e:
            self.fail(f"hist() raised {e} unexpectedly")

    def test_hist_2d_data(self):
        """Test hist with 2D data (single frame)"""
        # Create simple 2D test data
        data = np.random.randint(0, 65535, size=(20, 20))
        # Expand to 3D for hist function
        data_3d = np.expand_dims(data, axis=0)

        try:
            plotting.hist(data_3d, frame=0, vmin=1000, vmax=50000)
            plt.close("all")
        except Exception as e:
            self.fail(f"hist() raised {e} unexpectedly")

    def test_hist_creates_figure(self):
        """Test that hist creates a matplotlib figure"""
        data = np.random.randint(0, 65535, size=(5, 20, 20))

        # Clear any existing figures
        plt.close("all")

        plotting.hist(data, frame=2, vmin=1000, vmax=50000)

        # Check that figure 99 was created
        self.assertIn(99, plt.get_fignums())

        plt.close("all")

    def test_hist_GUI_3d_data(self):
        """Test hist_GUI with 3D data"""
        data = np.random.randint(0, 65535, size=(10, 20, 20))

        # Should return an interactive widget
        try:
            widget = plotting.hist_GUI(data, vmin=1000, vmax=50000)
            # Widget should have some attributes
            self.assertIsNotNone(widget)
        except Exception as e:
            self.fail(f"hist_GUI() raised {e} unexpectedly")

    def test_hist_GUI_2d_data(self):
        """Test hist_GUI with 2D data"""
        data = np.random.randint(0, 65535, size=(20, 20))

        try:
            widget = plotting.hist_GUI(data, vmin=1000, vmax=50000)
            self.assertIsNotNone(widget)
        except Exception as e:
            self.fail(f"hist_GUI() raised {e} unexpectedly")

    def test_hist_with_various_ranges(self):
        """Test hist with different vmin/vmax ranges"""
        data = np.random.randint(0, 65535, size=(5, 20, 20))

        # Test with different ranges
        ranges = [
            (0, 65535),
            (1000, 50000),
            (10000, 20000),
        ]

        for vmin, vmax in ranges:
            try:
                plotting.hist(data, frame=2, vmin=vmin, vmax=vmax)
                plt.close("all")
            except Exception as e:
                self.fail(f"hist() with range ({vmin}, {vmax}) raised {e}")


if __name__ == "__main__":
    unittest.main()
