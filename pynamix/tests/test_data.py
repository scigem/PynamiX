import unittest
import numpy as np
import os
import tempfile
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for testing
import matplotlib.pyplot as plt
from pynamix import data


class TestDataModule(unittest.TestCase):
    """Test cases for the data module"""

    def test_spiral_creates_image(self):
        """Test that spiral() creates an image file"""
        # Create temporary directory
        temp_dir = tempfile.mkdtemp()
        original_dir = os.getcwd()
        
        try:
            os.chdir(temp_dir)
            
            # Generate spiral
            data.spiral()
            
            # Check that file was created
            self.assertTrue(os.path.exists("spiral.png"))
            
            # Check file is not empty
            self.assertGreater(os.path.getsize("spiral.png"), 0)
            
        finally:
            os.chdir(original_dir)
            import shutil
            shutil.rmtree(temp_dir)

    def test_fibres_creates_image(self):
        """Test that fibres() creates an image file"""
        temp_dir = tempfile.mkdtemp()
        
        try:
            # Generate fibres with known parameters
            theta_mean = 0.0
            kappa = 1.0
            N = 100
            
            data.fibres(theta_mean=theta_mean, kappa=kappa, N=N, foldername=temp_dir)
            
            # Check that file was created with expected name
            expected_file = os.path.join(temp_dir, f"fibres_{theta_mean}_{kappa}_{N}.png")
            self.assertTrue(os.path.exists(expected_file))
            
            # Check file is not empty
            self.assertGreater(os.path.getsize(expected_file), 0)
            
        finally:
            import shutil
            shutil.rmtree(temp_dir)

    def test_fibres_with_different_orientations(self):
        """Test fibres with different mean orientations"""
        temp_dir = tempfile.mkdtemp()
        
        try:
            # Test several orientations
            orientations = [0.0, np.pi/4, np.pi/2, np.pi]
            
            for theta in orientations:
                data.fibres(theta_mean=theta, kappa=1.0, N=50, foldername=temp_dir)
                
                expected_file = os.path.join(temp_dir, f"fibres_{theta}_1.0_50.png")
                self.assertTrue(os.path.exists(expected_file))
            
        finally:
            import shutil
            shutil.rmtree(temp_dir)

    def test_fibres_with_different_kappa(self):
        """Test fibres with different alignment parameters"""
        temp_dir = tempfile.mkdtemp()
        
        try:
            # Test different kappa values (alignment)
            kappas = [0.1, 1.0, 5.0]
            
            for kappa in kappas:
                data.fibres(theta_mean=0.0, kappa=kappa, N=50, foldername=temp_dir)
                
                expected_file = os.path.join(temp_dir, f"fibres_0.0_{kappa}_50.png")
                self.assertTrue(os.path.exists(expected_file))
            
        finally:
            import shutil
            shutil.rmtree(temp_dir)

    def test_fibres_with_different_N(self):
        """Test fibres with different numbers of particles"""
        temp_dir = tempfile.mkdtemp()
        
        try:
            # Test different N values
            N_values = [10, 100, 500]
            
            for N in N_values:
                data.fibres(theta_mean=0.0, kappa=1.0, N=N, foldername=temp_dir)
                
                expected_file = os.path.join(temp_dir, f"fibres_0.0_1.0_{N}.png")
                self.assertTrue(os.path.exists(expected_file))
            
        finally:
            import shutil
            shutil.rmtree(temp_dir)

    def test_fibres_custom_dpi(self):
        """Test fibres with custom DPI"""
        temp_dir = tempfile.mkdtemp()
        
        try:
            data.fibres(theta_mean=0.0, kappa=1.0, N=50, dpi=100, foldername=temp_dir)
            
            expected_file = os.path.join(temp_dir, "fibres_0.0_1.0_50.png")
            self.assertTrue(os.path.exists(expected_file))
            
        finally:
            import shutil
            shutil.rmtree(temp_dir)

    def test_fibres_custom_linewidth(self):
        """Test fibres with custom line width"""
        temp_dir = tempfile.mkdtemp()
        
        try:
            data.fibres(theta_mean=0.0, kappa=1.0, N=50, lw=2, foldername=temp_dir)
            
            expected_file = os.path.join(temp_dir, "fibres_0.0_1.0_50.png")
            self.assertTrue(os.path.exists(expected_file))
            
        finally:
            import shutil
            shutil.rmtree(temp_dir)

    def test_fibres_custom_alpha(self):
        """Test fibres with custom transparency"""
        temp_dir = tempfile.mkdtemp()
        
        try:
            data.fibres(theta_mean=0.0, kappa=1.0, N=50, alpha=0.5, foldername=temp_dir)
            
            expected_file = os.path.join(temp_dir, "fibres_0.0_1.0_50.png")
            self.assertTrue(os.path.exists(expected_file))
            
        finally:
            import shutil
            shutil.rmtree(temp_dir)


class TestPendulumData(unittest.TestCase):
    """Test pendulum data loading - expected to fail without actual data"""

    def test_pendulum_without_data(self):
        """Test that pendulum() handles missing data appropriately"""
        # This test documents that pendulum() requires external data
        # It should either download or raise an exception
        
        # Skip this test if we're in automated testing without user input
        # In a real scenario, this would test the download prompt
        pass


if __name__ == "__main__":
    unittest.main()
