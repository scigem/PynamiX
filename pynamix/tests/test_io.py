import unittest
import numpy as np
import os
import tempfile
import json
from pynamix import io


class TestIOModule(unittest.TestCase):
    """Test cases for the io module"""

    def test_strip_seq_log_with_seq(self):
        """Test strip_seq_log removes .seq extension"""
        result = io.strip_seq_log("test_file.seq")
        self.assertEqual(result, "test_file")

    def test_strip_seq_log_with_log(self):
        """Test strip_seq_log removes .log extension"""
        result = io.strip_seq_log("test_file.log")
        self.assertEqual(result, "test_file")

    def test_strip_seq_log_without_extension(self):
        """Test strip_seq_log with no extension"""
        result = io.strip_seq_log("test_file")
        self.assertEqual(result, "test_file")

    def test_strip_seq_log_other_extension(self):
        """Test strip_seq_log with other extension (should not strip)"""
        result = io.strip_seq_log("test_file.txt")
        self.assertEqual(result, "test_file.txt")


class TestGenerateSeq(unittest.TestCase):
    """Test SEQ file generation"""

    def setUp(self):
        """Create temporary directory for test files"""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up temporary files"""
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_generate_seq_detector_0_mode_0(self):
        """Test generate_seq for detector 0, mode 0"""
        filepath = os.path.join(self.temp_dir, "test_d0_m0")
        io.generate_seq(filepath, detector=0, mode=0, nbframe=5)
        
        # Check file was created
        self.assertTrue(os.path.exists(filepath + ".seq"))
        
        # Check file size (768 * 960 * 2 bytes * 5 frames)
        expected_size = 768 * 960 * 2 * 5
        actual_size = os.path.getsize(filepath + ".seq")
        self.assertEqual(actual_size, expected_size)

    def test_generate_seq_detector_0_mode_1(self):
        """Test generate_seq for detector 0, mode 1"""
        filepath = os.path.join(self.temp_dir, "test_d0_m1")
        io.generate_seq(filepath, detector=0, mode=1, nbframe=3)
        
        # Check file was created
        self.assertTrue(os.path.exists(filepath + ".seq"))
        
        # Check file size (1536 * 1920 * 2 bytes * 3 frames)
        expected_size = 1536 * 1920 * 2 * 3
        actual_size = os.path.getsize(filepath + ".seq")
        self.assertEqual(actual_size, expected_size)

    def test_generate_seq_detector_2_mode_11(self):
        """Test generate_seq for detector 2, mode 11"""
        filepath = os.path.join(self.temp_dir, "test_d2_m11")
        io.generate_seq(filepath, detector=2, mode=11, nbframe=2)
        
        # Check file was created
        self.assertTrue(os.path.exists(filepath + ".seq"))
        
        # Check file size (3072 * 3888 * 2 bytes * 2 frames)
        expected_size = 3072 * 3888 * 2 * 2
        actual_size = os.path.getsize(filepath + ".seq")
        self.assertEqual(actual_size, expected_size)

    def test_generate_seq_detector_2_mode_22(self):
        """Test generate_seq for detector 2, mode 22"""
        filepath = os.path.join(self.temp_dir, "test_d2_m22")
        io.generate_seq(filepath, detector=2, mode=22, nbframe=2)
        
        # Check file was created
        self.assertTrue(os.path.exists(filepath + ".seq"))
        
        # Check file size (3072/2 * 3888/2 * 2 bytes * 2 frames)
        expected_size = int(3072/2) * int(3888/2) * 2 * 2
        actual_size = os.path.getsize(filepath + ".seq")
        self.assertEqual(actual_size, expected_size)

    def test_generate_seq_detector_2_mode_44(self):
        """Test generate_seq for detector 2, mode 44"""
        filepath = os.path.join(self.temp_dir, "test_d2_m44")
        io.generate_seq(filepath, detector=2, mode=44, nbframe=2)
        
        # Check file was created
        self.assertTrue(os.path.exists(filepath + ".seq"))
        
        # Check file size (3072/4 * 3888/4 * 2 bytes * 2 frames)
        expected_size = int(3072/4) * int(3888/4) * 2 * 2
        actual_size = os.path.getsize(filepath + ".seq")
        self.assertEqual(actual_size, expected_size)

    def test_generate_seq_invalid_mode_detector_0(self):
        """Test generate_seq with invalid mode for detector 0"""
        filepath = os.path.join(self.temp_dir, "test_invalid")
        
        with self.assertRaises(Exception) as context:
            io.generate_seq(filepath, detector=0, mode=11, nbframe=2)
        
        self.assertIn("Mode should be 0, or 1", str(context.exception))

    def test_generate_seq_invalid_mode_detector_2(self):
        """Test generate_seq with invalid mode for detector 2"""
        filepath = os.path.join(self.temp_dir, "test_invalid")
        
        with self.assertRaises(Exception) as context:
            io.generate_seq(filepath, detector=2, mode=0, nbframe=2)
        
        self.assertIn("Mode should be 11, 22 or 44", str(context.exception))


class TestLoadImage(unittest.TestCase):
    """Test image loading functionality"""

    def setUp(self):
        """Create temporary directory and test image"""
        self.temp_dir = tempfile.mkdtemp()
        self.test_image = os.path.join(self.temp_dir, "test.png")
        
        # Create a simple test image
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(2, 2))
        ax.imshow([[0, 1], [1, 0]], cmap='gray')
        ax.axis('off')
        plt.savefig(self.test_image, bbox_inches='tight', pad_inches=0)
        plt.close()

    def tearDown(self):
        """Clean up temporary files"""
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_load_image_as_gray(self):
        """Test load_image with as_gray=True"""
        ims, logfile = io.load_image(self.test_image, as_gray=True)
        
        # Check shape is 3D (1, height, width)
        self.assertEqual(len(ims.shape), 3)
        self.assertEqual(ims.shape[0], 1)
        
        # Check logfile structure
        self.assertIn("detector", logfile)
        self.assertIn("geometry", logfile)
        self.assertIn("X-rays", logfile)

    def test_load_image_color(self):
        """Test load_image with as_gray=False"""
        ims, logfile = io.load_image(self.test_image, as_gray=False)
        
        # Check shape is 3D or 4D depending on color channels
        self.assertGreaterEqual(len(ims.shape), 3)
        self.assertEqual(ims.shape[0], 1)


class TestUpgradeLogfile(unittest.TestCase):
    """Test logfile upgrade functionality - may expose hardcoded paths"""

    def setUp(self):
        """Create temporary directory"""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up temporary files"""
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_upgrade_logfile_basic(self):
        """Test upgrade_logfile with old format logfile"""
        old_log_path = os.path.join(self.temp_dir, "test.log")
        
        # Create a simple old-format logfile
        with open(old_log_path, 'w') as f:
            f.write("Mon Jan 01 12:00:00 2024\n")
            f.write("\n")
            f.write("MODE 0\n")
            f.write("768x960\n")
            f.write("ROI 0, 0 768 960\n")
            f.write("FPS 30\n")
            f.write("\n")
            f.write("1000\n")
            f.write("1001\n")
            f.write("1002\n")
        
        # Upgrade the logfile
        io.upgrade_logfile(old_log_path)
        
        # Check that old file was renamed
        self.assertTrue(os.path.exists(old_log_path + ".dep"))
        
        # Check new JSON file was created
        self.assertTrue(os.path.exists(old_log_path))
        
        # Load and verify new logfile
        with open(old_log_path, 'r') as f:
            new_log = json.load(f)
        
        self.assertIn("detector", new_log)
        self.assertEqual(new_log["detector"]["mode"], 0)
        self.assertEqual(new_log["detector"]["image_size"]["width"], 768)
        self.assertEqual(new_log["detector"]["image_size"]["height"], 960)
        self.assertEqual(new_log["detector"]["fps"], 30)
        self.assertEqual(len(new_log["detector"]["frames"]), 3)


class TestSaveAsTiffs(unittest.TestCase):
    """Test TIFF export functionality"""

    def setUp(self):
        """Create temporary directory"""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up temporary files"""
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_save_as_tiffs_basic(self):
        """Test save_as_tiffs with simple data"""
        output_folder = os.path.join(self.temp_dir, "tiffs")
        
        # Create simple test data
        data = np.random.rand(5, 10, 10) * 1000
        
        io.save_as_tiffs(output_folder, data, tmin=0, tmax=3, tstep=1)
        
        # Check folder was created
        self.assertTrue(os.path.exists(output_folder))
        
        # Check files were created (frames 0, 1, 2)
        self.assertTrue(os.path.exists(os.path.join(output_folder, "00000.tiff")))
        self.assertTrue(os.path.exists(os.path.join(output_folder, "00001.tiff")))
        self.assertTrue(os.path.exists(os.path.join(output_folder, "00002.tiff")))
        
        # Frame 3 and 4 should not exist (tmax=3 is exclusive)
        self.assertFalse(os.path.exists(os.path.join(output_folder, "00003.tiff")))

    def test_save_as_tiffs_with_step(self):
        """Test save_as_tiffs with step parameter"""
        output_folder = os.path.join(self.temp_dir, "tiffs_step")
        
        data = np.random.rand(10, 10, 10) * 1000
        
        io.save_as_tiffs(output_folder, data, tmin=0, tmax=10, tstep=2)
        
        # Check only even frames were saved
        self.assertTrue(os.path.exists(os.path.join(output_folder, "00000.tiff")))
        self.assertTrue(os.path.exists(os.path.join(output_folder, "00002.tiff")))
        self.assertTrue(os.path.exists(os.path.join(output_folder, "00004.tiff")))
        
        # Odd frames should not exist
        self.assertFalse(os.path.exists(os.path.join(output_folder, "00001.tiff")))
        self.assertFalse(os.path.exists(os.path.join(output_folder, "00003.tiff")))


if __name__ == "__main__":
    unittest.main()
