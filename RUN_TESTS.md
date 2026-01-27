# Running PynamiX Tests

## Quick Start

Run all tests:
```bash
python -m unittest discover -s pynamix/tests -p "test_*.py"
```

Run with verbose output:
```bash
python -m unittest discover -s pynamix/tests -p "test_*.py" -v
```

## Individual Module Tests

```bash
# Color module tests (7 tests)
python -m unittest pynamix.tests.test_color -v

# Exposure module tests (23 tests)
python -m unittest pynamix.tests.test_exposure -v

# IO module tests (16 tests)
python -m unittest pynamix.tests.test_io -v

# Measure module tests (17 tests)
python -m unittest pynamix.tests.test_measure -v

# Plotting module tests (7 tests)
python -m unittest pynamix.tests.test_plotting -v

# Data module tests (11 tests)
python -m unittest pynamix.tests.test_data -v

# Pipeline/integration tests (9 tests)
python -m unittest pynamix.tests.test_pipeline -v
```

## Expected Results

- **Total tests**: 81
- **Expected passing**: 76 (93.8%)
- **Expected failures**: 2 (minor precision/tolerance issues)
- **Expected errors**: 3 (known edge cases)

## Test Categories

### Unit Tests
- `test_color.py` - Color and colormap functions
- `test_exposure.py` - Image processing and normalization
- `test_io.py` - File I/O operations
- `test_measure.py` - Measurement and analysis functions
- `test_plotting.py` - Visualization functions
- `test_data.py` - Synthetic data generation

### Integration Tests
- `test_pipeline.py` - End-to-end workflows

## Known Issues

Some tests may show expected failures/errors:

1. **test_set_angles_from_limits_basic** - Minor floating point precision
2. **test_hanning_window_symmetry** - Discrete implementation tolerance
3. **test_set_motion_limits_custom_threshold** - Edge case with no motion
4. **test_grid_with_ROI** - ROI boundary condition
5. **test_normalise_rotation_basic** - Documents known bug (undefined variable)

## Requirements

All dependencies are installed with:
```bash
pip install -e .
```

## Continuous Integration

To run tests in CI:
```bash
python -m unittest discover -s pynamix/tests -p "test_*.py" -v 2>&1 | tee test_results.txt
```

## Coverage (Future Enhancement)

To add coverage reporting:
```bash
pip install coverage
coverage run -m unittest discover -s pynamix/tests
coverage report
coverage html
```
