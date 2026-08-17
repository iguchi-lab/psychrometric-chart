import tempfile
import unittest
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")

from psychrometric_chart import (  # noqa: E402
    CP_DA,
    ChartConfig,
    chart_temperature,
    create_chart,
    humidity_ratio,
    physical_state_mask,
    saturation_pressure,
)


class PsychrometricChartTests(unittest.TestCase):
    def test_water_is_the_default_below_freezing(self):
        self.assertAlmostEqual(
            saturation_pressure(-10.0),
            saturation_pressure(-10.0, "water"),
        )
        self.assertGreater(
            saturation_pressure(-10.0, "water"),
            saturation_pressure(-10.0, "ice"),
        )

    def test_water_and_ice_are_equal_at_and_above_zero(self):
        for temperature in (0.0, 20.0):
            self.assertAlmostEqual(
                saturation_pressure(temperature, "water"),
                saturation_pressure(temperature, "ice"),
            )

    def test_lower_pressure_increases_humidity_ratio(self):
        sea_level = humidity_ratio(25.0, 50.0, 101.325)
        high_altitude = humidity_ratio(25.0, 50.0, 80.0)
        self.assertGreater(high_altitude, sea_level)

    def test_chart_temperature_is_finite_at_old_singularity(self):
        value = chart_temperature(50.0, 0.0)
        self.assertTrue(np.isfinite(value))
        self.assertEqual(value, 50.0)

    def test_dry_bulb_line_leans_left_below_reference_temperature(self):
        temperature = 20.0
        low_humidity = chart_temperature(temperature, 0.005)
        high_humidity = chart_temperature(temperature, 0.020)
        self.assertLess(high_humidity, low_humidity)

    def test_chart_transform_matches_original_formula(self):
        temperature = 26.0
        humidity = 0.010
        enthalpy_value = CP_DA * temperature + humidity * (2501.0 + 1.86 * temperature)
        x_at_50 = (enthalpy_value - CP_DA * 50.0) / (2501.0 + 1.86 * 50.0)
        t_at_zero_humidity = enthalpy_value / CP_DA
        slope = x_at_50 / (50.0 - t_at_zero_humidity)
        original = 50.0 - (x_at_50 - humidity) / slope
        self.assertAlmostEqual(chart_temperature(temperature, humidity), original)

    def test_physical_mask_rejects_supersaturation(self):
        config = ChartConfig()
        temperature = np.array([20.0, 20.0])
        saturation = humidity_ratio(20.0, 100.0)
        humidity = np.array([saturation * 0.9, saturation * 1.1])
        np.testing.assert_array_equal(
            physical_state_mask(temperature, humidity, config),
            np.array([True, False]),
        )

    def test_svg_is_generated_for_both_subzero_modes(self):
        with tempfile.TemporaryDirectory() as directory:
            for surface in ("water", "ice"):
                figure, _ = create_chart(ChartConfig(subzero_surface=surface))
                output = Path(directory) / f"chart-{surface}.svg"
                figure.savefig(output)
                self.assertGreater(output.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
