"""Tests for the curve-fitting functions vendored from pytoxr.

These two functions were imported from pytoxr until the dependency was removed and the code
copied into thoipapy.utils. pytoxr had no test suite, and thoipapy's only caller is
paper_figures/retrospective.py, which no test exercises. They therefore arrived here with no
coverage at all.

What is worth pinning down is the constraint that gives sine_perfect_helix its name: a and b are
hardcoded so the curve keeps alpha-helical periodicity of 3.6 residues per turn, and only the
phase and the vertical offset are fitted. A change to either constant would still produce a
plausible-looking curve on the figure while no longer representing an alpha helix.
"""

import numpy as np
from scipy.optimize import leastsq
from thoipapy.utils import residuals, sine_perfect_helix


def test_residuals_are_zero_for_a_perfect_fit():
    def straight_line(constants, x):
        return constants[0] * x

    x = np.array([1.0, 2.0, 3.0])
    y = np.array([2.0, 4.0, 6.0])
    np.testing.assert_allclose(residuals([2.0], straight_line, x, y), [0.0, 0.0, 0.0])


def test_residuals_return_the_signed_distance_from_the_fitted_curve():
    def straight_line(constants, x):
        return constants[0] * x

    x = np.array([1.0, 2.0, 3.0])
    y = np.array([2.0, 4.0, 6.0])
    # a slope of 1 underestimates every point, so every residual is positive
    np.testing.assert_allclose(residuals([1.0], straight_line, x, y), [1.0, 2.0, 3.0])


def test_sine_perfect_helix_has_a_periodicity_of_3_6_residues():
    """The whole point of the function: b is fixed at 1.745, so 2 * pi / b is 3.6 residues."""
    x = np.linspace(0, 40, 40001)
    y = sine_perfect_helix([0.0, 0.0], x)
    # distance between successive maxima, in residues
    peaks = x[1:-1][(y[1:-1] > y[:-2]) & (y[1:-1] > y[2:])]
    periods = np.diff(peaks)
    assert len(periods) > 5
    np.testing.assert_allclose(periods, 3.6, atol=0.01)


def test_sine_perfect_helix_amplitude_and_offset_are_fixed():
    """a is fixed at 0.2, giving a curve 0.4 in height, centred on the fitted offset d."""
    x = np.linspace(0, 40, 4001)
    for offset in (0.0, 0.5, -2.0):
        y = sine_perfect_helix([0.0, offset], x)
        np.testing.assert_allclose(y.max() - y.min(), 0.4, atol=1e-4)
        np.testing.assert_allclose((y.max() + y.min()) / 2, offset, atol=1e-4)


def test_only_the_phase_and_offset_shift_the_curve():
    x = np.linspace(0, 40, 4001)
    baseline = sine_perfect_helix([0.0, 0.0], x)
    # a full period of phase shift returns the identical curve
    np.testing.assert_allclose(sine_perfect_helix([2 * np.pi, 0.0], x), baseline, atol=1e-12)
    # the offset is a pure vertical translation
    np.testing.assert_allclose(sine_perfect_helix([0.0, 1.5], x), baseline + 1.5, atol=1e-12)


def test_leastsq_fit_used_by_the_retrospective_figure():
    """Locks the fit at the one call site, so the plotted sine wave cannot drift silently."""
    x_helical = np.ogrid[0:23.4:14j]
    y_helical = [1, 0] * 7
    sine_constants_guess_perfhelix = [1.575, 0.5]

    fitted = leastsq(
        residuals, sine_constants_guess_perfhelix, args=(sine_perfect_helix, x_helical, y_helical), full_output=1
    )[0]

    np.testing.assert_allclose(fitted, [1.57464857, 0.5], atol=1e-6)
