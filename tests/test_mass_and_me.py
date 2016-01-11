#!/usr/bin/env python3

import mass_and_me
from nonstdlib import approx

def test_load_expected_masses():
    csv_path = 'sample_data/expected_masses.csv'
    assert mass_and_me.load_expected_masses(csv_path) == [1303, 1283]

def test_load_mass_spectra():
    csv_path = 'sample_data/mass_spectrum.csv'
    ms = mass_and_me.load_mass_spectra([csv_path])[0]

    assert ms.Mass[0] == approx(396.009089)
    assert ms.Intensity[0] == approx(0)
    assert ms.Mass[12] == approx(400.015755)
    assert ms.Intensity[12] == approx(49.810793)

def test_find_peak_indices():
    csv_path = 'sample_data/mass_spectrum.csv'
    ms = mass_and_me.load_mass_spectra([csv_path])[0]
    all_peaks = mass_and_me.find_peak_indices(ms)
    big_peaks = mass_and_me.find_peak_indices(ms, 100)

    assert all_peaks[0] == 12
    assert all_peaks[1] == 17;  assert big_peaks[0] == 17
    assert all_peaks[2] == 32
    assert all_peaks[3] == 46
    assert all_peaks[4] == 60
    assert all_peaks[5] == 73
    assert all_peaks[6] == 91
    assert all_peaks[7] == 102; assert big_peaks[1] == 102

