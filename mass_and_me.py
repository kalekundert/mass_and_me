#!/usr/bin/env python

"""\
Search for expected peaks in mass spectrometry data.

Usage:
    mass_plus_me <expected_masses_csv> <mass_spectra_csv>... [options]
    mass_plus_me plot <mass_spectra_csv>...

Arguments:
    <expected_masses_csv>
        A *.csv file (as can be exported from Excel) containing a list of the 
        expected masses in its first column.  This program will generate a 
        report detailing how well each of these masses is supported by the 
        given data.

    <mass_spectra_csv>
        One or more *.csv files (as can be exported from the instrument) 
        containing masses in a column labeled 'Mass' and the corresponding 
        intensities in a column labeled 'Intensity'.

Options:
    -p, --peak-threshold=<intensity>        [default: 0]
        The minimum intensity required for a peak to be considered.

    -v, --visualize
        Plot the mass spectra being analyzed in a new window.
"""

import itertools
import numpy as np
import pandas as pd
from nonstdlib import *

__version__ = '0.0.0'

class Peak:
    """
    An object representing a single peak in the mass spec data.  Each peak 
    contains two pieces of data: a mass and an intensity.
    """

    def __init__(self, mass, intensity):
        self.mass = mass
        self.intensity = intensity

    def __repr__(self):
        """
        Return a brief string describing this object.
        """
        return 'Peak({self.mass}, {self.intensity})' | fmt


class ExpectedMass:
    """
    An object representing a mass that the user is expecting to find in the 
    data.  Given a mass, this object can calculate quality metrics indicating 
    how well that mass is supported by the data.  More specifically, these 
    quality metrics comprise the intensities of peaks corresponding to the same 
    mass with different charges, to the addition of one or more protons, and to 
    modification by a methyl group.
    """

    def __init__(self, mass):
        self.mass = mass
        self.mz_peaks = {}
        self.ion_peaks = {}
        self.methyl_peaks = {}

    def __repr__(self):
        """
        Return a brief string describing this object.
        """
        return 'ExpectedMass({self.mass})' | fmt

    def __eq__(self, other):
        """
        Two ``ExpectedMass`` objects are equal if they have the same mass.
        """
        return self.mass == other.mass

    def __lt__(self, other):
        """
        One ``ExpectedMass`` is regarded as "less than" another if it's better 
        supported by the data.  This is useful for sorting, because the default 
        sorting algorithms in python put the smallest items first.
        """
        return self.mz_peaks[1].intensity > other.mz_peaks[1].intensity

    def calculate_quality_metrics(self, ms, peaks):
        """
        Look for the ancillary peaks that would indicate that expected mass is 
        supported by the data.
        """
        self._calculate_mz_peaks(ms, peaks)
        self._calculate_ion_peaks(ms, peaks)
        self._calculate_methyl_peaks(ms, peaks)

    def _calculate_mz_peaks(self, ms, peaks):
        """
        Find peaks that correspond to this mass but that have different 
        charges (z).  The code here keeps looking for peaks corresponding to 
        higher charges until it stops finding peaks.
        """
        self.mz_peaks = {1: 0}

        for z in itertools.count(1):
            peak = find_peak(ms, peaks, self.mass, z)
            if peak.intensity > 0:
                self.mz_peaks[z] = peak
            else:
                break

    def _calculate_ion_peaks(self, ms, peaks):
        """
        Find peaks that correspond to the addition of one or more protons to 
        the expected mass.  There are typically 4-5 such peaks, each one 
        decreasing in intensity.  The code here will keep looking for peaks
        that are heavier by the mass of a proton until it finds one that 
        increases in intensity, which is assumed to be a new species.
        """
        self.ion_peaks = {
                z: []
                for z in self.mz_peaks
        }
        for z in self.mz_peaks:
            prev_peak = self.mz_peaks[z]
            while True:
                this_peak = find_peak(ms, peaks, prev_peak.mass + 1, z)
                if this_peak.intensity >= prev_peak.intensity:
                    break
                else:
                    self.ion_peaks[z].append(this_peak)
                    prev_peak = this_peak

    def _calculate_methyl_peaks(self, ms, peaks):
        """
        Find peaks that correspond to the addition of a methyl group to the 
        expected mass.  Methyl has a mass or 14 (1 carbon + 3 hydrogens).
        """
        self.methyl_peaks = {
                z: find_peak(ms, peaks, self.mass + 14, z)
                for z in self.mz_peaks
        }

    def report_quality_metrics(self):
        """
        Print a brief summary of the previously calculated quality metric, to 
        assist the user in deciding which peaks to manually verify.
        """
        print_color('mass: {self.mass}' | fmt, 'magenta', 'bold')

        print(' m/z peaks:')
        for z, peak in self.mz_peaks.items():
            print('   z={}, m={:7.2f}, {}'.format(
                z, peak.mass, format_intensity(peak.intensity)))

        print(' ion series:')
        for z, peaks in self.ion_peaks.items():
            print('   z={}, m={:7.2f}(+n/{}), {}'.format(
                z, self.mass/z, z,
                ', '.join(format_intensity(x.intensity) for x in peaks)))

        print(' CH₃ peaks:')
        for z, peak in self.methyl_peaks.items():
            print('   z={}, m={:7.2f}, {}'.format(
                z, peak.mass, format_intensity(peak.intensity)))

        print()



def main():
    """
    Help identify promising peaks in mass spectrometry data.  Parse the options 
    provided by the user on the command line and then carry out the requested 
    analysis.
    """
    import docopt
    args = docopt.docopt(__doc__)
    expected_masses = load_expected_masses(args['<expected_masses_csv>'])
    mass_spectra = load_mass_spectra(args['<mass_spectra_csv>'])
    generate_report(
            mass_spectra,
            expected_masses,
            peak_threshold=int(args['--peak-threshold']),
    )
    if args['--visualize']:
        plot_mass_spectra(mass_spectra)

def load_expected_masses(csv_path):
    """
    Read the given ``*.csv`` file and return a list of all the neutral masses 
    found in the first column.  The file must contain nothing but masses (i.e. 
    no headers or anything).  Masses outside the first column will be ignored.
    """
    df = pd.read_csv(csv_path, header=None)
    return list(df[0])

def load_mass_spectra(csv_paths):
    """
    Read the given mass spectrometer ``*.csv`` output files and convert each
    one into a ``pandas.DataFrame`` with two columns: "Mass" and "Intensity".  
    The ``*.csv`` file must contain masses in the first column and intensities 
    in the second.  The data must start on the 8th row, and the seventh row 
    must contain the headers "Mass" and "Intensity" in the first and second 
    columns, respectively.
    """
    return [pd.read_csv(x, header=7) for x in csv_paths]

def generate_report(mass_spectra, expected_masses, peak_threshold=None):
    """
    For each mass spectra, search for peaks corresponding to the expected 
    masses and print out a report summarizing how well each of those masses is 
    supported by the data.
    """
    for ms in mass_spectra:
        peaks = find_peak_indices(ms, peak_threshold)
        masses = find_expected_masses(ms, peaks, expected_masses)
        for mass in masses:
            mass.report_quality_metrics()

def plot_mass_spectra(mass_spectra):
    """
    Display the data being analyzed in a separate window.
    """
    import os, pylab
    import multiprocessing as mp

    for ms in mass_spectra:
        ms.plot(x='Mass', y='Intensity')

    pylab.show()

def find_peak_indices(ms, peak_threshold=None):
    """
    Return a numpy array containing the indices of all the local maxima in the 
    given mass spectrum ``ms``, excluding those maxima with intensities less 
    ``peak_threshold``.  If ``peak_threshold`` isn't specified, all local 
    maxima are included.
    """
    from scipy.signal import argrelextrema
    peaks = argrelextrema(ms.Intensity.values, np.greater)[0]
    return peaks[ms.Intensity.values[peaks] > (peak_threshold or 0)]

def find_expected_masses(ms, peaks, expected_masses):
    """
    Search for the given expected masses in the data.  Return a list of hits 
    sorted such that the most well-supported masses are first.
    """
    found_masses = [ExpectedMass(x) for x in expected_masses]
    for found_mass in found_masses:
        found_mass.calculate_quality_metrics(ms, peaks)
    return sorted(found_masses)

def find_peak(ms, peaks, mass, z=1):
    """
    Find the peak that corresponds to the given mass. 
    """

    # I thought about using a binary search to implement this function, which 
    # would've been O(log(N)) because the MS data is pre-sorted, but I decided 
    # that it's probably faster to just use numpy.

    charged_mass = mass / z
    mass_diffs = np.abs(charged_mass - ms.Mass.values[peaks])
    peak_i = np.argmin(mass_diffs)

    # If the closest peak is further than 1 M/Z unit from the desired mass, 
    # then the peak we're looking for isn't in the data.  Indicate this by 
    # returning an intensity of zero.

    if mass_diffs[peak_i] > 1 / z:
        return Peak(charged_mass, 0)
    else:
        i = peaks[peak_i]
        return Peak(ms.Mass[i], ms.Intensity[i])

def format_intensity(intensity):
    """
    Make the given intensity easy for the user to understand.  First, convert 
    it to scientific notation following the convention for mass spectrometry 
    intensities.  Second, color the number by its exponent to make it easy for 
    the user to distinguish the largest values.
    """
    if intensity > 1e5: hilite = 'red'
    elif intensity > 1e4: hilite = 'yellow'
    else: hilite = 'normal'
    return color('{intensity:.3e}' | fmt, hilite)


if __name__ == '__main__':
    main()


