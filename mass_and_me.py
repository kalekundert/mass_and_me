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

__version__ = '0.0.0'

## Inputs
# ======
# List of expected masses

## Goal
# ====
# Rank those by how well they are supported by the data.

## Questions
# =========
# 1. Is the first column mass or mass/charge?
#
#    It's mass/charge.
#
# 2. How does the HPLC data play into things?
#
#    Fitzy can ask the mass spec to bin the data into particular time windows.  
#    The sample data she sent me is a ten minute window, I think.  It would be 
#    nice if this program could handle the binning, though, to save the manual 
#    effort of creating 5-6 windows.  To do that, each mass would need to be 
#    annotated with a time.  That may or may not be possible, though.
#

## Algorithm
# =========
# 1. Identify peaks in the mass spec data:
#   
#   a. Look for down/up/down patterns
#
#   b. Take intensity from "up" positions

# 2. For each expected mass, note the intensity of the nearest peak (or 0 if 
#    there are no peaks nearby).
#
# 3. For each expected mass, note if m/2, m/3, etc. peaks are present.
#
# 4. For each expected mass, note if an ion series is present (+1 +2 +3 for 
#    m/1, +1/2 +2/2 +3/2 for m/2, etc.)
#
# 5. For each expected mass, note if a methyl modification is present (+14 for 
#    m/1, +7 for m/2, +3.5 for m/3, etc.)
#
# 6. Sort each expected mass by how well it fulfills the above criteria.

Z_SERIES = 1, 2, 3
ION_SERIES = 1, 2, 3, 4


class ExpectedMass:

    def __init__(self, mass):
        self.mass = mass
        self.mz_intensities = {}
        self.ion_intensities = {}
        self.methyl_intensities = {}

    def __repr__(self):
        return 'ExpectedMass({})'.format(self.mass)

    def __str__(self):
        return """\
mass: {}
m/z intensities:
    {}
CH₃ intensities:
    {}
""".format(
            self.mass,
            '\n    '.join('z={}: {:12.3f}'.format(k,v) for k,v in self.mz_intensities.items()),
            '\n    '.join('z={}: {:12.3f}'.format(k,v) for k,v in self.methyl_intensities.items()),
        )

    def __eq__(self, other):
        return self.mass == other.mass

    def __lt__(self, other):
        """
        Return True is this mass appears to be better supported by the data 
        than the other mass.
        """
        return self.mz_intensities[1] > other.mz_intensities[1]

    def calculate_quality_metrics(self, ms, peaks):
        # Find peaks that correspond to this mass but that have different 
        # charges (z).  The code here keeps looking for peaks corresponding to 
        # higher charges until it stops finding peaks.

        self.mz_intensities = {1: 0}

        for z in itertools.count(1):
            intensity = find_intensity(ms, peaks, self.mass, z)
            if intensity > 0:
                self.mz_intensities[z] = intensity
            else:
                break

        # Find peaks that correspond to the addition of one or more protons to 
        # the expected mass.  There are typically 4-5 such peaks, each one 
        # decreasing in intensity.  The code here will keep looking for peaks
        # that are heavier by the mass of a proton until it finds one that 
        # increases in intensity, which is assumed to be a new species.

        self.ion_intensities = {
                z: []
                for z in self.mz_intensities
        }
        for z in self.mz_intensities:
            prev_mass = self.mass
            prev_intensity = self.mz_intensities[z]

            while True:
                next_mass = prev_mass + 1
                next_intensity = find_intensity(ms, peaks, next_mass, z)

                if next_intensity >= prev_intensity:
                    break
                else:
                    self.ion_intensities[z].append(next_intensity)
                    prev_mass = next_mass
                    prev_intensity = next_intensity

        # Find peaks that correspond to the addition of a methyl group to the 
        # expected mass.  Methyl has a mass or 14 (1 carbon + 3 hydrogens).
                
        self.methyl_intensities = {
                z: find_intensity(ms, peaks, self.mass + 14, z)
                for z in self.mz_intensities
        }



def main():
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
    df = pd.read_csv(csv_path, header=None)
    return list(df[0])

def load_mass_spectra(csv_paths):
    return [pd.read_csv(x, header=7) for x in csv_paths]

def generate_report(mass_spectra, expected_masses, peak_threshold):
    for ms in mass_spectra:
        peaks = find_peaks(ms, peak_threshold)
        masses = find_masses(ms, peaks, expected_masses)
        for mass in masses:
            print(mass)

def plot_mass_spectra(mass_spectra):
    import os, pylab
    import multiprocessing as mp

    for ms in mass_spectra:
        ms.plot(x='Mass', y='Intensity')

    pylab.show()

def find_peaks(ms, peak_threshold):
    from scipy.signal import argrelextrema
    peaks = argrelextrema(ms.Intensity.values, np.greater)[0]
    return peaks[ms.Intensity.values[peaks] > peak_threshold]

def find_masses(ms, peaks, expected_masses):
    found_masses = [ExpectedMass(x) for x in expected_masses]
    for found_mass in found_masses:
        found_mass.calculate_quality_metrics(ms, peaks)
    return sorted(found_masses)

def find_intensity(ms, peaks, mass, z=1):
    """
    Find the peak that corresponds to the given mass. 
    """
    
    # I thought using a binary search to implement this function, which 
    # would've been O(log(N)) because the MS data is pre-sorted, but I decided 
    # that it's probably faster to just use numpy.

    mass_diffs = np.abs(mass / z - ms.Mass.values[peaks])
    peak_i = np.argmin(mass_diffs)

    # If the closest peak is further than 1 M/Z unit from the desired mass, 
    # then the peak we're looking for isn't in the data.  Indicate this by 
    # returning an intensity of zero.

    if mass_diffs[peak_i] > 1 / z:
        return 0
    else:
        return ms.Intensity[peaks[peak_i]]


if __name__ == '__main__':
    main()


