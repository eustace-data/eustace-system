"""Compare fields for ocean."""

from comparefields import comparefields
import argparse
import datetime
import os
import eustaceconfig

BASELINE_DIRECTORY=os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/MS8_version2/Ocean/data/v.5.0')
"""Location of comparison data."""

def comparefields_ocean(baseline_filename,  primary_filename, ancillary_filename):
    """Run comparison given the named files."""

    # List for any issues
    issues = [ ]

    # field map comparison -> baseline field name
    fieldmap_primary = {
        'tas': 'tas',
    }
    fieldmap_ancillary = {
    'tas_unc_rand': 'unc_rand_tas',
    'tas_unc_corr_sat': 'unc_corr_tas',
    'tas_unc_sys': 'unc_syst_tas',
    'tas_unc_corr_mod': 'unc_corr2_tas',
    'tas_unc_sys_mod': 'unc_syst2_tas',
    'tas_unc_parameter_0': 'unc_parameter_0_tas',
    'tas_unc_parameter_1': 'unc_parameter_1_tas',
    'tas_unc_parameter_2': 'unc_parameter_2_tas',
    'tas_unc_parameter_3': 'unc_parameter_3_tas',
    'tas_unc_parameter_4': 'unc_parameter_4_tas'
    }

    # Run comparison (primary file)
    issues += comparefields(baseline_filename,
                            primary_filename,
                            fieldmap_primary,
                            offsets={ 'tas': 273.15 })

    # Run comparison (ancillary file)
    issues += comparefields(baseline_filename,
                            ancillary_filename,
                            fieldmap_ancillary)

    return issues


def main():
    """Run on satstace output from command line."""

    # Parse commandline
    parser = argparse.ArgumentParser('comparefields_ocean')
    parser.add_argument('satstacepath', help='base path of satstace output')
    parser.add_argument('datestring', help='date to process YYYYmmdd')
    args = parser.parse_args()

    # Date string from args
    datestring = args.datestring

    # Year string
    yearstring = args.datestring[0:4]

    # Directories for specified year
    baseline_year_directory = os.path.join(BASELINE_DIRECTORY, yearstring)
    comparison_year_directory = os.path.join(args.satstacepath, os.path.join(os.path.join('ocean', yearstring)))

    # Baseline filename
    baseline_filename = os.path.join(baseline_year_directory, 'eustace_AATSR_satellite_ocean_{datestring}.nc'.format(datestring=datestring))

    # Comparison filename for primary variables
    primary_filename = os.path.join(comparison_year_directory, 'tas_ocean_eustace_0_{datestring}.nc'.format(datestring=datestring))

    # Comparison filename for ancillary variables
    ancillary_filename = os.path.join(comparison_year_directory, 'ancillary_ocean_eustace_0_{datestring}.nc'.format(datestring=datestring))

    # Run comparison
    issues = comparefields_ocean(baseline_filename, primary_filename, ancillary_filename)

    # Results
    if issues:
        for issue in issues:
            print issue
    else:
        print 'OK'

if __name__=='__main__':
    main()
