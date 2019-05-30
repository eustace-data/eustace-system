"""Compare fields for land."""

from comparefields import comparefields
import argparse
import datetime
import os
import eustaceconfig

BASELINE_DIRECTORY=os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst/test_0.25_0.25')
"""Location of comparison data."""

def comparefields_land(baseline_filename,  primary_filename, ancillary_filename):
    """Run comparison given the named files."""

    issues = [ ]

    # field map comparison -> baseline field name
    fieldmap_primary = {
        'tasmin': 'tasmin',
        'tasmax': 'tasmax',
    }
    fieldmap_ancillary = {
        'tasmin_unc_rand' : 'unc_rand_tasmin',
        'tasmin_unc_corr_atm' : 'unc_corr_atm_tasmin',
	'tasmin_unc_corr_sfc' : 'unc_corr_sfc_tasmin',
        'tasmin_unc_sys' : 'unc_sys_tasmin',
        'tasmax_unc_rand': 'unc_rand_tasmax',
        'tasmax_unc_corr_atm': 'unc_corr_atm_tasmax',
	'tasmax_unc_corr_sfc': 'unc_corr_sfc_tasmax',
        'tasmax_unc_sys': 'unc_sys_tasmax',
        'tasmin_model_number': 'tasmin_model_number',
        'tasmax_model_number': 'tasmax_model_number',
    }


    issues += comparefields(baseline_filename,
                            primary_filename,
                            fieldmap_primary)

    # Run comparison (ancillary file)
    issues += comparefields(baseline_filename,
                            ancillary_filename,
                            fieldmap_ancillary)
    return issues

def main():
    """Run from command line."""

    # Parse commandline
    parser = argparse.ArgumentParser('comparefields_land')
    parser.add_argument('satstacepath', help='base path of satstace output')
    parser.add_argument('datestring', help='date to process YYYYmmdd')
    args = parser.parse_args()
    
        # Date string from args
    datestring = args.datestring

    # Year string
    yearstring = args.datestring[0:4]

    # Directories for specified year
    comparison_year_directory = os.path.join(args.satstacepath, os.path.join(os.path.join('land', yearstring)))

    # Baseline filename
    baseline_filename = os.path.join(BASELINE_DIRECTORY, 'eustace_satellite_4.100_{datestring}.nc'.format(datestring=datestring))

    # Comparison filename for primary variables
    primary_filename = os.path.join(comparison_year_directory, 'tas_land_eustace_0_{datestring}.nc'.format(datestring=datestring))

    # Comparison filename for ancillary variables
    ancillary_filename = os.path.join(comparison_year_directory, 'ancillary_land_eustace_0_{datestring}.nc'.format(datestring=datestring))

    issues = comparefields_land(baseline_filename, primary_filename, ancillary_filename)

    # Results
    if issues:
        for issue in issues:
            print issue
    else:
        print 'OK'

if __name__=='__main__':
    main()

