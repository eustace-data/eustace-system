"""Compare fields for ice."""

from comparefields import comparefields
import argparse
import datetime
import os
import eustaceconfig

BASELINE_DIRECTORY=os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/iat_from_ist/EUSTACE_satstace_v1')
"""Location of comparison data."""

def main():
    """Run from command line."""

    # Parse commandline
    parser = argparse.ArgumentParser('comparefields_ice')
    parser.add_argument('satstacepath', help='base path of satstace output')
    parser.add_argument('datestring', help='date to process YYYYmmdd')
    args = parser.parse_args()

    # List for any issues
    issues = [ ]

    # field map comparison -> baseline field name
    fieldmap_primary = {
        'tas': 'tas',
        'tasmin': 'tasmin',
        'tasmax': 'tasmax',
        'tasuncertainty': 'TotUC',
        'tasminuncertainty': 'TotUCmin',
        'tasmaxuncertainty': 'TotUCmax',
    }
    fieldmap_ancillary = {
        'tas_unc_no_cloud': 'TotU',
        'tas_unc_rand':'RU',
        'tas_unc_corr_local': 'SSU',
        'tas_unc_sys': 'LSU',
        'tas_unc_cloud': 'CU',
        'tasmin_unc_no_cloud': 'TotUmin',
        'tasmin_unc_rand':'RUmin',
        'tasmin_unc_corr_local': 'SSUmin',
        'tasmin_unc_sys': 'LSUmin',
        'tasmin_unc_cloud': 'CUmin',
        'tasmax_unc_no_cloud': 'TotUmax',
        'tasmax_unc_rand':'RUmax',
        'tasmax_unc_corr_local': 'SSUmax',
        'tasmax_unc_sys': 'LSUmax',
        'tasmax_unc_cloud': 'CUmax',
    }

    # Date string from args
    datestring = args.datestring

    # Year string
    yearstring = args.datestring[0:4]

    # Baseline directory for specified year
    baseline_year_directory = os.path.join(BASELINE_DIRECTORY, yearstring)   

    # Comparison directory for specified year
    comparison_year_directory = os.path.join(args.satstacepath, os.path.join(os.path.join('ice', yearstring)))

    # Comparison filename
    primary_filename = os.path.join(comparison_year_directory, 'tas_ice_eustace_0_{datestring}.nc'.format(datestring=datestring))
    ancillary_filename = os.path.join(comparison_year_directory, 'ancillary_ice_eustace_0_{datestring}.nc'.format(datestring=datestring))

    # Run comparison (southern hemisphere)
    issues += comparefields(os.path.join(baseline_year_directory, 'eustace_satstace_ice_sh_{datestring}.nc'.format(datestring=datestring)),
                            primary_filename,
                            fieldmap_primary,
                            rowend=160)
    issues += comparefields(os.path.join(baseline_year_directory, 'eustace_satstace_ice_sh_{datestring}.nc'.format(datestring=datestring)),
                            ancillary_filename,
                            fieldmap_ancillary,
                            rowend=160)

    # Run comparison (northern hemisphere)
    issues += comparefields(os.path.join(baseline_year_directory, 'eustace_satstace_ice_nh_{datestring}.nc'.format(datestring=datestring)),
                            primary_filename,
                            fieldmap_primary,
                            rowstart=560)
    issues += comparefields(os.path.join(baseline_year_directory, 'eustace_satstace_ice_nh_{datestring}.nc'.format(datestring=datestring)),
                            ancillary_filename,
                            fieldmap_ancillary,
                            rowstart=560)

    # Results
    if issues:
        for issue in issues:
            print issue
    else:
        print 'OK'

if __name__=='__main__':
    main()
