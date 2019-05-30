"""Extreme and areal value QC checks on EUSTACE analysis products"""

import numpy
import iris
from dateutil import parser
import datetime
from dateutil.relativedelta import relativedelta
import os.path

import argparse
import flags


"""

Percentile calculation

"""

def cell_percentile_calculation(input_directory, output_directory, iteration, start_year, end_year, calendar_day):
    """location based threshold calculation using a fixed period on given calendar day"""
    
    refererence_leap_year = 2000
    year_shift_start = start_year - refererence_leap_year
    year_shift_end   = end_year - refererence_leap_year
    reference_time = datetime.datetime(refererence_leap_year, 1, 1) +relativedelta(days=calendar_day)
    
    dates_to_load = [reference_time +relativedelta(years = yearshift) for yearshift in range(year_shift_start, year_shift_end+1)]
    print dates_to_load
    
    files_to_load = [os.path.join(input_directory, 'eustace_analysis_'+str(iteration), '{:04d}'.format(processdate.year), 'eustace_analysis_{:d}_{:04d}{:02d}{:02d}.nc'.format(iteration, processdate.year, processdate.month, processdate.day)) for processdate in dates_to_load]
    
    outfile = os.path.join(output_directory, 'eustace_analysis_{:d}_percentiles_{:02d}{:02d}'.format(iteration, reference_time.month, reference_time.day))
    
    cubes = iris.load(files_to_load, u'air_temperature')
    #print cubes
    #a = [cube.data for cube in cubes]
    analysis_data = numpy.stack([numpy.squeeze( cube.data ) for cube in cubes], axis=2)
    lower_quartile_block, median_block, upper_quartile_block = numpy.percentile(analysis_data, [25, 50, 75], axis=2)
    
    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))
    
    numpy.savez(outfile, lower_quartile_block=lower_quartile_block, median_block=median_block, upper_quartile_block=upper_quartile_block)


def cell_percentile_calculation_month(input_directory, output_directory, iteration, start_year, end_year, month_number):
    """location based threshold calculation for all days in month"""
    
    refererence_leap_year = 2000
    reference_date = datetime.datetime(refererence_leap_year, 1, 1)
    
    start_date = datetime.datetime(refererence_leap_year, month_number, 1)
    end_date = datetime.datetime(refererence_leap_year, month_number+1, 1)-relativedelta(days=1) if month_number != 12 else datetime.datetime(refererence_leap_year+1, 1, 1)-relativedelta(days=1)
    
    for calendar_day in range( (start_date - reference_date).days, (end_date - reference_date).days + 1):
        print "calendar_day", calendar_day
        cell_percentile_calculation(input_directory, output_directory, iteration, start_year, end_year, calendar_day)

def main():

    print 'Submission of advanced standard analysis jobs'
    
    parser = argparse.ArgumentParser(description='Monthly batches of masking operations')
    
    parser.add_argument('--month_number',  type=int, default=1, help='the number of the month to compute')
    parser.add_argument('--input_directory', type=str, default="/work/scratch/cmorice/advanced_standard/", help='directory in which masking information can be found')
    parser.add_argument('--iteration',  type=int, default=9, help='operation index at which the analysis grid produced')
    parser.add_argument('--output_directory',  type=str, default="/gws/nopw/j04/eustace_vol2/masking/", help='operation index at which the analysis grid produced')
    parser.add_argument('--start_year',  type=int, default=1880, help='first analysis year to load for constraint tests')
    parser.add_argument('--end_year',  type=int, default=2015, help='last analysis year to load for constraint tests')
    
    args = parser.parse_args()

    cell_percentile_calculation_month(args.input_directory, args.output_directory, args.iteration, args.start_year, args.end_year, args.month_number)

if __name__ == '__main__':

    main()
