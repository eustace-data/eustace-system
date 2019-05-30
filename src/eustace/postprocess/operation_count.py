import argparse
import datetime
from dateutil.relativedelta import relativedelta

def operation_dates(reference_time, operation_index):
    """Get datetime objects for a given operation index"""
    
    first_day_date = reference_time + relativedelta(months = operation_index)
    last_day_date  = first_day_date + relativedelta(months = 1) - relativedelta(days = 1)
    
    date_list = [first_day_date + relativedelta(days = dayindex) for dayindex in range( (last_day_date-first_day_date).days+1 ) ]
    
    return date_list

def operation_count(start_year, end_year):
    """Get the number of operations between two years"""
    
    first_day_date = datetime.date(start_year, 1, 1)
    last_day_date = datetime.date(end_year, 12, 31)
    
    timedelta = relativedelta(last_day_date + relativedelta(days = 1), first_day_date)
    
    n_operations = timedelta.years*12 + timedelta.months
    
    return n_operations
     
def main():

    parser = argparse.ArgumentParser(description='Counts monthly operations for a range of years')
    parser.add_argument('start_year', type=int)
    parser.add_argument('end_year', type=int)
    
    args = parser.parse_args()
    print operation_count(args.start_year, args.end_year)
    
if __name__=='__main__':

    main()
