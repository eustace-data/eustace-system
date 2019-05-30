"""Apply observation count to multiple sources."""

import argparse
import json
from obscount import ObsCounter

def main():

    # Parse arguments
    parser = argparse.ArgumentParser('obscount')
    parser.add_argument('--datapath', required=True, help='base path of data')
    parser.add_argument('--sourcelist', required=True, help='local filename of source list')
    args = parser.parse_args()

    # Load source list
    sourcelist = json.load(open(args.sourcelist))

    # List for results
    resultslist = { }

    # Iterate over sources
    for sourcedetails in sourcelist:

        # Results for this source
        sourceresults = { }

        # Make entry in results list
        resultslist[sourcedetails['source']] = sourceresults

        # Each source should have one or more observables
        for observable in sourcedetails['observables']:

            try:

                # Build counter object
                counter = ObsCounter(
                    args.datapath,
                    str( sourcedetails['source'] ),
                    str( observable ),
                    str( sourcedetails['startdate'] ),
                    str( sourcedetails['enddate'] ))

                # Attempt count
                sourceresults[observable] = counter.summary_statistics().__dict__

            except Exception as whatwentwrong:

                # Store error if fails
                sourceresults[observable] = { 'error': str(whatwentwrong) }
                
    # Print results
    print json.dumps(resultslist, indent=4)

if __name__ == '__main__':
    main()
