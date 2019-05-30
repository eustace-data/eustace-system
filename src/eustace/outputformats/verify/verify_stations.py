"""Verify against station data spec."""

import argparse
import sys
from verify import Verifier
from specification_stations import SpecificationStationTemperature
from specification_stations import SpecificationStationStatus


class VerifyStationData(Verifier):
    """Verifier class for station data."""

    def check_all(self, args):
        """Verify that specified file complies with EUSTACE specification for station data."""

        # Open two kinds of file (should both exist)
        # Note file names are not cross-checked the file name pattern against
        # specification document

        pathname_temperature = (args.pathnameprefix + '_temperature.nc') if args.pathnameprefix else None
        self.check_file(pathname_temperature, SpecificationStationTemperature(), 'Open temperature file')

        pathname_status = (args.pathnameprefix + '_status.nc') if args.pathnameprefix else None
        self.check_file(pathname_status, SpecificationStationStatus(), 'Open status file')


def main():
    """Entry point for running stand-alone."""

    # read program arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('pathnameprefix')
    args = parser.parse_args()

    # do verification
    verified = VerifyStationData().run(args, sys.stdout)

    # give zero exit value if verified ok, non-zero otherwise
    sys.exit(0 if verified else 1)

if __name__ == "__main__":
    main()
