"""Verify against global field spec."""

import sys
import argparse
from verify import Verifier
from specification_globalfield import SpecificationGlobalField


class VerifyGlobalField(Verifier):
    """Verifier class for global field files. """

    def check_all(self, args):
        """Verify that specified file complies with EUSTACE specification for global fields."""

        self.check_file(args.pathname, SpecificationGlobalField(), 'Open field file')


def main():
    """Entry point for running stand-alone."""

    # read program arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('pathname')
    args = parser.parse_args()

    # do verification
    verified = VerifyGlobalField().run(args, sys.stdout)

    # give zero exit value if verified ok, non-zero otherwise
    sys.exit(0 if verified else 1)

if __name__ == "__main__":
    main()
