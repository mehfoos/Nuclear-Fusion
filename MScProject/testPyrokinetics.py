#!/usr/bin/env python3
"""
Test/Play with Pyrokinetics + GS2
"""

from freegs import _geqdsk
from pyrokinetics import Pyro

def main():
    """ Main entry point of the app """
    print("hello world")

    pyro = Pyro(
        gk_code="/gs2/tests/linear_tests/cyclone_itg/cyclone_itg_low_res_base.in",
    )

    pyro.write_gk_file("input.cgyro", gk_code="CGYRO")


if __name__ == "__main__":
    """ This is executed when run from the command line """
    main()
