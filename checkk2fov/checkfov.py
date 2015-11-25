#! /usr/bin/env python

### a simple script that will call K2fov
from __future__ import division, print_function
import sys

import K2onSilicon
import definefov
import fov
import greatcircle
import projection
import rotate
import run_list
	
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run K2onSilicon to find which targets in a \
        list call on active silicon for a given K2 campaign.')
    parser.add_argument('csv_file',
        help='Name of input csv file with targets, column are \
        Ra_degrees, Dec_degrees, Kepmag', type=str)
    parser.add_argument('campaign',
        help='K2 Campaign number', type=int)
    parser.add_argument('epicfile',
        help='EPIC catalog file', type=str)
    args = parser.parse_args()

    K2onSilicon.K2onSilicon(args.csv_file,args.campaign,args.epicfile)


