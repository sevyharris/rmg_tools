#!/bin/bash

# Script to run select autopep8 fixes (less aggressive than default)
autopep8 -i --aggressive --select=E201,E202,E221,E225,E226,E227,E228,E231,E242,E253,E261,E262,E265,W291,W292,W293,W504 $1

# usage:
# >>> select_autopep8.sh my_file.py

