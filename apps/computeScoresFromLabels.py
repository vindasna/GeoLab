#!/usr/bin/python3

import sys, os, subprocess, shutil

import numpy as np
import scipy.stats
import scipy.special
import matplotlib.pyplot as plt

import argparse
import textwrap
from argparse import RawTextHelpFormatter

import multiprocessing



#---- For printing colors in terminal ----#
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#-----------------------------------------#

DOC = """
---------------------------

Compute scores for a segmentation

Command example:
----------------
python3 computeScoresFromLabels.py \
-pl predicted.labels \
-pd predicted.dict \
-tl predicted.labels \
-td predicted.dict \
-v 1 \

"""
