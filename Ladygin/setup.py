#!/usr/bin/env python

from distutils.core import setup

setup (name = 'ConvTrack',
       version = '0.0.1',
       description = 'Tracker of kpoint convergence',
       packages=['ConvergenceTracker'],
       scripts=['ConvTrack'],
       install_requires=['ase'])