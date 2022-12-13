#!/usr/bin/env python

from setuptools import find_packages, setup

setup(name="wham_derivatives",
      version="0.1",
      description="A code that calculates wham derivatives",
      author="Zeke Piskulich",
      author_email='piskulichz@gmail.com',
      platforms=["any"],  # or more specific, e.g. "win32", "cygwin", "osx"
      license="BSD",
      url="https://github.com/piskuliche/wham_derivatives",
      packages=find_packages('wham_derivatives'),
      )