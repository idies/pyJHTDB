#! /usr/bin/env python

from setuptools import setup
import datetime

now = datetime.datetime.now()

#TODO: change format string at some point before the year 9999
date_name = '{0:0>4}{1:0>2}{2:0>2}'.format(now.year, now.month, now.day)

setup(
        name = 'pyJHTDB',
        version = date_name,
        packages = ['pyJHTDB',])

