#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Package: ECCpy, for EC50 calculation in python
Author: Mark Teese
License: ECCpy is free software, distributed under the GNU Lesser General Public License 3 (LGPLv3)
"""
import thoipapy.Sine_Curve.eccpy.curvefit
import thoipapy.Sine_Curve.eccpy.analysis
import thoipapy.Sine_Curve.eccpy.judgefit
import thoipapy.Sine_Curve.eccpy.settings
import thoipapy.Sine_Curve.eccpy.tools
#import thoipapy.Sine_Curve.eccpy.curvefit.run_curvefit
from thoipapy.Sine_Curve.eccpy.analysis import run_analysis
from thoipapy.Sine_Curve.eccpy.analysis import compare_rawdata

# from os.path import dirname, basename, isfile
# import glob
# modules = glob.glob(dirname(__file__)+"/*.py")
# __all__ = [ basename(f)[:-3] for f in modules if isfile(f)]