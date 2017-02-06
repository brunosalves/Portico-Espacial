# -*- coding: utf-8 -*-
"""
Created on Sat Feb 04 17:50:06 2017

@author: Bruno1

Erros
"""
from __future__ import division
import numpy as np

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class CalculationError(Error):
    """Exception raised for errors found in the calculations.

    Attributes:
        expr -- input expression in which the error occurred
        msg  -- explanation of the error
    """

    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg