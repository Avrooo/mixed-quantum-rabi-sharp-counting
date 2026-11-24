"""Recomputes Table (exact states vs Fock diagonalisation): level index, residual ||(H-E)psi||, overlap, unitary-image check."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from math import factorial, sqrt
from scipy.linalg import eigh
import exact as ex
from make_figures2 import ref, fockH, Eexact, tf, g2, Om, Dev  # reuse refined states (re-runs figure script)
