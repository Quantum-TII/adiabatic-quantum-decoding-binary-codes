import numpy as np
from sympy import symbols, expand, Matrix, Identity
from sympy.matrices.dense import matrix2numpy
from qibo import matrices, hamiltonians, models, callbacks
import matplotlib.pyplot as plt
import collections


def num2list(n, b):
    """Turn a number into a list of its binary bits.
    Args:
        n (int): number to turn into binary.
        b (int): number of bits.
        
    Returns:
        x (list): list of binary bits that represent number n.
    """
    n = "{0:0{bits}b}".format(n, bits=b)
    x = []
    for i in range(b):
        x.append(int(n[i]))
    return x


def f_3(x1, x2, x3):
    """Function for the 3-7 nonlinear code.
    
    """
    f = np.array([x1*x2*x3 - x1*x2 - x1*x3 + x1 - x2*x3 + x2, 
                  x1*x2 - x1*x3 + x3, 
                  x1*x2*x3 - x1*x2 - 2*x2*x3 + x2 + x3,  
                  x1*x2*x3 - x1*x3 - x2*x3 + 1, 
                  x1*x2*x3 + x1*x3 - x1 - x2*x3 + 1, 
                  2*x1*x2*x3 - 2*x1*x2 + x1 - x2*x3 + x2, 
                  -2*x1*x2*x3 + x1*x3 + 2*x2*x3 - x2 - x3 + 1])
    return f


def f_8(x1, x2, x3, x4, x5, x6, x7, x8):
    """Function for the 8-16 nonlinear code.
    
    """
    f = np.array([-2*x1*x2 + x1 + x2,
                  -2*x2*x3 + x2 + x3,
                  -2*x3*x4 + x3 + x4,
                  -2*x4*x5 + x4 + x5,
                  2*x5*x6*x7*x8 - x5 - x6*x7*x8 + 1,
                  -8*x2*x4*x5*x6*x7*x8 + 8*x2*x4*x5*x6*x8 + 16*x2*x4*x5*x7*x8 - 8*x2*x4*x5*x7 - 8*x2*x4*x5*x8 + 4*x2*x4*x5 + 4*x2*x4*x6*x7*x8 - 4*x2*x4*x6*x8 - 8*x2*x4*x7*x8 + 4*x2*x4*x7 + 4*x2*x4*x8 - 2*x2*x4 + 4*x2*x5*x6*x7*x8 - 4*x2*x5*x6*x8 - 8*x2*x5*x7*x8 + 4*x2*x5*x7 + 4*x2*x5*x8 - 2*x2*x5 - 2*x2*x6*x7*x8 + 2*x2*x6*x8 +4*x2*x7*x8 - 2*x2*x7 - 2*x2*x8 + x2 + 4*x4*x5*x6*x7*x8 - 4*x4*x5*x6*x8 - 8*x4*x5*x7*x8 + 4*x4*x5*x7 + 4*x4*x5*x8 - 2*x4*x5 - 2*x4*x6*x7*x8 + 2*x4*x6*x8 +4*x4*x7*x8 - 2*x4*x7 - 2*x4*x8 + x4 - 2*x5*x6*x7*x8 + 2*x5*x6*x8 + 4*x5*x7*x8 - 2*x5*x7 - 2*x5*x8 + x5 + x6*x7*x8 - x6*x8 - 2*x7*x8 + x7 + x8,
                  8*x1*x2*x4*x6*x7*x8 - 8*x1*x2*x4*x6*x8 - 8*x1*x2*x4*x7*x8 + 8*x1*x2*x4*x7 + 8*x1*x2*x4*x8 - 4*x1*x2*x4 - 4*x1*x2*x6*x7*x8 + 4*x1*x2*x6*x8 + 4*x1*x2*x7*x8 - 4*x1*x2*x7 - 4*x1*x2*x8 + 2*x1*x2 - 4*x1*x4*x6*x7*x8 + 4*x1*x4*x6*x8 + 4*x1*x4*x7*x8 - 4*x1*x4*x7 - 4*x1*x4*x8 + 2*x1*x4 + 2*x1*x6*x7*x8 - 2*x1*x6*x8 -2*x1*x7*x8 + 2*x1*x7 + 2*x1*x8 - x1 - 4*x2*x4*x6*x7*x8 + 4*x2*x4*x6*x8 + 4*x2*x4*x7*x8 - 4*x2*x4*x7 - 4*x2*x4*x8 + 2*x2*x4 + 2*x2*x6*x7*x8 - 2*x2*x6*x8 -2*x2*x7*x8 + 2*x2*x7 + 2*x2*x8 - x2 + 2*x4*x6*x7*x8 - 2*x4*x6*x8 - 2*x4*x7*x8 + 2*x4*x7 + 2*x4*x8 - x4 - x6*x7*x8 + x6*x8 + x7*x8 - x7 - x8 + 1,
                  96*x1*x2*x3*x4*x5*x6*x7*x8 - 32*x1*x2*x3*x4*x5*x6*x7 - 32*x1*x2*x3*x4*x5*x6*x8 -32*x1*x2*x3*x4*x5*x7*x8 + 16*x1*x2*x3*x4*x5 - 48*x1*x2*x3*x4*x6*x7*x8 + 16*x1*x2*x3*x4*x6*x7 + 16*x1*x2*x3*x4*x6*x8 + 16*x1*x2*x3*x4*x7*x8 - 8*x1*x2*x3*x4 - 48*x1*x2*x3*x5*x6*x7*x8 + 16*x1*x2*x3*x5*x6*x7 + 16*x1*x2*x3*x5*x6*x8 + 16*x1*x2*x3*x5*x7*x8 - 8*x1*x2*x3*x5 + 24*x1*x2*x3*x6*x7*x8 - 8*x1*x2*x3*x6*x7 - 8*x1*x2*x3*x6*x8 - 8*x1*x2*x3*x7*x8 + 4*x1*x2*x3 - 48*x1*x2*x4*x5*x6*x7*x8 + 16*x1*x2*x4*x5*x6*x7 + 16*x1*x2*x4*x5*x6*x8 + 16*x1*x2*x4*x5*x7*x8 - 8*x1*x2*x4*x5 + 24*x1*x2*x4*x6*x7*x8 - 8*x1*x2*x4*x6*x7 - 8*x1*x2*x4*x6*x8 - 8*x1*x2*x4*x7*x8 + 4*x1*x2*x4 + 24*x1*x2*x5*x6*x7*x8 - 8*x1*x2*x5*x6*x7 - 8*x1*x2*x5*x6*x8 - 8*x1*x2*x5*x7*x8 + 4*x1*x2*x5 - 12*x1*x2*x6*x7*x8 + 4*x1*x2*x6*x7 + 4*x1*x2*x6*x8 + 4*x1*x2*x7*x8 - 2*x1*x2 - 48*x1*x3*x4*x5*x6*x7*x8 + 16*x1*x3*x4*x5*x6*x7 + 16*x1*x3*x4*x5*x6*x8 + 16*x1*x3*x4*x5*x7*x8 - 8*x1*x3*x4*x5 + 24*x1*x3*x4*x6*x7*x8 - 8*x1*x3*x4*x6*x7 - 8*x1*x3*x4*x6*x8 - 8*x1*x3*x4*x7*x8 + 4*x1*x3*x4 + 24*x1*x3*x5*x6*x7*x8 - 8*x1*x3*x5*x6*x7 - 8*x1*x3*x5*x6*x8 - 8*x1*x3*x5*x7*x8 + 4*x1*x3*x5 - 12*x1*x3*x6*x7*x8 + 4*x1*x3*x6*x7 + 4*x1*x3*x6*x8 + 4*x1*x3*x7*x8 - 2*x1*x3 + 24*x1*x4*x5*x6*x7*x8 -8*x1*x4*x5*x6*x7 - 8*x1*x4*x5*x6*x8 - 8*x1*x4*x5*x7*x8 + 4*x1*x4*x5 - 12*x1*x4*x6*x7*x8 + 4*x1*x4*x6*x7 + 4*x1*x4*x6*x8 + 4*x1*x4*x7*x8 - 2*x1*x4 - 12*x1*x5*x6*x7*x8 + 4*x1*x5*x6*x7 + 4*x1*x5*x6*x8 + 4*x1*x5*x7*x8 - 2*x1*x5 + 6*x1*x6*x7*x8 - 2*x1*x6*x7 - 2*x1*x6*x8 - 2*x1*x7*x8 + x1 - 48*x2*x3*x4*x5*x6*x7*x8 + 16*x2*x3*x4*x5*x6*x7 + 16*x2*x3*x4*x5*x6*x8 + 16*x2*x3*x4*x5*x7*x8 - 8*x2*x3*x4*x5 + 24*x2*x3*x4*x6*x7*x8 - 8*x2*x3*x4*x6*x7 -8*x2*x3*x4*x6*x8 - 8*x2*x3*x4*x7*x8 + 4*x2*x3*x4 + 24*x2*x3*x5*x6*x7*x8 - 8*x2*x3*x5*x6*x7 - 8*x2*x3*x5*x6*x8 - 8*x2*x3*x5*x7*x8 + 4*x2*x3*x5 - 12*x2*x3*x6*x7*x8 + 4*x2*x3*x6*x7 + 4*x2*x3*x6*x8 + 4*x2*x3*x7*x8 - 2*x2*x3 + 24*x2*x4*x5*x6*x7*x8 - 8*x2*x4*x5*x6*x7 - 8*x2*x4*x5*x6*x8 - 8*x2*x4*x5*x7*x8 + 4*x2*x4*x5 - 12*x2*x4*x6*x7*x8 + 4*x2*x4*x6*x7 + 4*x2*x4*x6*x8 + 4*x2*x4*x7*x8 -2*x2*x4 - 12*x2*x5*x6*x7*x8 + 4*x2*x5*x6*x7 + 4*x2*x5*x6*x8 + 4*x2*x5*x7*x8 - 2*x2*x5 + 6*x2*x6*x7*x8 - 2*x2*x6*x7 - 2*x2*x6*x8 - 2*x2*x7*x8 + x2 + 24*x3*x4*x5*x6*x7*x8 - 8*x3*x4*x5*x6*x7 - 8*x3*x4*x5*x6*x8 - 8*x3*x4*x5*x7*x8 + 4*x3*x4*x5 - 12*x3*x4*x6*x7*x8 + 4*x3*x4*x6*x7 + 4*x3*x4*x6*x8 + 4*x3*x4*x7*x8 -2*x3*x4 - 12*x3*x5*x6*x7*x8 + 4*x3*x5*x6*x7 + 4*x3*x5*x6*x8 + 4*x3*x5*x7*x8 - 2*x3*x5 + 6*x3*x6*x7*x8 - 2*x3*x6*x7 - 2*x3*x6*x8 - 2*x3*x7*x8 + x3 - 12*x4*x5*x6*x7*x8 + 4*x4*x5*x6*x7 + 4*x4*x5*x6*x8 + 4*x4*x5*x7*x8 - 2*x4*x5 + 6*x4*x6*x7*x8 - 2*x4*x6*x7 - 2*x4*x6*x8 - 2*x4*x7*x8 + x4 + 6*x5*x6*x7*x8 - 2*x5*x6*x7 - 2*x5*x6*x8 - 2*x5*x7*x8 + x5 - 3*x6*x7*x8 + x6*x7 + x6*x8 + x7*x8,
                  -8*x2*x5*x6*x7*x8 + 4*x2*x5*x6*x7 + 4*x2*x5*x8 - 2*x2*x5 + 4*x2*x6*x7*x8 - 2*x2*x6*x7 - 2*x2*x8 + x2 + 4*x5*x6*x7*x8 - 2*x5*x6*x7 - 2*x5*x8 + x5 - 2*x6*x7*x8 + x6*x7 + x8,
                  -2*x1*x6*x7*x8 + 4*x1*x6*x7 + 2*x1*x6*x8 - 2*x1*x6 + 2*x1*x7*x8 - 2*x1*x7 - 2*x1*x8 + x1 + x6*x7*x8 - 2*x6*x7 - x6*x8 + x6 - x7*x8 + x7 + x8,
                  -4*x1*x4*x6*x7 - 4*x1*x4*x6*x8 + 4*x1*x4*x6 + 4*x1*x4*x7 - 2*x1*x4 + 2*x1*x6*x7 + 2*x1*x6*x8 - 2*x1*x6 - 2*x1*x7 + x1 + 2*x4*x6*x7 + 2*x4*x6*x8 - 2*x4*x6 - 2*x4*x7 + x4 - x6*x7 - x6*x8 + x6 + x7,
                  8*x1*x3*x5*x6*x7*x8 - 8*x1*x3*x5*x6*x8 + 8*x1*x3*x5*x6 - 8*x1*x3*x5*x7*x8 + 8*x1*x3*x5*x8 - 4*x1*x3*x5 - 4*x1*x3*x6*x7*x8 + 4*x1*x3*x6*x8 - 4*x1*x3*x6 + 4*x1*x3*x7*x8 - 4*x1*x3*x8 + 2*x1*x3 - 4*x1*x5*x6*x7*x8 + 4*x1*x5*x6*x8 - 4*x1*x5*x6 + 4*x1*x5*x7*x8 - 4*x1*x5*x8 + 2*x1*x5 + 2*x1*x6*x7*x8 - 2*x1*x6*x8 +2*x1*x6 - 2*x1*x7*x8 + 2*x1*x8 - x1 - 4*x3*x5*x6*x7*x8 + 4*x3*x5*x6*x8 - 4*x3*x5*x6 + 4*x3*x5*x7*x8 - 4*x3*x5*x8 + 2*x3*x5 + 2*x3*x6*x7*x8 - 2*x3*x6*x8 +2*x3*x6 - 2*x3*x7*x8 + 2*x3*x8 - x3 + 2*x5*x6*x7*x8 - 2*x5*x6*x8 + 2*x5*x6 - 2*x5*x7*x8 + 2*x5*x8 - x5 - x6*x7*x8 + x6*x8 - x6 + x7*x8 - x8 + 1,
                  8*x2*x3*x4*x6*x7*x8 + 8*x2*x3*x4*x6*x8 - 8*x2*x3*x4*x6 - 8*x2*x3*x4*x7*x8 + 4*x2*x3*x4 - 4*x2*x3*x6*x7*x8 - 4*x2*x3*x6*x8 + 4*x2*x3*x6 + 4*x2*x3*x7*x8 - 2*x2*x3 - 4*x2*x4*x6*x7*x8 - 4*x2*x4*x6*x8 + 4*x2*x4*x6 + 4*x2*x4*x7*x8 - 2*x2*x4 + 2*x2*x6*x7*x8 + 2*x2*x6*x8 - 2*x2*x6 - 2*x2*x7*x8 + x2 - 4*x3*x4*x6*x7*x8 - 4*x3*x4*x6*x8 + 4*x3*x4*x6 + 4*x3*x4*x7*x8 - 2*x3*x4 + 2*x3*x6*x7*x8 + 2*x3*x6*x8 - 2*x3*x6 - 2*x3*x7*x8 + x3 + 2*x4*x6*x7*x8 + 2*x4*x6*x8 - 2*x4*x6 - 2*x4*x7*x8 + x4 - x6*x7*x8 - x6*x8 + x6 + x7*x8,
                  -16*x1*x3*x4*x5*x6*x8 + 16*x1*x3*x4*x5*x7*x8 - 16*x1*x3*x4*x5*x7 + 8*x1*x3*x4*x5+ 8*x1*x3*x4*x6*x8 - 8*x1*x3*x4*x7*x8 + 8*x1*x3*x4*x7 - 4*x1*x3*x4 + 8*x1*x3*x5*x6*x8 - 8*x1*x3*x5*x7*x8 + 8*x1*x3*x5*x7 - 4*x1*x3*x5 - 4*x1*x3*x6*x8+ 4*x1*x3*x7*x8 - 4*x1*x3*x7 + 2*x1*x3 + 8*x1*x4*x5*x6*x8 - 8*x1*x4*x5*x7*x8 + 8*x1*x4*x5*x7 - 4*x1*x4*x5 - 4*x1*x4*x6*x8 + 4*x1*x4*x7*x8 - 4*x1*x4*x7 + 2*x1*x4 - 4*x1*x5*x6*x8 + 4*x1*x5*x7*x8 - 4*x1*x5*x7 + 2*x1*x5 + 2*x1*x6*x8 - 2*x1*x7*x8 + 2*x1*x7 - x1 + 8*x3*x4*x5*x6*x8 - 8*x3*x4*x5*x7*x8 + 8*x3*x4*x5*x7 - 4*x3*x4*x5 - 4*x3*x4*x6*x8 + 4*x3*x4*x7*x8 - 4*x3*x4*x7 + 2*x3*x4 - 4*x3*x5*x6*x8 + 4*x3*x5*x7*x8 - 4*x3*x5*x7 + 2*x3*x5 + 2*x3*x6*x8 - 2*x3*x7*x8 +2*x3*x7 - x3 - 4*x4*x5*x6*x8 + 4*x4*x5*x7*x8 - 4*x4*x5*x7 + 2*x4*x5 + 2*x4*x6*x8- 2*x4*x7*x8 + 2*x4*x7 - x4 + 2*x5*x6*x8 - 2*x5*x7*x8 + 2*x5*x7 - x5 - x6*x8 + x7*x8 - x7 + 1,
                  -32*x1*x2*x3*x5*x6*x7*x8 + 32*x1*x2*x3*x5*x6*x8 - 16*x1*x2*x3*x5*x6 + 16*x1*x2*x3*x5*x7*x8 - 16*x1*x2*x3*x5*x8 + 8*x1*x2*x3*x5 + 16*x1*x2*x3*x6*x7*x8 - 16*x1*x2*x3*x6*x8 + 8*x1*x2*x3*x6 - 8*x1*x2*x3*x7*x8 + 8*x1*x2*x3*x8 - 4*x1*x2*x3 + 16*x1*x2*x5*x6*x7*x8 - 16*x1*x2*x5*x6*x8 + 8*x1*x2*x5*x6 - 8*x1*x2*x5*x7*x8 + 8*x1*x2*x5*x8 - 4*x1*x2*x5 - 8*x1*x2*x6*x7*x8 + 8*x1*x2*x6*x8 - 4*x1*x2*x6 + 4*x1*x2*x7*x8 - 4*x1*x2*x8 + 2*x1*x2 + 16*x1*x3*x5*x6*x7*x8 - 16*x1*x3*x5*x6*x8 + 8*x1*x3*x5*x6 - 8*x1*x3*x5*x7*x8 + 8*x1*x3*x5*x8 - 4*x1*x3*x5 - 8*x1*x3*x6*x7*x8 + 8*x1*x3*x6*x8 - 4*x1*x3*x6 + 4*x1*x3*x7*x8 - 4*x1*x3*x8 + 2*x1*x3 - 8*x1*x5*x6*x7*x8 + 8*x1*x5*x6*x8 - 4*x1*x5*x6 + 4*x1*x5*x7*x8 - 4*x1*x5*x8 + 2*x1*x5 + 4*x1*x6*x7*x8 - 4*x1*x6*x8 + 2*x1*x6 - 2*x1*x7*x8 + 2*x1*x8 - x1 + 16*x2*x3*x5*x6*x7*x8 - 16*x2*x3*x5*x6*x8 + 8*x2*x3*x5*x6 - 8*x2*x3*x5*x7*x8 + 8*x2*x3*x5*x8 - 4*x2*x3*x5 - 8*x2*x3*x6*x7*x8+ 8*x2*x3*x6*x8 - 4*x2*x3*x6 + 4*x2*x3*x7*x8 - 4*x2*x3*x8 + 2*x2*x3 - 8*x2*x5*x6*x7*x8 + 8*x2*x5*x6*x8 - 4*x2*x5*x6 + 4*x2*x5*x7*x8 - 4*x2*x5*x8 + 2*x2*x5 + 4*x2*x6*x7*x8 - 4*x2*x6*x8 + 2*x2*x6 - 2*x2*x7*x8 + 2*x2*x8 - x2 - 8*x3*x5*x6*x7*x8 + 8*x3*x5*x6*x8 - 4*x3*x5*x6 + 4*x3*x5*x7*x8 - 4*x3*x5*x8 + 2*x3*x5 + 4*x3*x6*x7*x8 - 4*x3*x6*x8 + 2*x3*x6 - 2*x3*x7*x8 + 2*x3*x8 - x3 + 4*x5*x6*x7*x8 - 4*x5*x6*x8 + 2*x5*x6 - 2*x5*x7*x8 + 2*x5*x8 - x5 - 2*x6*x7*x8 + 2*x6*x8 - x6 + x7*x8 - x8 + 1,
                  2*x3*x6*x7*x8 - 2*x3*x6*x7 + 2*x3*x6 - 2*x3*x7*x8 + 2*x3*x7 - x3 - x6*x7*x8 + x6*x7 - x6 + x7*x8 - x7 + 1])
    return f


def r_3():
    """Target state for the 3-7 nonlinear code.
    
    """
    return np.array([0, 1, 0, 1, 0, 0, 0])


def r_8():
    """Target state for the 8-16 nonlinear code.
    
    """
    return np.array([1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0])


def to_ancilla(x1, x2, x3, a):
    """How a 3-body interaction between x1, x2 and x3 is decomposed using an ancilla a.
    
    """
    return 2*x2*a - 2*x2 + 2*x3*a - 2*x3 - x1*a + x1 -3*a + 3 + x2*x3


def substitutions(sym):
    """Substitutions needed as 0 or 1 squared are equal to themselves.
    
    """
    s = []
    for x in sym:
        s.append((x**2, x))
    return s


def sym_dict_num(sym):
    """Dictionary with each symbol and their corresponding qubit.
    """
    sym_num = {}
    for i, s in enumerate(sym):
        sym_num[s] = i
    return sym_num


def to_hamiltonian(f, r, sym):
    """Tranform the non-linear functions into a hamiltonian that minimizes to find r and simplifies it.
    
    """
    h = np.sum((f-r)**2)
    h = expand(h)
    h = h.subs(substitutions(sym))
    return h


def to_two_body(h, ancillas, bits):
    """Turns 3-body hamiltonian into a 2-body one by use of an ancilla.
    
    """
    anc = 0
    for i in reversed(range(3, bits+1)):
        print()
        print(f'Dealing with {i}-term interactions.')
        terms_s, terms_c, overall_constant = symbolic_to_data(h)
        print(f'Total terms in the hamiltonian: {len(terms_s)}')
        for term in terms_s:
            if len(term) == i:
                xtemp1 = 1
                #for j in range(i-2):
                for j in range(int(i/3)):
                #for j in range(int(i/2)):
                    xtemp1 *= term[j]
                xtemp2 = 1
                #for j in range(i-2, i-1):
                for j in range(int(i/3), 2*int(i/3)):
                #for j in range(int(i/2), i-1):
                    xtemp2 *= term[j]
                xtemp3 = 1
                #for j in range(i-1, i):
                for j in range(2*int(i/3), i):
                #for j in range(i-1, i):
                    xtemp3 *= term[j]
                h = expand(h.subs(xtemp1*xtemp2*xtemp3, to_ancilla(xtemp1, xtemp2, xtemp3, ancillas[anc])))
                anc += 1
                print(f'Added ancilla number {anc}')
    print('Total ancillas used: {}\n'.format(anc))
    return h, anc


def to_pauli(h, sym, sym_z):
    """Turn the hamiltonian into pauli matrices.
    
    """
    s = []
    for i in range(len(sym)):
        s.append((sym[i], (1-sym_z[i])/2))
    hz = h.subs(s)
    hz = expand(hz)
    return hz


def m(qubits, n1, n2, m1, m2):
    """Create the full matrix of a 2-body interaction between qubits n1 and n2.
    
    """
    m = 1
    for i in range(qubits):
        if i == n1:
            m = np.kron(m, m1)
        elif i == n2:
            m = np.kron(m, m2)
        else:
            m = np.kron(m, matrices.I)
    return m


def sub_matrix(sym_z):
    """Create the substitution scheme for the matrices. Important to substitute the 2-body terms first.
    
    """
    s = []
    n = len(sym_z)
    for i in range(n):
        for j in range(i+1, n):
            s.append((sym_z[i]*sym_z[j], Matrix(m(n, i, j, matrices.Z, matrices.Z))))
    for i in range(n):
        s.append((sym_z[i], Matrix(m(n, i, i, matrices.Z, matrices.Z))))
    return s


def to_matrix(h, sym_z):
    """Turn the symbol hamiltonian into a numpy array.
    
    """
    terms_s, terms_c, overall_constant = symbolic_to_data(h)
    n = len(sym_z)
    hm = h.subs(sub_matrix(sym_z))
    hm = hm.subs(overall_constant, overall_constant*Identity(2**n).as_mutable())
    hm = hm.evalf()
    hm = matrix2numpy(hm, dtype=np.complex128)
    return hm


def h0(qubits):
    """Initial hamiltonian for adiabatic evolution.
    Args:
        qubits (int): # of total qubits in the instance.

    Return:
        h0 (np.array): 2**qubits x 2**qubits initial Hamiltonian.
    """
    h0 = 0
    for i in range(qubits):
        h0 += 0.5*(np.eye(2**qubits)-m(qubits, i, i, matrices.X, matrices.X))
    return h0


def plot(qubits, ground, first, gap, dt, T):
    """Get the first two eigenvalues and the gap energy
    Args:
        qubits (int): # of total qubits in the instance.
        ground (list): ground state energy during the evolution.
        first (list): first excited state during the evolution.
        gap (list): gap energy during the evolution.
        T (float): Final time for the schedue.
        dt (float): time interval for the evolution.

    Returns:
        {}_qubits_energy.png: energy evolution of the ground and first excited state.
        {}_qubits_gap_energy.png: gap evolution during the adiabatic process.
    """
    fig, ax = plt.subplots()
    times = np.arange(0, T+dt, dt)
    ax.plot(times, ground, label='ground state', color='C0')
    ax.plot(times, first, label='first excited state', color='C1')
    plt.ylabel('energy')
    plt.xlabel('schedule')
    plt.title('Energy during adiabatic evolution')
    ax.legend()
    fig.tight_layout()
    plt.show()
    fig.savefig('{}_bits_energy.png'.format(qubits), dpi=300, bbox_inches='tight')
    fig, ax = plt.subplots()
    ax.plot(times, gap, label='gap energy', color='C0')
    plt.ylabel('energy')
    plt.xlabel('schedule')
    plt.title('Gap energy during adiabatic evolution')
    ax.legend()
    fig.tight_layout()
    plt.show()
    fig.savefig('{}_bits_gap.png'.format(qubits), dpi=300, bbox_inches='tight')
    
    
def symbolic_to_data(symbolic_hamiltonian):
    """Transforms a symbolic Hamiltonian to lists of every term.
        
    Args:
        symbolic_hamiltonian: The full Hamiltonian written with symbols.
    
    Returns:
        A dictionary that maps pairs of targets to the corresponding 4x4
        matrix that acts on this pair in the given ``symbolic_hamiltonian``.
    """ 
    terms_s = []
    terms_c = []
    overall_constant = 0
    for term in symbolic_hamiltonian.args:
        if not term.args:
            expression = (term,)
        else:
            expression = term.args

        symbols = [x for x in expression if x.is_symbol]
        numbers = [x for x in expression if not x.is_symbol]

        if len(numbers) > 1:
            raise ValueError("Hamiltonian must be expanded before using this method.")
        elif numbers:
            constant = float(numbers[0])
        else:
            constant = 1

        if not symbols:
            overall_constant += constant
            
        terms_s.append(symbols)
        terms_c.append(constant)
    
    return terms_s, terms_c, overall_constant
    
    
def symbolic_to_dwave(symbolic_hamiltonian, symbol_mnum):
    """Transforms a symbolic Hamiltonian to a dictionary of targets and matrices.
    
    Works for Hamiltonians with one and two qubit terms only.
    
    Args:
        symbolic_hamiltonian: The full Hamiltonian written with symbols.
        symbol_map: Dictionary that maps each symbol that appears in the 
            Hamiltonian to a pair of (target, matrix).
            For example {z1: (0, matrices.Z)}.
    
    Returns:
        A dictionary that maps pairs of targets to the corresponding 4x4
        matrix that acts on this pair in the given ``symbolic_hamiltonian``.
    """ 
    Q = {}
    overall_constant = 0
    for term in symbolic_hamiltonian.args:
        if not term.args:
            expression = (term,)
        else:
            expression = term.args

        symbols = [x for x in expression if x.is_symbol]
        numbers = [x for x in expression if not x.is_symbol]

        if len(numbers) > 1:
            raise ValueError("Hamiltonian must be expanded before using this method.")
        elif numbers:
            constant = float(numbers[0])
        else:
            constant = 1

        if not symbols:
            overall_constant += constant

        elif len(symbols) == 1:
            target = symbol_num[symbols[0]]
            Q[(target, target)] = constant

        elif len(symbols) == 2:
            target1 = symbol_num[symbols[0]]
            target2 = symbol_num[symbols[1]]
            Q[(target1, target2)] = constant

        else:
            raise ValueError("Only one and two qubit terms are allowed.")
    
    return Q, overall_constant
    
    
def symbolic_to_dict(symbolic_hamiltonian, symbol_map):
    """Transforms a symbolic Hamiltonian to a dictionary of targets and matrices.
    
    Works for Hamiltonians with one and two qubit terms only.
    The two qubit terms should be sufficiently many so that every 
    qubit appears as the first target at least once.
    
    Args:
        symbolic_hamiltonian: The full Hamiltonian written with symbols.
        symbol_map: Dictionary that maps each symbol that appears in the 
            Hamiltonian to a pair of (target, matrix).
            For example {z1: (0, matrices.Z)}.
    
    Returns:
        A dictionary that maps pairs of targets to the corresponding 4x4
        matrix that acts on this pair in the given ``symbolic_hamiltonian``.
    """ 
    one_qubit_terms = dict()
    two_qubit_terms = dict()
    first_targets = dict()
    overall_constant = 0
    for term in symbolic_hamiltonian.args:
        if not term.args:
            expression = (term,)
        else:
            expression = term.args

        symbols = [x for x in expression if x.is_symbol]
        numbers = [x for x in expression if not x.is_symbol]

        if len(numbers) > 1:
            raise ValueError("Hamiltonian must be expanded before using this method.")
        elif numbers:
            constant = float(numbers[0])
        else:
            constant = 1

        if not symbols:
            overall_constant += constant

        elif len(symbols) == 1:
            target, matrix = symbol_map[symbols[0]]
            one_qubit_terms[target] = constant * matrix

        elif len(symbols) == 2:
            target1, matrix1 = symbol_map[symbols[0]]
            target2, matrix2 = symbol_map[symbols[1]]
            if target1 in first_targets and target2 not in first_targets:
                two_qubit_terms[(target2, target1)] = constant * np.kron(matrix2, matrix1)
                first_targets[target2] = (target2, target1)
            else:
                two_qubit_terms[(target1, target2)] = constant * np.kron(matrix1, matrix2)
                first_targets[target1] = (target1, target2)

        else:
            raise ValueError("Only one and two qubit terms are allowed.")

    all_terms = dict(two_qubit_terms)
    for target in one_qubit_terms.keys():
        if target not in first_targets:
            raise ValueError(f"Qubit {target} has not been used as the first target.")
        pair = first_targets[target]
        all_terms[pair] = np.kron(one_qubit_terms[target], matrices.I) + two_qubit_terms[pair]
    
    return all_terms, overall_constant


def split_keys(full_dict):
    """Splits a dictionary of terms to multiple dictionaries.
    
    Each qubit should not appear in more that one term in each 
    dictionary to ensure commutation relations in the definition
    of `TrotterHamiltonian`.
    """
    all_pairs = set(full_dict.keys())
    group_pairs = [set()]
    group_singles = [set()]
    for pair in all_pairs:
        q0, q1 = pair
        flag = True
        for g, s in zip(group_pairs, group_singles):
            if q0 not in s and q1 not in s:
                s.add(q0)
                s.add(q1)
                g.add(pair)
                flag = False
                break
        if flag:
            group_pairs.append({pair})
            group_singles.append({q0, q1})
    return [{k: full_dict[k] for k in g}
            for g in group_pairs]
            

def symbolic_to_trotter(symbolic_hamiltonian, symbol_map):
    """Transforms a symbolic Hamiltonian to a Trotter Hamiltonian.
    
    Works for Hamiltonians with one and two qubit terms only.
    The two qubit terms should be sufficiently many so that every 
    qubit appears as the first target at least once.
    
    Args:
        symbolic_hamiltonian: The full Hamiltonian written with symbols.
        symbol_map: Dictionary that maps each symbol that appears in the 
            Hamiltonian to a pair of (target, matrix).
            For example {z1: (0, matrices.Z)}.
    
    Returns:
        A `TrotterHamiltonian` that implements the given `symbolic_hamiltonian`.
    """ 
    all_terms, constant = symbolic_to_dict(symbolic_hamiltonian, symbol_map)
    ham_terms = {k: hamiltonians.Hamiltonian(2, v, numpy=True) for k, v in all_terms.items()}
    groups = split_keys(ham_terms)
    return hamiltonians.TrotterHamiltonian(*groups) + constant


def h0t():
    """Generate the 2 qubit Hamiltonian for Trotter evolution, equivalent to the intial Hamiltonian.
    Returns:
        h0t (Hamiltonian): 4 x 4 Hamiltonian used as a base for the initial Hamiltonian.
    """
    m0 = 0.5 * np.kron(matrices.I - matrices.X, matrices.I)
    return hamiltonians.Hamiltonian(2, m0, numpy=True)


def h_null():
    """Generate the 2 qubit Hamiltonian for Trotter evolution, equivalent to the intial Hamiltonian.
    Returns:
        h0t (Hamiltonian): 4 x 4 Hamiltonian used as a base for the initial Hamiltonian.
    """
    m0 = np.kron(matrices.I-matrices.I, matrices.I-matrices.I)
    return hamiltonians.Hamiltonian(2, m0, numpy=True)


def h0_parts(hp_parts):
    done = []
    h0_p = []
    for dic in hp_parts:
        dic2 = dic.copy()
        for key in set(dic2.keys()):
            if key[0] not in done:
                dic2[key] = h0t()
                done.append(key[0])
            else:
                dic2[key] = h_null()
        h0_p.append(dic2)
    return tuple(h0_p)
                
    
def ground_state(nqubits):
    """Returns |++...+> state to be used as the ground state of the easy Hamiltonian."""
    import tensorflow as tf
    from qibo.config import DTYPES
    s = np.ones(2 ** nqubits) / np.sqrt(2 ** nqubits)
    return tf.cast(s, dtype=DTYPES.get('DTYPECPX'))