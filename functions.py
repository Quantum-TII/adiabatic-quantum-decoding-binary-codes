import numpy as np
from sympy import symbols, expand, Matrix, Identity
from sympy.matrices.dense import matrix2numpy
#from qibo import matrices, hamiltonians, models, callbacks
import matplotlib.pyplot as plt
import collections
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from dwave.embedding.chain_strength import uniform_torque_compensation
from greedy import SteepestDescentSolver
import dimod
import dwave.inspector as insp
import matplotlib.pyplot as plt


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
    """Function for the 3-7 nonlinear code. (ANF)
    
    Args:
        x1, x2, x3 (symbol): original bits for the function.
        
    Returns:
        f (np.array): value for wach output bit of the nonlinear function.
        
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
    """Function for the 8-16 nonlinear code. (ANF)
    
    Args:
        x1, ..., x8 (symbol): original bits for the function.
        
    Returns:
        f (np.array): value for wach output bit of the nonlinear function.
        
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


def gadget(x, y, a):
    """Penalization function for the gadget substitution (x*y -> a).
    
    Args:
        x, y (symbol): qubits to be substituted.
        a (symbol): ancilla for the substitution.
        
    Returns:
        penalization function.
    
    """
    return 3*a + x*y - 2*x*a - 2*y*a


def substitutions(sym):
    """Substitutions needed as 0 or 1 squared are equal to themselves.
    
    Args:
        sym (list): list of qubits to reduce due to binary values.
        
    Returns:
        s (list): instructions for the substitutions.
    
    """
    s = []
    for x in sym:
        s.append((x**2, x))
    return s


def to_hamiltonian(f, r, sym):
    """Tranform the non-linear functions into a hamiltonian that minimizes the energy to find r and simplifies it.
    
    Args:
        f (np.array): functions defining the problem.
        r (np.array): target bitstring.
        sym (list): nstructions for the substitutions.
    
    Returns:
        h (symbol): hamiltonian encoding the solution of the problem.
    
    """
    h = np.sum((f-r)**2)
    h = expand(h)
    h = h.subs(substitutions(sym))
    return h


def separate_2_body(h, h_2):
    """Separate a hamiltonian between 1 and 2-order terms and the rest.
    
    Args:
        h (symbol): hamiltonian with more than 2-qubit terms.
        h_2 (symbol): hamiltonian with up to 2-qubit terms.
        
    Returns:
        h (symbol): updated hamiltonian with more than 2-qubit terms.
        h_2 (symbol): updated hamiltonian with up to 2-qubit terms.
        
    """
    terms_s, terms_c, overall_constant = symbolic_to_data(h)
    for i in range(len(terms_s)):
        if len(terms_s[i]) == 2:
            h -= terms_c[i]*terms_s[i][0]*terms_s[i][1]
            h_2 += terms_c[i]*terms_s[i][0]*terms_s[i][1]
        elif len(terms_s[i]) == 1:
            h -= terms_c[i]*terms_s[i][0]
            h_2 += terms_c[i]*terms_s[i][0]
    return h, h_2
    

def add_gadget(h, h_2, x1, x2, xa, control):
    """Add the gadget ancillas.
    
    Args:
        h (symbol): hamiltonian with more than 2-qubit terms.
        h_2 (symbol): hamiltonian with up to 2-qubit terms.
        x1, x2 (symbol): qubits to substitute.
        xa (symbol): ancilla qubit introduced.
        control (Bool): whether minimum control is desired.
        
    Returns:
        h (symbol): hamiltonian with more than 2-qubit terms.
        h_2 (symbol): hamiltonian with up to 2-qubit terms.
    
    """
    h, h_2 = separate_2_body(h, h_2)
    terms_s, terms_c, overall_constant = symbolic_to_data(h)
    h = expand(h.subs(x1*x2, xa))
    if control:
        d_m = 0
        d_p = 0
        for j in range(len(terms_s)):
            if x1 in terms_s[j] and x2 in terms_s[j]:
                if terms_c[j] < 0:
                    d_m += -terms_c[j]
                else:
                    d_p += terms_c[j]
        h += (1+max(d_m, d_p))*gadget(x1, x2, xa)
    else:
        for j in range(len(terms_s)):
            if x1 in terms_s[j] and x2 in terms_s[j]:
                h += (1+np.abs(terms_c[j]))*gadget(x1, x2, xa)
    return h, h_2


def to_gadget(h, x, ancillas, control):
    """Function that decomposes a multi-qubit interaction hamiltonian into 2-qubit using ancillas.
    
    Args:
        h (symbol): hamiltonian to decompose.
        x (list): original variables.
        ancillas (list): ancillas to use.
        control (Bool): whether minimum control is desired.
        
    Returns:
        h (symbol): Hamiltonian with up to 2-qubit interactions.
        anc (int): number of ancillas used.
        
    """
    h_2 = 0
    anc = 0
    h, h_2 = add_gadget(h, h_2, x[0], x[1], ancillas[anc], control)
    anc += 1
    if len(ancillas) == 1:
        return h+h_2, anc
    h, h_2 = add_gadget(h, h_2, x[2], x[3], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[0], x[2], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[1], x[3], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[0], x[3], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[1], x[2], ancillas[anc], control)
    anc += 1
    
    h, h_2 = add_gadget(h, h_2, x[0+4], x[1+4], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[2+4], x[3+4], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[0+4], x[2+4], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[1+4], x[3+4], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[0+4], x[3+4], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[1+4], x[2+4], ancillas[anc], control)
    anc += 1
    
    h, h_2 = add_gadget(h, h_2, x[0], x[9], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[1], x[9], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[2], x[8], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[3], x[8], ancillas[anc], control)
    anc += 1
    
    h, h_2 = add_gadget(h, h_2, x[0+4], x[9+6], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[1+4], x[9+6], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[2+4], x[8+6], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[3+4], x[8+6], ancillas[anc], control)
    anc += 1
    
    h, h_2 = add_gadget(h, h_2, x[8], x[9], ancillas[anc], control)
    anc += 1
    h, h_2 = add_gadget(h, h_2, x[8+6], x[9+6], ancillas[anc], control)
    anc += 1
    
    return h+h_2, anc


def check_interactions(h, high = False):
    """Check the number of terms for all different order of interactions in a hamiltonian.
    
    Args:
        h (symbol): hamiltonian to analyse the interactions.
        high (bool): flag to print the highest order interactions.
        
    Returns:
        terms (dict): dictionary with the different terms for each order of interactions.
    
    """
    terms_s, terms_c, overall_constant = symbolic_to_data(h)
    terms = {}
    for i in range(len(terms_s)):
        if len(terms_s[i]) not in terms.keys():
            terms[len(terms_s[i])] = 0
        terms[len(terms_s[i])] += 1
    if high:
        h = max(terms.keys())
        for i in range(len(terms_s)):
            if len(terms_s[i]) == h:
                print(terms_s[i])
        print()
    return terms


def check_two_body(h):
    """Check if given Hamiltonian only contains up to 2-body interactions.
    
    Args:
        h (symbol): Hamiltonian to be analysed.
        
    Returns:
        ValueError: if there are terms with more than two body interactions.
    
    """
    terms_s, terms_c, overall_constant = symbolic_to_data(h)
    length = []
    for i in range(len(terms_s)):
        if len(terms_s[i]) > 2:
            raise ValueError(f'Error! Hamiltonian with more than two-body interactions. Found at least a {len(terms_s[i])}-body interaction.')

    
def symbolic_to_data(symbolic_hamiltonian):
    """Transforms a symbolic Hamiltonian to lists of every term.
        
    Args:
        symbolic_hamiltonian: The full Hamiltonian written with symbols.
    
    Returns:
        matching list of the symbols in each term of the hamiltonian and the corresponding constant.
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
    
    
def symbolic_to_dwave(symbolic_hamiltonian, symbol_num):
    """Transforms a symbolic Hamiltonian to a dictionary of targets and matrices.
    
    Works for Hamiltonians with one and two qubit terms only.
    
    Args:
        symbolic_hamiltonian: The full Hamiltonian written with symbols.
        symbol_num: Dictionary that maps each symbol that appears in the 
            Hamiltonian to its target.
    
    Returns:
       Q (dict): Dictionary with the interactions to send to the DWAVE machine.
       overall_constant (int): Constant that cannot be given to DWAVE machine.
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
            #print(symbols[0], target)
            Q[(target, target)] = constant

        elif len(symbols) == 2:
            target1 = symbol_num[symbols[0]]
            target2 = symbol_num[symbols[1]]
            Q[(target1, target2)] = constant

        else:
            raise ValueError("Only one and two qubit terms are allowed.")
    
    return Q, overall_constant


def dwave(h, h2, sym, symbol_num, bits, T, chainstrength, numruns, inspect, repeat, iterations):
    """Function to run the problem with the dwave hardware.
    
    Args:
        h (symbol): hamiltonian where the solution is encoded.
        h2 (symbol): hamiltonian to recover final energy.
        sym (list): symbols to perfom the substitution.
        symbol_num (dict): dictionary with the pairings between symbols and qubits.
        bits (int): number of bits of the initial bitstring.
        T (float): Total annealing time.
        chainstrength (float): Chainstrangth for device. Leave None if calculated automatically.
        numruns (int): Number of samples to take from quantum computer.
        inspect (Bool): Open the dwave inspector.
        repeat (Bool): Fix ancillas using repetition.
        iterations (int): number of repetitions for the ancilla fixing.
        
    Returns:
        best_sample (list): bit positions of the best found solution.
        best_energy (float): energy value associated to the found minimum.
        
    """
    check_two_body(h)
    Q, constant = symbolic_to_dwave(h, symbol_num)
    model = dimod.BinaryQuadraticModel.from_qubo(Q, offset = 0.0)
    if not chainstrength:
        chainstrength = 0
        for i in Q.values():
            chainstrength += abs(i)
        chainstrength *= 3/len(Q)
        print(f'Automatic chain strength: {chainstrength}\n')
    else:
        print(f'Chosen chain strength: {chainstrength}\n')
    sampler = EmbeddingComposite(DWaveSampler())
    response = sampler.sample(model, chain_strength=chainstrength, num_reads=numruns, annealing_time=T, answer_mode='histogram')
    record = response.record
    order = np.argsort(record['energy'])
    best_sample = record.sample[order[0]]
    best_energy = record.energy[order[0]]+constant
    print(f'Best result found:\n')
    print(f'Relevant bits: {best_sample[:bits]}\n')
    print(f'Ancillas: {best_sample[bits:]}\n')
    print(f'With energy: {best_energy}\n')
    print(f'Occurences: {record.num_occurrences[order[0]]}\n')
    print(f'The best {min(len(record.sample), bits)} samples found in the evolution are:\n')
    for i in range(min(len(record.sample), bits)):
        print(f'Bits: {record.sample[order[i]][:bits]}    Ancillas: {record.sample[order[i]][bits:]}    with energy: {record.energy[order[i]]+constant}    num. occurences: {record.num_occurrences[order[i]]}\n')
    if inspect:
        insp.show(response)
    energy = []
    w = 0
    for i in range(len(record.energy)):
        for j in range(record.num_occurrences[order[i]]):
            energy.append(record.energy[order[i]]+constant)
    if repeat:
        best_sample, best_energy, w = dwave_iterative(h, h2, sym, record, order, symbol_num, bits, T, chainstrength, numruns, iterations)
        return best_sample, best_energy, w
    else:
        return best_sample, best_energy, w


def dwave_iterative(h, h2, sym, record, order, symbol_num, bits, T, chainstrength, numruns, iterations):
    """Iterative method that fixes qubits that are thought to be found in the best position.
    
    Args:
        h (symbol): hamiltonian where the solution is encoded.
        h2 (symbol): hamiltonian to recover final energy.
        sym (list): symbols to perfom the substitution.
        record: last result from dwave.
        order (list): order from small to large energy.
        symbol_num (dict): dictionary with the pairings between symbols and qubits.
        bits (int): number of bits of the initial bitstring.
        T (float): Total annealing time.
        chainstrength (float): Chainstrangth for device. Leave None if calculated automatically.
        numruns (int): Number of samples to take from quantum computer.
        iterations (int): number of repetitions for the ancilla fixing.
        
    Returns:
        result (list): reconstructed best solution found after iterating.
        h (float): energy of theh system using the result output.
        w (int): iteration when the first solution is found.
        
    """
    fix = []
    out = []
    for w in range(iterations):
        c = bits
        for j in range(bits, len(sym)):
            if j in fix:
                continue
            else:
                a = True
                b = record.sample[order[0]][c]
                for k in range(min(len(record.sample), bits)):#, int(np.ceil(np.log2(len(sym)))))):
                    if b != record.sample[order[k]][c]:
                        a = False
                if a:
                    fix.append(j)
                    out.append(b)
                c += 1
        print(f'The same value was found in positions {fix} \n')
        print(f'with values {out}.\n')
        for j in range(len(fix)):
            h = h.subs(sym[fix[j]], out[j])
        terms = check_interactions(h, high=False)
        print(f'Total number of qubits needed for the next step: {len(sym)-len(fix)}.\n')
        print('Number of terms for each k-body interactions after gadget aplication.\n')
        print(terms, '\n')
        Q, constant = symbolic_to_dwave(h, symbol_num)
        model = dimod.BinaryQuadraticModel.from_qubo(Q, offset = 0.0)
        if not chainstrength:
            chainstrength = 0
            for i in Q.values():
                chainstrength += abs(i)
            chainstrength *= 3/len(Q)
            print(f'Automatic chain strength: {chainstrength}\n')
        else:
            print(f'Chosen chain strength: {chainstrength}\n')
        sampler = EmbeddingComposite(DWaveSampler())
        response = sampler.sample(model, chain_strength=chainstrength, num_reads=numruns, annealing_time=T, answer_mode='histogram')
        record = response.record
        order = np.argsort(record['energy'])
        best_sample = record.sample[order[0]]
        best_energy = record.energy[order[0]]
        print(f'Best result found: {best_sample}\n')
        print(f'With energy: {best_energy+constant}\n')
        print(f'Occurences: {record.num_occurrences[order[0]]}\n')
        print(f'The best {min(len(record.sample), int(np.ceil(np.log2(len(sym)))))} samples found in the evolution are:\n')
        for i in range(min(len(record.sample), bits)):#, int(np.ceil(np.log2(len(sym)))))):
            print(f'Result: {record.sample[order[i]]}    with energy: {record.energy[order[i]]+constant}    num. occurences: {record.num_occurrences[order[i]]}\n')
        energy = []
        for i in range(len(record.energy)):
            for j in range(record.num_occurrences[order[i]]):
                energy.append(record.energy[order[i]]+constant)
        if best_energy+constant == 2:
            print('Solution found!\n')
            break
    print('Reconstructing state...\n')
    c = 0
    result = []
    for i in range(len(sym)):
        if i in fix:
            result.append(out[fix.index(i)])
        else:
            result.append(best_sample[c])
            c += 1
    print(f'Reconstructed result:\n')
    print(f'Relevant bits: {result[:bits]}\n')
    print(f'Ancillas: {result[bits:]}\n')
    for i in range(bits):
        h2 = h2.subs(sym[i], result[i])
    print(f'With total energy: {h2}\n')
    return result, h2, w
