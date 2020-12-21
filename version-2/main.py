#!/usr/bin/env python
from sympy import symbols
import numpy as np
import functions
import argparse


def main(bits, p_sol, p_cnot, T, chainstrength, numruns, solution, dwave, greedy, inspect, repeat, iterations):
    """Adiabatic evolution for nonlinear functions.

    Args:
        bits (int): Number of bits of the example.
        p_sol (float): Penalization of wrong solution.
        p_cnot (float): Penalization of wrong gate.
        T (float): Total annealing time.
        chainstrength (float): Chainstrangth for device. Leave None if calculated automatically.
        numruns (int): Number of samples to take from quantum computer.
        solution (Bool): Keep track and print the solution and its energy.
        dwave (Bool): run everything on dwave.
        greedy (Bool): process results using the greedy algorithm.
        inspect (Bool): Open the dwave inspector.
        repeat (Bool): Fix ancillas using repetition.
        iterations (int): number of repetitions for the ancilla fixing.

    -------------- if minimization is run via dwave -------------------
    Returns:
        best_sample (list): bit positions of the best found solution.
        best_energy (float): energy value associated to the found minimum.
        
    """
    x = symbols(' '.join((f'x{i}' for i in range(1, bits+1))))
    if bits == 3:
        f = functions.f_3(*x)
        r = functions.r_3()
        outputs = functions.create_outputs(f, r)
        results = (1, 1, 0)
    elif bits == 8:
        f = functions.f_8(*x)
        r = functions.r_8()
        outputs = functions.create_outputs(f, r)
        results = (0, 1, 0, 0, 0, 1, 0, 1)
    else:
        raise ValueError('Only instances for 3 and 8 bits are supported.\n')
    print(f'Solving a {bits} bit instance of a nonlinear function using Adiabatic Evolution.\n')
    print('Creating problem Hamiltonian.\n')
    h, anc, results_o, results_a = functions.create_hamiltonian(f, r, outputs, x, p_sol, p_cnot)
    terms = functions.check_interactions(h)
    print('Number of terms for each k-body interactions.\n')
    print(terms, '\n')
    sym = x
    for i in range(len(outputs)):
        sym += tuple(outputs[i])
        results += tuple(results_o[i])
    sym += tuple(anc)
    results += tuple(results_a)
    symbol_num = {}
    for i in range(len(sym)):
        symbol_num[sym[i]] = i
    print(f'Total number of qubits needed for the 2-local hamiltonian encoding {bits} bits: {len(sym)}.\n')
    if solution:
        s = []
        for j in range(len(sym)):
            s.append((sym[j], results[j]))
        groundstate = h.subs(s)
        print('Checking the known groundstate with the ancilla qubits:\n')
        print(f'    State: {val}')
        print(f'    Energy: {groundstate}\n')
    if dwave:
        best_sample, best_energy = dwave(h, symbol_num, bits, T, chainstrength, numruns, greedy, inspect, repeat, iterations)
        return best_sample, best_energy


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", default=3, type=int)
    parser.add_argument("--p_sol", default=1, type=float)
    parser.add_argument("--p_cnot", default=3, type=float)
    parser.add_argument("--T", default=20, type=float)
    parser.add_argument("--chainstrength", default=None, type=float)
    parser.add_argument("--numruns", default=500, type=int)
    parser.add_argument("--solution", action="store_true")
    parser.add_argument("--dwave", action="store_true")
    parser.add_argument("--greedy", action="store_true")
    parser.add_argument("--inspect", action="store_true")
    parser.add_argument("--repeat", action="store_true")
    parser.add_argument("--iterations", default=5, type=int)
    args = vars(parser.parse_args())
    main(**args)
