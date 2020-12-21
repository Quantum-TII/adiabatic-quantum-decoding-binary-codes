#!/usr/bin/env python
from sympy import symbols
import numpy as np
import functions
import argparse
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from dwave.embedding.chain_strength import uniform_torque_compensation
from greedy import SteepestDescentSolver
import dimod
import dwave.inspector


def main(bits, T, chainstrength, numruns, greedy, inspect):
    """

    Args:
        nqubits (int): number of qubits for the file that contains the
            information of an Exact Cover instance.
        instance (int): intance used for the desired number of qubits.
        T (float): 
        

    Returns:
        
    """
    x = symbols(' '.join((f'x{i}' for i in range(1, bits+1))))
    if bits == 3:
        ancillas = (symbols(' '.join((f'xa{i}' for i in range(1, 2)))), )
        f = functions.f_3(*x)
        r = functions.r_3()
    elif bits == 8:
        ancillas = symbols(' '.join((f'xa{i}' for i in range(1, 23))))
        f = functions.f_8(*x)
        r = functions.r_8()
    else:
        raise ValueError('Only instances for 3 and 8 bits are supported.\n')
        
    sym = x+ancillas
    symbol_num = {}
    for i in range(len(sym)):
        symbol_num[sym[i]] = i
    print(f'Solving a {bits} bit instance of a nonlinear function using Adiabatic Evolution.\n')
    print('Creating problem Hamiltonian.\n')
    h = functions.to_hamiltonian(f, r, sym)
    terms = functions.check_interactions(h)
    print('Number of terms for each k-body interactions.\n')
    print(terms, '\n')
    print('Reducing the k-body interactions using ancillas.\n')
    if control:
        h2, anc = functions.to_gadget_ruge2(h, sym, ancillas, bits)
    else:
        h2, anc = functions.to_gadget_ruge1(h, sym, ancillas, bits)
        
        
    terms = functions.check_interactions(h2, high=False)
    print(f'Total number of qubits needed for the 2-local hamiltonian encoding {bits} bits: {anc+bits}.\n')
    print('Number of terms for each k-body interactions after gadget aplication.\n')
    print(terms, '\n')
    
    
    functions.check_two_body(h2)
    Q, constant = functions.symbolic_to_dwave(h2, symbol_num)

    model = dimod.BinaryQuadraticModel.from_qubo(Q, offset = 0.0)
    if not chainstrength:
        chainstrength = dwave.embedding.chain_strength.uniform_torque_compensation(model)
        print(f'Automatic chain strength: {chainstrength}\n')
    else:
        print(f'Chosen chain strength: {chainstrength}\n')
        
    sampler = EmbeddingComposite(DWaveSampler())
    if greedy:
        solver_greedy = SteepestDescentSolver()
        sampleset = sampler.sample(model, chain_strength=chainstrength, num_reads=numruns, annealing_time=T, answer_mode='raw')
        response = solver_greedy.sample(model, initial_states=sampleset)
    else:
        response = sampler.sample(model, chain_strength=chainstrength, num_reads=numruns, annealing_time=T, answer_mode='histogram')
    
    
    best_sample = response.record.sample[0]
    best_energy = response.record.energy[0]
    print(f'Best result found: {best_sample}\n')
    print(f'With energy: {best_energy+constant}\n')
    good_samples = response.record.sample[:min(len(response.record.sample), nqubits)]
    good_energies = response.record.energy[:min(len(response.record.energy), nqubits)]
    print(f'The best {len(good_samples)} samples found in the evolution are:\n')
    for i in range(len(good_samples)):
        print(f'Sample: {good_samples[i]}    with energy: {good_energies[i]+constant}\n')

    if inspect:
        dwave.inspector.show(response)

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", default=8, type=int)
    parser.add_argument("--T", default=20, type=float)
    parser.add_argument("--chainstrength", default=None, type=float)
    parser.add_argument("--numruns", default=100, type=int)
    parser.add_argument("--greedy", action="store_true")
    parser.add_argument("--inspect", action="store_true")
    args = vars(parser.parse_args())
    main(**args)