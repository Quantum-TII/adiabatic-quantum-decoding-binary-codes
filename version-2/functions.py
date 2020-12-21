import numpy as np
from sympy import symbols, expand, Matrix, Identity
from sympy.matrices.dense import matrix2numpy
from qibo import matrices, hamiltonians, models, callbacks
import matplotlib.pyplot as plt
import collections
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from dwave.embedding.chain_strength import uniform_torque_compensation
from greedy import SteepestDescentSolver
import dimod
import dwave.inspector


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
    f = np.array([x1*x2*x3 + x1*x2 + x1*x3 + x1 + x2*x3 + x2, 
                  x1*x2 + x1*x3 + x3,
                  x1*x2*x3 + x1*x2 + x2 + x3,  
                  x1*x2*x3 + x1*x3 + x2*x3 + 1, 
                  x1*x2*x3 + x1*x3 + x1 + x2*x3 + 1, 
                  x1 + x2*x3 + x2, 
                  x1*x3 + x2 + x3 + 1])
    return f


def f_8(x1, x2, x3, x4, x5, x6, x7, x8):
    """Function for the 8-16 nonlinear code. (ANF)
    
    Args:
        x1, ..., x8 (symbol): original bits for the function.
        
    Returns:
        f (np.array): value for wach output bit of the nonlinear function.
        
    """
    f = np.array([x1 + x2,
                  x2 + x3,
                  x3 + x4,
                  x4 + x5,
                  x5 + x6*x7*x8 + 1,
                  x2 + x4 + x5 + x6*x7*x8 + x6*x8 + x7 + x8,
                  x1 + x2 + x4 + x6*x7*x8 + x6*x8 + x7*x8 + x7 + x8 + 1,
                  x1 + x2 + x3 + x4 + x5 + x6*x7*x8 + x6*x7 + x6*x8 + x7*x8,
                  x2 + x5 + x6*x7 + x8,
                  x1 + x6*x7*x8 + x6*x8 + x6 + x7*x8 + x7 + x8,
                  x1 + x4 + x6*x7 + x6*x8 + x6 + x7,
                  x1 + x3 + x5 + x6*x7*x8 + x6*x8 + x6 + x7*x8 + x8 + 1,
                  x2 + x3 + x4 + x6*x7*x8 + x6*x8 + x6 + x7*x8,
                  x1 + x3 + x4 + x5 + x6*x8 + x7*x8 + x7 + 1,
                  x1 + x2 + x3 + x5 + x6 + x7*x8 + x8 + 1,
                  x3 + x6*x7*x8 + x6*x7 + x6 + x7*x8 + x7 + 1])
    return f


def r_3():
    """Target state for the 3-7 nonlinear code.
    
    """
    return np.array([0, 1, 0, 1, 0, 0, 0])


def r_8():
    """Target state for the 8-16 nonlinear code.
    
    """
    return np.array([1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0])


def cnot_penalization(xt, xc, xr, xa, p_cnot):
    """Penalization function for a CNOT gate. 
    
    Args:
        xt (symbol): target qubit.
        xc (symbol): control qubit.
        xr (symbol): qubit where the result of the CNOT is stored.
        xa (symbol): ancilla qubit to break 3-body interactions.
        p_cnot (float): constant multiplying the penalization term.
        
    Returns:
        Term penalizing combinations that do not represent a CNOT gate.
        
    """
    return p_cnot*(2*xc*xt - 2*xc*xr - 2*xt*xr - 4*xc*xa - 4*xt*xa + 4*xr*xa + xc + xt + xr + 4*xa)


def and_penalization(x1, x2, x12, p_cnot):
    """Penalization function for an AND gate.
    
    Args:
        x1 (symbol): first qubit.
        x2 (symbol): second qubit.
        x12 (symbol): multiplication of the first and second qubits.
        p_cnot (float): constant multiplying the penalization term.
        
    Returns:
        Term penalizing combinations that do not represent an AND gate.
        
    """
    return p_cnot*(x1*x2 - 2*x1*x12 - 2*x2*x12 + 3*x12)


def twoffoli_penalization(xt, xc1, xc2, xr, xb, xa, p_cnot):
    """Penalization function for a Toffoli gate.
    
    Args:
        xt (symbol): target qubit.
        xc1 (symbol): control qubit.
        xc2 (symbol): control qubit.
        xr (symbol): qubit where the result of the Toffoli is stored.
        xb (symbol): multiplication of the two control qubits.
        xa (symbol): ancilla qubit to break 3-body interactions.
        p_cnot (float): constant multiplying the penalization term.
        
    Returns:
        Term penalizing combinations that do not represent a Toffoli gate.
        
    """
    pen = cnot_penalization(xt, xb, xr, xa, p_cnot)
    pen += and_penalization(xc1, xc2, xb, p_cnot)
    return pen


def const_penalization(xt, xr, p_cnot):
    """Penalization function for an X gate.
    
    Args:
        xt (symbol): target qubit.
        xr (symbol): qubit after the X gate is performed.
        xa (symbol): ancilla qubit to break 3-body interactions.
        p_cnot (float): constant multiplying the penalization term.
        
    Returns:
        Term penalizing combinations that do not represent an X gate.
        
    """
    return p_cnot*(2*xt*xr - xt - xr + 1)


def threeffoli_penalization(xt, xc1, xc2, xc3, xc, xr, xb, xa, p_cnot):
    """Penalization function for a Toffoli gate with 3 controls.
    
    Args:
        xt (symbol): target qubit.
        xc1 (symbol): control qubit.
        xc2 (symbol): control qubit.
        xc2 (symbol): control qubit.
        xc (symbol): multiplication of the first two controls.
        xr (symbol): qubit where the result of the Toffoli with 3 controls is stored.
        xb (symbol): accounts for the product of all controls.
        xa (symbol): ancilla qubit to break 3-body interactions.
        p_cnot (float): constant multiplying the penalization term.
        
    Returns:
        Term penalizing combinations that do not represent a Toffoli gate with 3 controls.
        
    """
    pen = and_penalization(xc1, xc2, xc, p_cnot)
    pen += and_penalization(xc, xc3, xb, p_cnot)
    pen += cnot_penalization(xt, xb, xr, xa, p_cnot)
    return pen


def create_outputs(f, r):
    """Assign the variables that are going to be used as outputs throughout the computation. 
    
    Args:
        f (np.array): functions that map the original to the target bitstring space.
        r (np.array): target bitstring.
        
    Returns:
        outputs (list): extra qubits needed to store the results of applying gates throughout the computation.
        
    """
    outputs = []
    o = 1
    for i in range(len(f)):
        out = []
        xo = symbols(f'xo{o}')
        out.append(xo)
        o += 1
        for j in range(len(f[i].args)):
            xo = symbols(f'xo{o}')
            out.append(xo)
            o += 1
        if r[i] == 1:
            xo = symbols(f'xo{o}')
            out.append(xo)
            o += 1
        outputs.append(out)
    print(f'Total output qubits used: {o}\n')
    return outputs


def create_hamiltonian(f, r, outputs, x, p_sol, p_cnot):
    """Creation of the Hamiltonian using a polynomal number of ancilla qubits with the number of terms in the function.
    
    Args:
        f (np.array): functions that map the original to the target bitstring space.
        r (np.array): target bitstring.
        outputs (list): qubits used for outputs when applying the gate penalization functions.
        x (tuple): qubits that encode the solution of the problem.
        p_sol (float): penalization applied to the checking of the solution.
        p_cnot (float): penalization to the gates penalization functions. 
    
    Returns:
        h (symbol): hamiltonian that encodes in its ground state the solution of the problem.
        anc (list): ancillas used to create the hamiltonian.
        results_o (list): expected result for the result qubits given initial state.
        results_a (list): expected result for the ancillas given initial state.
        
    """
    h = 0
    anc = []
    ancillas = {}
    c = 1
    if len(x) == 3:
        results = {x[0]: 0, x[1]: 0, x[2]: 0}
    elif len(x) == 8:
        results = {x[0]: 0, x[1]: 1, x[2]: 0, x[3]: 0, x[4]: 0, x[5]: 1, x[6]: 0, x[7]: 1}
    results_o = create_outputs(f, r)
    results_a = []
    for i, row in enumerate(f):
        results_o[i][0] = 0
        h += p_cnot*outputs[i][0]
        for j, term in enumerate(row.args):
            if not term.args:
                expression = (term,)
            else:
                expression = term.args
            symbol = [x for x in expression if x.is_symbol]
            numbers = [x for x in expression if not x.is_symbol]

            if len(numbers) > 1:
                raise ValueError("Hamiltonian must be expanded before using this method.")
            elif numbers:
                constant = float(numbers[0])
            else:
                constant = 1
            
            if len(symbol) == 0:
                h += const_penalization(outputs[i][j], outputs[i][j+1], p_cnot)
                results_o[i][j+1] = (results_o[i][j]+1)%2
            elif len(symbol) == 1:
                if ancillas.get(symbol[0]*outputs[i][j]) is None:
                    a = symbols(f'xa{c}')
                    anc.append(a)
                    ancillas[symbol[0]*outputs[i][j]] = a
                    results_a.append(results[symbol[0]]*results_o[i][j])
                    c += 1
                h += cnot_penalization(outputs[i][j], symbol[0], outputs[i][j+1], ancillas[symbol[0]*outputs[i][j]], p_cnot)
                results_o[i][j+1] = (results_o[i][j]+results[symbol[0]])%2
            elif len(symbol) == 2:
                if ancillas.get(symbol[0]*symbol[1]) is None:
                    a = symbols(f'xa{c}')
                    anc.append(a)
                    ancillas[symbol[0]*symbol[1]] = a
                    results_a.append(results[symbol[0]]*results[symbol[1]])
                    c += 1
                if ancillas.get(ancillas[symbol[0]*symbol[1]]*outputs[i][j]) is None:
                    a = symbols(f'xa{c}')
                    anc.append(a)
                    ancillas[ancillas[symbol[0]*symbol[1]]*outputs[i][j]] = a
                    results_a.append(results[symbol[0]]*results[symbol[1]]*results_o[i][j])
                    c += 1
                h += twoffoli_penalization(outputs[i][j], symbol[0], symbol[1], outputs[i][j+1], ancillas[symbol[0]*symbol[1]], ancillas[ancillas[symbol[0]*symbol[1]]*outputs[i][j]], p_cnot)
                results_o[i][j+1] = (results_o[i][j]+(results[symbol[0]]*results[symbol[1]]))%2
            elif len(symbol) == 3:
                if ancillas.get(symbol[0]*symbol[1]) is None:
                    a = symbols(f'xa{c}')
                    anc.append(a)
                    ancillas[symbol[0]*symbol[1]] = a
                    results_a.append(results[symbol[0]]*results[symbol[1]])
                    c += 1
                if ancillas.get(ancillas[symbol[0]*symbol[1]]*symbol[2]) is None:
                    a = symbols(f'xa{c}')
                    anc.append(a)
                    ancillas[ancillas[symbol[0]*symbol[1]]*symbol[2]] = a
                    results_a.append(results[symbol[0]]*results[symbol[1]]*results[symbol[2]])
                    c += 1
                if ancillas.get(ancillas[ancillas[symbol[0]*symbol[1]]*symbol[2]]*outputs[i][j]) is None:
                    a = symbols(f'xa{c}')
                    anc.append(a)
                    ancillas[ancillas[ancillas[symbol[0]*symbol[1]]*symbol[2]]*outputs[i][j]] = a
                    results_a.append(results[symbol[0]]*results[symbol[1]]*results[symbol[2]]*results_o[i][j])
                    c += 1
                h += threeffoli_penalization(outputs[i][j], symbol[0], symbol[1], symbol[2], ancillas[symbol[0]*symbol[1]], outputs[i][j+1], ancillas[ancillas[symbol[0]*symbol[1]]*symbol[2]], ancillas[ancillas[ancillas[symbol[0]*symbol[1]]*symbol[2]]*outputs[i][j]], p_cnot)
                results_o[i][j+1] = (results_o[i][j]+(results[symbol[0]]*results[symbol[1]]*results[symbol[2]]))%2
            else:
                print('ouch')
        if r[i] == 1:
            h += const_penalization(outputs[i][j+1], outputs[i][j+2], p_cnot)
            results_o[i][j+2] = (results_o[i][j+1]+1)%2
            print(outputs[i][-1])
            h += p_sol*outputs[i][j+1]
        else:
            print(outputs[i][-1])
            h += p_sol*outputs[i][j+1]
    print(f'Total number of ancillas added: {c}\n')
    return expand(h), anc, results_o, results_a


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
            #print(symbols[0], target1)
            #print(symbols[1], target2)
            Q[(target1, target2)] = constant

        else:
            raise ValueError("Only one and two qubit terms are allowed.")
    
    return Q, overall_constant


def dwave(h, symbol_num, bits, T, chainstrength, numruns, greedy, inspect, repeat, iterations):
    """Function to run the problem with the dwave hardware.
    
    Args:
        h (symbol): hamiltonian where the solution is encoded.
        symbol_num (dict): dictionary with the pairings between symbols and qubits.
        bits (int): number of bits of the initial bitstring.
        T (float): Total annealing time.
        chainstrength (float): Chainstrangth for device. Leave None if calculated automatically.
        numruns (int): Number of samples to take from quantum computer.
        greedy (Bool): process results using the greedy algorithm.
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
    if greedy:
        solver_greedy = SteepestDescentSolver()
        sampleset = sampler.sample(model, chain_strength=chainstrength, num_reads=numruns, annealing_time=T, answer_mode='raw')
        response = solver_greedy.sample(model, initial_states=sampleset)
    else:
        response = sampler.sample(model, chain_strength=chainstrength, num_reads=numruns, annealing_time=T, answer_mode='histogram')
    record = response.record
    order = np.argsort(record['energy'])
    best_sample = record.sample[order[0]]
    best_energy = record.energy[order[0]]+constant
    print(f'Best result found:\n')
    print(f'Relevant bits: {best_sample[:bits]}\n')
    print(f'Ancillas: {best_sample[bits:]}\n')
    print(f'With energy: {best_energy}\n')
    print(f'The best {min(len(record.sample), bits)} samples found in the evolution are:\n')
    for i in range(min(len(record.sample), bits)):
        print(f'Bits: {record.sample[order[i]][:bits]}    Ancillas: {record.sample[order[i]][bits:]}    with energy: {record.energy[order[i]]+constant}\n')
    if inspect:
        dwave.inspector.show(response)
    if repeat:
        best_sample, best_energy = dwave_iterative(h, symbol_num, bits, T, chainstrength, numruns, iterations)
    return best_sample, best_energy


def dwave_iterative(h, symbol_num, bits, T, chainstrength, numruns, iterations):
    """Iterative method that fixes qubits that are thought to be found in the best position.
    
    Args:
        h (symbol): hamiltonian where the solution is encoded.
        symbol_num (dict): dictionary with the pairings between symbols and qubits.
        bits (int): number of bits of the initial bitstring.
        T (float): Total annealing time.
        chainstrength (float): Chainstrangth for device. Leave None if calculated automatically.
        numruns (int): Number of samples to take from quantum computer.
        iterations (int): number of repetitions for the ancilla fixing.
        
    Returns:
        result (list): reconstructed best solution found after iterating.
        h (float): energy of theh system using the result output.
        
    """
    fix = []
    out = []
    for i in range(iterations):
        c = bits
        for j in range(bits, len(sym)):
            if j in fix:
                continue
            else:
                a = True
                b = record.sample[order[0]][c]
                for k in range(min(len(record.sample), np.ceil(np.log2(len(sym))))):
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
        terms = functions.check_interactions(h, high=False)
        print(f'Total number of qubits needed for the next step: {len(sym)-len(fix)}.\n')
        print('Number of terms for each k-body interactions after gadget aplication.\n')
        print(terms, '\n')
        Q, constant = functions.symbolic_to_dwave(h, symbol_num)
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
        print(f'The best {min(len(record.sample), bits)} samples found in the evolution are:\n')
        for i in range(min(len(record.sample), bits)):
            print(f'Result: {record.sample[order[i]]}    with energy: {record.energy[order[i]]+constant}\n')
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
        h = h.subs(sym[i], result[i])
    print(f'With total energy: {h}\n')
    return result, h
