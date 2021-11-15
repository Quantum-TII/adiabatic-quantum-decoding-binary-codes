import numpy as np
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import dimod
import dwave.inspector as insp
import itertools


class dwave_annealing:
    def __init__(self, H, bits, sym, numruns=1000, T=20, chainstrength=None):
        """Set up the variables needed for the class.

        """
        self.H = H
        self.bits = bits
        self.sym = sym
        self.numruns = numruns
        self.T = T
        self.chainstrength = chainstrength
        self.energies_histogram = []

    def set_numruns(self, numruns):
        """Allows for change of the number of runs of the class.

        Args: 
            numruns (int): number of measurements that the device will take.

        """
        self.numruns = numruns

    def get_symbol_num(self, sym):
        """Get the dictionary that maps each symbol to the corresponding logical qubit.

        Args:
            sym (tuple): symbols that make up the Hamiltonian.

        Returns:
            symbol_num (dict): dictionary that assignes a logical qubit to each symbol variable.

        """
        symbol_num = {}
        for i, s in enumerate(sym):
            symbol_num[s] = i
        return symbol_num

    def symbolic_to_dwave(self, symbolic_hamiltonian, symbol_num):
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
                Q[(target, target)] = constant
            elif len(symbols) == 2:
                target1 = symbol_num[symbols[0]]
                target2 = symbol_num[symbols[1]]
                Q[(target1, target2)] = constant
            else:
                raise ValueError("Only one and two qubit terms are allowed.")
        
        return Q, overall_constant

    def run_dwave(self, H, sym):
        """Run a quantum annealing protocol using D-Wave's devices.

        Args:
            H (symbol): Hamiltonian that encodes the solution of a problem.
            sys (tuple): symbols that make up the Hamiltonian.

        Returns:
            response: object that contains the details of the experiment.

        """
        chainstrength = self.chainstrength
        numruns, T = self.numruns, self.T
        symbol_num = self.get_symbol_num(sym)
        Q, constant = self.symbolic_to_dwave(H, symbol_num)
        model = dimod.BinaryQuadraticModel.from_qubo(Q, offset = constant)

        if not chainstrength:
            chainstrength = 0
            for i in Q.values():
                chainstrength += abs(i)
            chainstrength *= 3/len(Q)
            #print(f'Automatic chain strength: {chainstrength}\n')
        #else:
            #print(f'Chosen chain strength: {chainstrength}\n')

        sampler = EmbeddingComposite(DWaveSampler())
        response = sampler.sample(model, chain_strength=chainstrength, num_reads=numruns, annealing_time=T, answer_mode='histogram')
        return response

    def get_energies(self, response):
        """Get data for the energies mesured after a dwave experiment.

        Args:
            response: result of a dwave run.

        Returns:
            best_sample (list): sample that yielded the least energy.
            best_energy (float): lowest energy measured in the experiment.

        """
        record = response.record
        order = np.argsort(record['energy'])

        best_sample = record.sample[order[0]]
        best_energy = record.energy[order[0]]

        energy = []
        for i in range(len(record.energy)):
            for j in range(record.num_occurrences[order[i]]):
                energy.append(record.energy[order[i]])

        self.energies_histogram.append(energy)

        #Printing results
        print(f'Best result found:\n')
        print(f'Relevant bits: {best_sample[:self.bits]}\n')
        print(f'Ancillas: {best_sample[self.bits:]}\n')
        print(f'With energy: {best_energy}\n')
        print(f'Occurences: {record.num_occurrences[order[0]]}\n')

        print(f'The best {min(len(record.sample), self.bits)} samples found in the evolution are:\n')
        for i in range(min(len(record.sample), self.bits)):
            print(f'Bits: {record.sample[order[i]][:self.bits]}    Ancillas: {record.sample[order[i]][self.bits:]}    with energy: {record.energy[order[i]]}    num. occurences: {record.num_occurrences[order[i]]}\n')

        return best_sample, best_energy

    def clear_energies_histogram(self):
        """Clear the saved energies to start a new collection.

        """
        self.energies_histogram = []

    def get_embedding(self, response):
        """Get data for the logical and physical qubits needed for the embeding of a single response.

        Args:
            response: result of a dwave run.

        Returns:
            log_q (float): logical qubits needed for the embedding.
            phys_q (float): physical qubits needed for the embedding.

        """
        embedding = response.info['embedding_context']['embedding']

        log_q = len(embedding.keys())
        phys_q = sum(len(chain) for chain in embedding.values())

        print(f"Number of logical variables: {log_q}")
        print(f"Number of physical qubits used in embedding: {phys_q}\n")

        return log_q, phys_q

    def get_phys_qubits(self, H, sym, average=1):
        """Get an average of the logical and physical qubits needed for the
        embeding of the problem.

        Args:
            H (symbol): symbolic Hamiltonian that wants to be solved.
            sym (tuple): symbols that make up the Hamiltonian.
            average (int): number of instances to average the qubit numbers.

        Returns:
            log_q (float): average of the logical qubits needed for the embedding.
            phys_q (float): average of the physical qubits needed for the embedding.

        """
        log_q = 0
        phys_q = 0
        for _ in range(average):
            response = self.run_dwave(H, sym)
            lq, pq = self.get_embedding(response)
            log_q += lq
            phys_q += pq

        return log_q/average, phys_q/average

    def dwave_single_run(self):
        """Function to run the problem with the dwave hardware.
        
        Args:
            H (symbol): hamiltonian where the solution is encoded.
            sym (list): symbols to perfom the substitution.
            
        Returns:
            best_sample (list): bit positions of the best found solution.
            best_energy (float): energy value associated to the found minimum.
            
        """
        response = self.run_dwave(self.H, self.sym)

        best_sample, best_energy = self.get_energies()

        if self.inspect:
            insp.show(response)

        if best_energy == 0:
            print("Solution found!\n")
        
        return best_sample, best_energy

    def dwave_iterative(self, iterations, threshold = 10, inspect=False):
        """Iterative method that fixes qubits that are thought to be found in the best position.
        
        Args:
            iterations (int): number of repetitions for the ancilla fixing.
            threshold (int): number of low energy solutions to check for same value.

        Returns:
            result (list): reconstructed best solution found after iterating.
            best_energy (float): energy of the system using the result output.
            w (int): iteration when the solution is first found.
            
        """
        bits = self.bits
        H = self.H
        sym = list(self.sym)

        fix = []
        out = []

        for w in range(iterations):

            response = self.run_dwave(H, self.sym)

            if inspect:
                insp.show(response)

            record = response.record
            order = np.argsort(record['energy'])

            best_sample, best_energy = self.get_energies(response)

            if best_energy == 0:
                print(f'Solution found! Final iteration: {w}\n')
                break
            
            print('Beginning fixing variables with the same value:\n')
            c = bits
            for j in range(bits, len(self.sym)):
                if j in fix:
                    continue
                else:
                    a = True
                    b = record.sample[order[0]][c]
                    for k in range(min(len(record.sample), threshold)):
                        if b != record.sample[order[k]][c]:
                            a = False
                            break
                    if a:
                        fix.append(j)
                        out.append(b)
                    c += 1
            print(f'The same value was found in positions {fix} \n')
            print(f'with values {out}.\n')

            for j in range(len(fix)):
                H = H.subs(self.sym[fix[j]], out[j])
                #sym.pop(fix[j])

            print(f'Total number of qubits needed for the next step: {len(self.sym)-len(fix)}.\n')

        print('Reconstructing state...\n')
        result = []
        c = 0
        for i in range(len(self.sym)):
            if i in fix:
                result.append(out[fix.index(i)])
            else:
                result.append(best_sample[c])
                c += 1

        print(f'Reconstructed result:\n')
        print(f'Relevant bits: {result[:bits]}\n')
        print(f'Ancillas: {result[bits:]}\n')

        best_energy = self.H.subs(((self.sym[i], result[i]) for i in range(len(self.sym))))
        print(f'With total energy: {best_energy}\n')

        return result, best_energy, w
