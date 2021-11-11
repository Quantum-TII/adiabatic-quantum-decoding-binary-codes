import numpy as np
from sympy import symbols, expand
import collections
import itertools
from scipy.special import binom


class truncated_embedding:
    """Class that implements a truncated embedding with cutoff = 4 of the
    solution of an MQ problem into the ground state of an up to 2-body 
    Hamiltonian

    """
    def __init__(self, bits, p, x, k=4):
        """Initialize the variables needed for the class.
        Args:
            bits (int): number of variables present in the MQ problem.
            p (list): equations to solve with sympy.Symbols.
            x (sympy.Symbols): symbols used in the equations p.
            k (int): cutoff for the transformation of equations (Only k=4 functional)

        """
        self.bits = bits
        self.p = p
        self.x = x
        self.a = symbols(' '.join((f'a{i}' for i in range(1, int(binom(bits, 2))+1))))
        self.H = 0
        self.k = k

    def gadget(self, x, y, a):
        """Penalization function for the gadget substitution (x*y -> a).
        
        Args:
            x, y (symbol): qubits to be substituted.
            a (symbol): ancilla for the substitution.
            
        Returns:
            penalization function.
        
        """
        return 3*a + x*y - 2*x*a - 2*y*a

    def reduce_to_one_body(self):
        """Add ancillary variables in order to only have single variable in the 
        equations. Add the necessary gadgets into the Hamiltonian.

        """
        x, a = self.x, self.a
        subs = [i for i in itertools.combinations(x, 2)]
        for i in range(len(a)):
            d = 0
            for j in range(len(self.p)):
                if subs[i][0]*subs[i][1] in self.p[j].args:
                    d += 1
                    self.p[j] = self.p[j].subs(subs[i][0]*subs[i][1], a[i])
            self.H += (1+d)*self.gadget(subs[i][0], subs[i][1], a[i])

    def cutoff_p(self):
        """Cutoff the original equations into up most k=4 variable terms.

        """
        k = self.k
        p = self.p
        l = 0
        for i in range(len(p)):
            l += int(np.ceil((len(p[i].args)-k)/(k-2)))
        b = symbols(' '.join((f'b{i}' for i in range(1, l+1))))
        self.b = b
        new_p = []
        z = 0
        for i in range(len(p)):
            if len(p[i].args) <= k:
                new_p.append(p[i])
                continue
            else:
                new_p.append(sum(p[i].args[:k-1])+b[z])
                z += 1
                for j in range(k-1, len(p[i].args), k-2):
                    if len(p[i].args[j:]) < k:
                        new_p.append(b[z-1]+sum(p[i].args[j:]))
                        break
                    else:
                        new_p.append(b[z-1]+sum(p[i].args[j:j+k-2])+b[z])
                        z += 1
        return new_p

    def encode_new_p(self, new_p):
        """Go from ANF into NNF and create the final Hamiltonian for the
        short equations.

        Args:
            new_p (list): list of equations of up to 4 variables. 
            (3 variables might appear by construction)

        Returns:
            self.H (Symbols): Hamiltonian that encodes the solution of the
            MQ problem into its ground state.

        """
        if self.k != 4:
            raise ValueError("Only cutoff of k=4 supported for the Hamiltonian construction.")
        l = 0
        for i in range(len(new_p)):
            if len(new_p[i].args) == 4:
                if new_p[i].args[0] == 1:
                    l += 1
                else:
                    l += 2
            elif len(new_p[i].args) == 3:
                if new_p[i].args[0] == 1:
                    l += 0
                else:
                    l += 1
        c = symbols(' '.join((f'c{i}' for i in range(1, l+1))))
        self.c = c

        z = 0
        for i in range(len(new_p)):
            if len(new_p[i].args) == 4:
                if new_p[i].args[0] == 1:
                    self.H += -new_p[i].args[1]-new_p[i].args[2]-new_p[i].args[3]+1+2*new_p[i].args[1]*new_p[i].args[2]+2*new_p[i].args[1]*new_p[i].args[3]+2*new_p[i].args[2]*new_p[i].args[3]-4*new_p[i].args[1]*c[z]
                    self.H += 5*self.gadget(new_p[i].args[1], new_p[i].args[2], c[z])
                    z += 1
                else:
                    self.H += new_p[i].args[2]+new_p[i].args[0]+new_p[i].args[1]+new_p[i].args[3]-2*new_p[i].args[2]*new_p[i].args[0]-2*new_p[i].args[2]*new_p[i].args[1]-2*new_p[i].args[2]*new_p[i].args[3]-2*new_p[i].args[0]*new_p[i].args[1]-2*new_p[i].args[0]*new_p[i].args[3]-2*new_p[i].args[1]*new_p[i].args[3]+4*new_p[i].args[2]*c[z]+4*c[z+1]*new_p[i].args[0]+4*c[z]*new_p[i].args[3]+4*c[z+1]*new_p[i].args[1]-8*c[z]*c[z+1]
                    self.H += 9*self.gadget(new_p[i].args[0], new_p[i].args[1], c[z])
                    self.H += 9*self.gadget(new_p[i].args[2], new_p[i].args[3], c[z+1])
                    z += 2
            elif len(new_p[i].args) == 3:
                if new_p[i].args[0] == 1:
                    self.H += -new_p[i].args[1]-new_p[i].args[2]+2*new_p[i].args[1]*new_p[i].args[2]+1
                else:
                    self.H += new_p[i].args[0]+new_p[i].args[1]+new_p[i].args[2]-2*new_p[i].args[0]*new_p[i].args[1]-2*new_p[i].args[0]*new_p[i].args[2]-2*new_p[i].args[1]*new_p[i].args[2]+4*new_p[i].args[0]*c[z]
                    self.H += 5*self.gadget(new_p[i].args[1], new_p[i].args[2], c[z])
                    z += 1
        return self.H

    def create_hamiltonian(self):
        """Combine all steps to create the full Hamiltonian.

        Returns:
            self.H (Symbols): Hamiltonian that encodes the solution of the
            MQ problem into its ground state.

        """
        self.reduce_to_one_body()
        new_p = self.cutoff_p()
        self.encode_new_p(new_p)
        return self.H

    def get_symbols(self):
        """Get the total number of qubits needed for the implementation.

        Returns:
            sym (tuple): all qubits used in the Hamiltonian embedding.

        """
        self.sym = self.x+self.a+self.b+self.c
        return self.sym

    def __call__(self):
        """Equivalent to `truncated_embedding.construct_hamiltonian`."""
        return self.create_hamiltonian()


class penalization_embedding:
    """Class that implements a penalization embedding inspired by a 
    circuit to adiabatic conversion to encode the solution of an 
    MQ problem into the ground state of an up to 2-body Hamiltonian

    """
    def __init__(self, bits, p, x, p_cnot=3, p_sol=2):
        """Initialize the variables needed for the class.
        Args:
            bits (int): number of variables present in the MQ problem.
            p (list): equations to solve with sympy.Symbols.
            x (sympy.Symbols): symbols used in the equations p.
            p_cnot (int): amount of penalization to the circuit gates.
            p_sol (int): ampunt of penalization to the solution.

        """
        self.bits = bits
        self.p = p
        self.x = x
        self.p_cnot = p_cnot
        self.p_sol = p_sol
        self.outputs = self.create_outputs()

    def cnot_penalization(self, xt, xc, xr, xa):
        """Penalization function for a CNOT gate. 
        
        Args:
            xt (symbol): target qubit.
            xc (symbol): control qubit.
            xr (symbol): qubit where the result of the CNOT is stored.
            xa (symbol): ancilla qubit to break 3-body interactions.
            
        Returns:
            Term penalizing combinations that do not represent a CNOT gate.
            
        """
        return self.p_cnot*(2*xc*xt - 2*xc*xr - 2*xt*xr - 4*xc*xa - 4*xt*xa + 4*xr*xa + xc + xt + xr + 4*xa)

    def and_penalization(self, x1, x2, x12):
        """Penalization function for an AND gate.
        
        Args:
            x1 (symbol): first qubit.
            x2 (symbol): second qubit.
            x12 (symbol): multiplication of the first and second qubits.
            
        Returns:
            Term penalizing combinations that do not represent an AND gate.
            
        """
        return self.p_cnot*(x1*x2 - 2*x1*x12 - 2*x2*x12 + 3*x12)

    def toffoli_penalization(self, xt, xc1, xc2, xr, xb, xa):
        """Penalization function for a Toffoli gate.
        
        Args:
            xt (symbol): target qubit.
            xc1 (symbol): control qubit.
            xc2 (symbol): control qubit.
            xr (symbol): qubit where the result of the Toffoli is stored.
            xb (symbol): multiplication of the two control qubits.
            xa (symbol): ancilla qubit to break 3-body interactions.
            
        Returns:
            Term penalizing combinations that do not represent a Toffoli gate.
            
        """
        pen = self.cnot_penalization(xt, xb, xr, xa)
        pen += self.and_penalization(xc1, xc2, xb)
        return pen

    def const_penalization(self, xt, xr):
        """Penalization function for an X gate.
        
        Args:
            xt (symbol): target qubit.
            xr (symbol): qubit after the X gate is performed.
            xa (symbol): ancilla qubit to break 3-body interactions.
            
        Returns:
            Term penalizing combinations that do not represent an X gate.
            
        """
        return self.p_cnot*(2*xt*xr - xt - xr + 1)

    def create_outputs(self):
        """Assign the variables that are going to be used as outputs throughout the computation. 
            
        Returns:
            outputs (list): extra qubits needed to store the results of applying gates throughout the computation.
            
        """
        outputs = []
        o = 1
        for i in range(len(self.p)):
            out = []
            xo = symbols(f'xo{o}')
            out.append(xo)
            o += 1
            for j in range(len(self.p[i].args)):
                xo = symbols(f'xo{o}')
                out.append(xo)
                o += 1
            outputs.append(out)
        print(f'Total output qubits used: {o}\n')
        return outputs

    def create_hamiltonian(self):
        """Creation of the Hamiltonian using a polynomal number of ancilla qubits with the number of terms in the function.
        
        Returns:
            H (symbol): hamiltonian that encodes in its ground state the solution of the problem.
            anc (list): ancillas used to create the hamiltonian.
            
        """
        outputs = self.outputs
        x = self.x
        H = 0
        anc = []
        ancillas = {}
        c = 1
        for i, row in enumerate(self.p):
            H += 2*self.p_cnot*outputs[i][0]
            for j, term in enumerate(row.args):
                if not term.args:
                    expression = (term,)
                else:
                    expression = term.args
                symbol = [x for x in expression if x.is_symbol]
                numbers = [x for x in expression if not x.is_symbol]

                if len(numbers) > 1:
                    raise ValueError("Expression must be expanded before using this method.")
                elif numbers:
                    constant = float(numbers[0])
                else:
                    constant = 1
                
                if len(symbol) == 0:
                    H += self.const_penalization(outputs[i][j], outputs[i][j+1])
                elif len(symbol) == 1:
                    if ancillas.get(symbol[0]*outputs[i][j]) is None:
                        a = symbols(f'xa{c}')
                        anc.append(a)
                        ancillas[symbol[0]*outputs[i][j]] = a
                        c += 1
                    H += self.cnot_penalization(outputs[i][j], symbol[0], outputs[i][j+1], ancillas[symbol[0]*outputs[i][j]])
                elif len(symbol) == 2:
                    if ancillas.get(symbol[0]*symbol[1]) is None:
                        a = symbols(f'xa{c}')
                        anc.append(a)
                        ancillas[symbol[0]*symbol[1]] = a
                        c += 1
                    if ancillas.get(ancillas[symbol[0]*symbol[1]]*outputs[i][j]) is None:
                        a = symbols(f'xa{c}')
                        anc.append(a)
                        ancillas[ancillas[symbol[0]*symbol[1]]*outputs[i][j]] = a
                        c += 1
                    H += self.toffoli_penalization(outputs[i][j], symbol[0], symbol[1], outputs[i][j+1], ancillas[symbol[0]*symbol[1]], ancillas[ancillas[symbol[0]*symbol[1]]*outputs[i][j]])
            H += self.p_sol*outputs[i][j+1]
        self.H = expand(H)
        self.anc = anc
        return self.H

    def get_symbols(self):
        """Get the total number of qubits needed for the implementation.

        Returns:
            sym (tuple): all qubits used in the Hamiltonian embedding.

        """
        sym = self.x
        for i in range(len(self.outputs)):
            sym += tuple(self.outputs[i])
        sym += tuple(self.anc)
        self.sym = sym
        return self.sym

    def __call__(self):
        """Equivalent to `penalization_embedding.create_hamiltonian`."""
        return self.create_hamiltonian()


class direct_embedding:
    """Class that implements a direct embedding of the solution of an MQ 
    problem into the ground state of an up to 2-body Hamiltonian

    """
    def __init__(self, bits, p, x, control=True):
        """Initialize the variables needed for the class.
        Args:
            bits (int): number of variables present in the MQ problem.
            p (list): equations to solve with sympy.Symbols.
            x (sympy.Symbols): symbols used in the equations p.
            control (Bool): whether minimum control is desired.

        """
        self.bits = bits
        self.p = p
        self.x = x
        if bits//2 == 0:
            self.a = symbols(' '.join((f'a{i}' for i in range(1, int(2**((bits+2)/2)-2-bits+1)))))
        else:
            self.a = symbols(' '.join((f'a{i}' for i in range(1,int( 3*2**((bits-1)/2)-2-bits+1)))))
        self.H = 0
        self.control = control
        self.sym = self.x+self.a

    def gadget(self, x, y, a):
        """Penalization function for the gadget substitution (x*y -> a).
        
        Args:
            x, y (symbol): qubits to be substituted.
            a (symbol): ancilla for the substitution.
            
        Returns:
            penalization function.
        
        """
        return 3*a + x*y - 2*x*a - 2*y*a

    def symbolic_to_data(self, symbolic_hamiltonian):
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

    def add_gadget(self, x1, x2, xa):
        """Add the gadget ancillas.
        
        Args:
            x1, x2 (symbol): qubits to substitute.
            xa (symbol): ancilla qubit introduced.

        """
        terms_s, terms_c, overall_constant = self.symbolic_to_data(self.H)
        self.H = expand(self.H.subs(x1*x2, xa))
        if self.control:
            d_m = 0
            d_p = 0
            for j in range(len(terms_s)):
                if x1 in terms_s[j] and x2 in terms_s[j]:
                    if terms_c[j] < 0:
                        d_m += -terms_c[j]
                    else:
                        d_p += terms_c[j]
            self.H += (1+max(d_m, d_p))*self.gadget(x1, x2, xa)
        else:
            for j in range(len(terms_s)):
                if x1 in terms_s[j] and x2 in terms_s[j]:
                    self.H += (1+np.abs(terms_c[j]))*gadget(x1, x2, xa)

    def to_gadget_5(self, x, ancillas):
        """Function that decomposes a multi-qubit interaction hamiltonian into 2-qubit using ancillas.
        
        Args:
            x (list): original variables.
            ancillas (list): ancillas to use.
            
        """
        self.add_gadget(x[0], x[1], ancillas[0])
        self.add_gadget(x[3], x[4], ancillas[1])
        self.add_gadget(x[2], x[5+1], ancillas[2])
        self.add_gadget(x[2], x[3], ancillas[3])
        self.add_gadget(x[2], x[4], ancillas[4])
        
    def to_gadget_9(self, x, ancillas):
        """Function that decomposes a multi-qubit interaction hamiltonian into 2-qubit using ancillas.
        
        Args:
            x (list): all variables.
            ancillas (list): ancillas to use.
            
        """
        self.add_gadget(x[0], x[1], ancillas[0])
        self.add_gadget(x[2], x[3], ancillas[1])
        self.add_gadget(x[9+0], x[9+1], ancillas[2]) #01 23
        self.add_gadget(x[4], x[9+2], ancillas[3]) #0123 4
        self.add_gadget(x[2], x[4], ancillas[4])
        self.add_gadget(x[9+0], x[9+4], ancillas[5]) #01 24
        self.add_gadget(x[3], x[4], ancillas[6])
        self.add_gadget(x[9+0], x[9+6], ancillas[7]) #01 34
        self.add_gadget(x[2], x[9+0], ancillas[8]) #2 01
        self.add_gadget(x[3], x[9+0], ancillas[9]) #3 01
        self.add_gadget(x[4], x[9+0], ancillas[10]) #4 01

        self.add_gadget(x[0], x[4], ancillas[11])
        self.add_gadget(x[9+1], x[9+11], ancillas[12]) #04 23
        self.add_gadget(x[1], x[4], ancillas[13])
        self.add_gadget(x[9+1], x[9+13], ancillas[14]) #14 23
        self.add_gadget(x[0], x[9+1], ancillas[15]) #0 23
        self.add_gadget(x[1], x[9+1], ancillas[16]) #1 23
        self.add_gadget(x[4], x[9+1], ancillas[17]) #4 23

        self.add_gadget(x[1], x[9+4], ancillas[18]) #1 24
        self.add_gadget(x[1], x[9+6], ancillas[19]) #1 34
        self.add_gadget(x[0], x[9+6], ancillas[20]) #0 34
        self.add_gadget(x[0], x[9+4], ancillas[21]) #0 24
        self.add_gadget(x[0], x[3], ancillas[22])
        self.add_gadget(x[1], x[3], ancillas[23])
        self.add_gadget(x[0], x[2], ancillas[24])
        self.add_gadget(x[1], x[2], ancillas[25])

        self.add_gadget(x[5], x[6], ancillas[26])
        self.add_gadget(x[7], x[8], ancillas[27])
        self.add_gadget(x[9+26], x[9+27], ancillas[28]) #56 78
        self.add_gadget(x[7], x[9+26], ancillas[29]) #7 56
        self.add_gadget(x[8], x[9+26], ancillas[30]) #8 56
        self.add_gadget(x[5], x[9+27], ancillas[31]) #5 78
        self.add_gadget(x[6], x[9+27], ancillas[32]) #6 78
        self.add_gadget(x[5], x[7], ancillas[33])
        self.add_gadget(x[5], x[8], ancillas[34])
        self.add_gadget(x[6], x[7], ancillas[35])
        self.add_gadget(x[6], x[8], ancillas[36]) 

    def create_hamiltonian(self):
        """Combine all steps to create the full Hamiltonian.

        Returns:
            self.H (Symbols): Hamiltonian that encodes the solution of the
            MQ problem into its ground state.

        """
        self.H += np.sum(self.p)
        self.H = expand(self.H)
        if self.bits == 5:
            self.to_gadget_5(self.sym, self.a)
        elif self.bits == 9:
            self.to_gadget_9(self.sym, self.a)
        else:
            raise ValueError(f"Not implemented for {self.bits} variables.")
        return expand(self.H)

    def get_symbols(self):
        """Get the total number of qubits needed for the implementation.

        Returns:
            sym (tuple): all qubits used in the Hamiltonian embedding.

        """
        return self.sym

    def __call__(self):
        """Equivalent to `direct_embedding.create_hamiltonian`."""
        return self.create_hamiltonian()
