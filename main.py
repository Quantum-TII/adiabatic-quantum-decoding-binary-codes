import numpy as np
from sympy import symbols
from qibo import matrices, hamiltonians, models, callbacks
import functions
import argparse


def main(bits, T, dt):
    """

    Args:
        

    Returns:
        
    """
    if bits == 3:
        x1, x2, x3, xa1 = symbols('x1 x2 x3 xa1')
        z1, z2, z3, za1 = symbols('z1 z2 z3 za1')
        sym = [x1, x2, x3, xa1]
        ancillas = [xa1]
        sym_z = [z1, z2, z3, za1]  
        symbol_map = {z1: (0, matrices.Z), z2: (1, matrices.Z), z3: (2, matrices.Z), za1: (3, matrices.Z)}
        f = functions.f_3(x1, x2, x3)
        r = functions.r_3()
    elif bits == 8:
        x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
        xa1, xa2, xa3, xa4, xa5, xa6, xa7, xa8, xa9, xa10, xa11, xa12, xa13, xa14, xa15, xa16, xa17, xa18, xa19, xa20 = symbols('xa1 xa2 xa3 xa4 xa5 xa6 xa7 xa8 xa9 xa10 xa11 xa12 xa13 xa14 xa15 xa16 xa17 xa18 xa19 xa20')
        xa21, xa22, xa23, xa24, xa25, xa26, xa27, xa28, xa29, xa30, xa31, xa32, xa33, xa34, xa35, xa36, xa37, xa38, xa39, xa40 = symbols('xa21 xa22 xa23 xa24 xa25 xa26 xa27 xa28 xa29 xa30 xa31 xa32 xa33 xa34 xa35 xa36 xa37 xa38 xa39 xa40')
        xa41, xa42, xa43, xa44, xa45, xa46, xa47, xa48, xa49, xa50, xa51, xa52, xa53, xa54, xa55, xa56, xa57, xa58, xa59, xa60 = symbols('xa41 xa42 xa43 xa44 xa45 xa46 xa47 xa48 xa49 xa50 xa51 xa52 xa53 xa54 xa55 xa56 xa57 xa58 xa59 xa60')
        xa61, xa62, xa63, xa64, xa65, xa66, xa67, xa68, xa69, xa70, xa71, xa72, xa73, xa74, xa75, xa76, xa77, xa78, xa79, xa80 = symbols('xa61 xa62 xa63 xa64 xa65 xa66 xa67 xa68 xa69 xa70 xa71 xa72 xa73 xa74 xa75 xa76 xa77 xa78 xa79 xa80')
        z1, z2, z3, z4, z5, z6, z7, z8 = symbols('z1 z2 z3 z4 z5 z6 z7 z8')
        za1, za2, za3, za4, za5, za6, za7, za8, za9, za10, za11, za12, za13, za14, za15, za16, za17, za18, za19, za20 = symbols('za1 za2 za3 za4 za5 za6 za7 za8 za9 za10 za11 za12 za13 za14 za15 za16 za17 za18 za19 za20')
        za21, za22, za23, za24, za25, za26, za27, za28, za29, za30, za31, za32, za33, za34, za35, za36, za37, za38, za39, za40 = symbols('za21 za22 za23 za24 za25 za26 za27 za28 za29 za30 za31 za32 za33 za34 za35 za36 za37 za38 za39 za40')
        za41, za42, za43, za44, za45, za46, za47, za48, za49, za50, za51, za52, za53, za54, za55, za56, za57, za58, za59, za60 = symbols('za41 za42 za43 za44 za45 za46 za47 za48 za49 za50 za51 za52 za53 za54 za55 za56 za57 za58 za59 za60')
        za61, za62, za63, za64, za65, za66, za67, za68, za69, za70, za71, za72, za73, za74, za75, za76, za77, za78, za79, za80 = symbols('za61 za62 za63 za64 za65 za66 za67 za68 za69 za70 za71 za72 za73 za74 za75 za76 za77 za78 za79 za80')
    else:
        raise ValueError('Only instaces for 3 and 8 bits are supported')
        
        sym = [x1, x2, x3, x4, x5, x6, x7, x8, xa1, xa2, xa3, xa4, xa5, xa6, xa7, xa8, xa9, xa10, xa11, xa12, xa13, xa14, xa15, xa16, xa17, xa18, xa19, xa20, xa21, xa22, xa23, xa24, xa25, xa26, xa27, xa28, xa29, xa30, xa31, xa32, xa33, xa34, xa35, xa36, xa37, xa38, xa39, xa40, xa41, xa42, xa43, xa44, xa45, xa46, xa47, xa48, xa49, xa50, xa51, xa52, xa53, xa54, xa55, xa56, xa57, xa58, xa59, xa60, xa61, xa62, xa63, xa64, xa65, xa66, xa67, xa68, xa69, xa70, xa71, xa72, xa73, xa74, xa75, xa76, xa77, xa78, xa79, xa80]
        ancillas = [xa1, xa2, xa3, xa4, xa5, xa6, xa7, xa8, xa9, xa10, xa11, xa12, xa13, xa14, xa15, xa16, xa17, xa18, xa19, xa20, xa21, xa22, xa23, xa24, xa25, xa26, xa27, xa28, xa29, xa30, xa31, xa32, xa33, xa34, xa35, xa36, xa37, xa38, xa39, xa40, xa41, xa42, xa43, xa44, xa45, xa46, xa47, xa48, xa49, xa50, xa51, xa52, xa53, xa54, xa55, xa56, xa57, xa58, xa59, xa60, xa61, xa62, xa63, xa64, xa65, xa66, xa67, xa68, xa69, xa70, xa71, xa72, xa73, xa74, xa75, xa76, xa77, xa78, xa79, xa80]
        sym_z = [z1, z2, z3, z4, z5, z6, z7, z8, za1, za2, za3, za4, za5, za6, za7, za8, za9, za10, za11, za12, za13, za14, za15, za16, za17, za18, za19, za20, za21, za22, za23, za24, za25, za26, za27, za28, za29, za30, za31, za32, za33, za34, za35, za36, za37, za38, za39, za40, za41, za42, za43, za44, za45, za46, za47, za48, za49, za50, za51, za52, za53, za54, za55, za56, za57, za58, za59, za60, za61, za62, za63, za64, za65, za66, za67, za68, za69, za70, za71, za72, za73, za74, za75, za76, za77, za78, za79, za80]
        f = functions.f_8(x1, x2, x3, x4, x5, x6, x7, x8)
        r = functions.r_8()
        
    h = functions.to_hamiltonian(f, r, sym)
    #print(h)
    #print()
    h2, anc = functions.to_two_body(h, ancillas, bits)
    #print(h2)
    #print()
    hz = functions.to_pauli(h2, sym, sym_z)
    #print(hz)
    #print()
    #h_p = functions.to_matrix(hz, sym_z)
    hp_trotter = functions.symbolic_to_trotter(hz, symbol_map)
    hp_parts = hp_trotter.parts
    parts0 = functions.h0_parts(hp_parts)
    gs = lambda: functions.ground_state(nqubits)
    h0_trotter = hamiltonians.TrotterHamiltonian(*parts0, ground_state=gs)
    nqubits = bits+anc
    s = lambda t: t
    ground = callbacks.Gap(0)
    excited = callbacks.Gap(1)
    gap = callbacks.Gap()
    evolve_trotter = models.AdiabaticEvolution(h0_trotter, hp_trotter, s, dt, solver='exp',
                                       callbacks=[gap, ground, excited])
    initial_state = np.ones(2 ** nqubits) / np.sqrt(2 ** nqubits)
    final_state = evolve_trotter(final_time=T, initial_state=initial_state)
    output_dec = (np.abs(final_state.numpy())**2).argmax()
    max_output = "{0:0{bits}b}".format(output_dec, bits = nqubits)
    max_prob = (np.abs(final_state.numpy())**2).max()
    print(f'Adiabatic evolution with total time {T}, evolution step {dt}.\n')
    print(f'Most common solution after adiabatic evolution: {max_output[:nqubits-1]}.\n')
    print(f'Value of the ancillary qubits: {max_output[nqubits-1:]}.\n')
    print(f'Found with probability: {max_prob}.\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", default=3, type=int)
    parser.add_argument("--T", default=20, type=float)
    parser.add_argument("--dt", default=0.1, type=float)
    args = vars(parser.parse_args())
    main(**args)
