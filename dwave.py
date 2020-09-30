import numpy as np
from sympy import symbols
import functions
import argparse
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite


def main(bits, chainstrength, nruns):
    """

    Args:
        

    Returns:
        
    """
    if bits == 3:
        x1, x2, x3, xa1 = symbols('x1 x2 x3 xa1')
        z1, z2, z3, za1 = symbols('z1 z2 z3 za1')
        sym = [x1, x2, x3, xa1]
        ancillas = [xa1]
        f = functions.f_3(x1, x2, x3)
        r = functions.r_3()
    elif bits == 8:
        x1, x2, x3, x4, x5, x6, x7, x8 = symbols('x1 x2 x3 x4 x5 x6 x7 x8')
        xa1, xa2, xa3, xa4, xa5, xa6, xa7, xa8, xa9, xa10, xa11, xa12, xa13, xa14, xa15, xa16, xa17, xa18, xa19, xa20 = symbols('xa1 xa2 xa3 xa4 xa5 xa6 xa7 xa8 xa9 xa10 xa11 xa12 xa13 xa14 xa15 xa16 xa17 xa18 xa19 xa20')
        xa21, xa22, xa23, xa24, xa25, xa26, xa27, xa28, xa29, xa30, xa31, xa32, xa33, xa34, xa35, xa36, xa37, xa38, xa39, xa40 = symbols('xa21 xa22 xa23 xa24 xa25 xa26 xa27 xa28 xa29 xa30 xa31 xa32 xa33 xa34 xa35 xa36 xa37 xa38 xa39 xa40')
        xa41, xa42, xa43, xa44, xa45, xa46, xa47, xa48, xa49, xa50, xa51, xa52, xa53, xa54, xa55, xa56, xa57, xa58, xa59, xa60 = symbols('xa41 xa42 xa43 xa44 xa45 xa46 xa47 xa48 xa49 xa50 xa51 xa52 xa53 xa54 xa55 xa56 xa57 xa58 xa59 xa60')
        xa61, xa62, xa63, xa64, xa65, xa66, xa67, xa68, xa69, xa70, xa71, xa72, xa73, xa74, xa75, xa76, xa77, xa78, xa79, xa80 = symbols('xa61 xa62 xa63 xa64 xa65 xa66 xa67 xa68 xa69 xa70 xa71 xa72 xa73 xa74 xa75 xa76 xa77 xa78 xa79 xa80')
        sym = [x1, x2, x3, x4, x5, x6, x7, x8, xa1, xa2, xa3, xa4, xa5, xa6, xa7, xa8, xa9, xa10, xa11, xa12, xa13, xa14, xa15, xa16, xa17, xa18, xa19, xa20, xa21, xa22, xa23, xa24, xa25, xa26, xa27, xa28, xa29, xa30, xa31, xa32, xa33, xa34, xa35, xa36, xa37, xa38, xa39, xa40, xa41, xa42, xa43, xa44, xa45, xa46, xa47, xa48, xa49, xa50, xa51, xa52, xa53, xa54, xa55, xa56, xa57, xa58, xa59, xa60, xa61, xa62, xa63, xa64, xa65, xa66, xa67, xa68, xa69, xa70, xa71, xa72, xa73, xa74, xa75, xa76, xa77, xa78, xa79, xa80]
        ancillas = [xa1, xa2, xa3, xa4, xa5, xa6, xa7, xa8, xa9, xa10, xa11, xa12, xa13, xa14, xa15, xa16, xa17, xa18, xa19, xa20, xa21, xa22, xa23, xa24, xa25, xa26, xa27, xa28, xa29, xa30, xa31, xa32, xa33, xa34, xa35, xa36, xa37, xa38, xa39, xa40, xa41, xa42, xa43, xa44, xa45, xa46, xa47, xa48, xa49, xa50, xa51, xa52, xa53, xa54, xa55, xa56, xa57, xa58, xa59, xa60, xa61, xa62, xa63, xa64, xa65, xa66, xa67, xa68, xa69, xa70, xa71, xa72, xa73, xa74, xa75, xa76, xa77, xa78, xa79, xa80]
        f = functions.f_8(x1, x2, x3, x4, x5, x6, x7, x8)
        r = functions.r_8()
        
    h = functions.to_hamiltonian(f, r, sym)
    #print(h)
    #print()
    h2, anc = functions.to_two_body(h, ancillas, bits)
    #print(h2)
    #print()
    symbol_num = sym_dict_num(sym)
    Q, constant = functions.symbolic_to_dwave(hz, symbol_num)


    # Run the QUBO on the solver from your config file
    sampler = EmbeddingComposite(DWaveSampler(solver={'qpu': True}))
    response = sampler.sample_qubo(Q, chain_strength=chainstrength, num_reads=nruns)

    print(response)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", default=3, type=int)
    parser.add_argument("--chainstrength", default=1, type=float)
    parser.add_argument("--nruns", default=10, type=int)
    args = vars(parser.parse_args())
    main(**args)