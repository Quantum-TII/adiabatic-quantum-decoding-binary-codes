from embedding_classes import direct_embedding
from dwave_classes import dwave_annealing
from examples_nnf import *
from sympy import symbols
import numpy as np
import argparse
import csv


def savedata(data, title):
    """Save data from a multidimensional list into csv files.

    Args:
        data (list): data you wish to save.
        title (str): name of the file to create.

    Returns:
        csv files for every dimension of the data.

    """
    for i in range(len(data)):
        Data=list(data[i])
        with open(f'{title}_{i}.csv', 'a') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(zip(Data))


def main(bits, it, th,ins):
    """Main program that collects measurements for a D-Wave annealing run for the
    MQ problem.

    Args:
        bits (int): number of bits of the instance to solve.
        it (int): maximum number of iterations for the iterative method.
        th (int): number of low energy solutions to check for the iterative method.

    Returns:
        csv files with the histogram data of the performed measurements.

    """
    x = symbols(" ".join((f"x{i}" for i in range(1, bits+1))))

    if bits == 5:
        p, sol = example_5(*x)
    elif bits == 9:
        p, sol = example_9(*x)
    else:
        raise ValueError(f'Not implemented for {bits} variables.')

    print(f'Actual solution(s) for the MQ problem: {sol}\n')

    direct = direct_embedding(bits, p, x)
    dir_H = direct.create_hamiltonian()
    dir_sym = direct.get_symbols()

    dir_dwave = dwave_annealing(dir_H, bits, dir_sym, numruns=1000, chainstrength=None)

    dir_result, dir_best_energy, dir_it = dir_dwave.dwave_iterative(iterations=it, threshold=th, inspect=ins)
    dir_en_hist = dir_dwave.energies_histogram

    savedata(dir_en_hist, title=f'direct_{bits}_bits_energies')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", default=5, type=int)
    parser.add_argument("--it", default=5, type=int)
    parser.add_argument("--th", default=10, type=int)
    parser.add_argument("--ins", action="store_true")
    args = vars(parser.parse_args())
    main(**args)