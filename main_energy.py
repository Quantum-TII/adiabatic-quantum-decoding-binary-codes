from embedding_classes import truncated_embedding, penalization_embedding
from dwave_classes import dwave_annealing
from examples import *
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


def main(bits, tru, pen, it, th,ins):
    """Main program that collects measurements for a D-Wave annealing run for the
    MQ problem.

    Args:
        bits (int): number of bits of the instance to solve.
        tru (Bool): set to True to use the truncated embedding.
        pen (Bool): set to True to use the penalization embedding.
        it (int): maximum number of iterations for the iterative method.
        th (int): number of low energy solutions to check for the iterative method.

    Returns:
        csv files with the histogram data of the performed measurements.

    """
    x = symbols(" ".join((f"x{i}" for i in range(1, bits+1))))

    if bits == 4:
        p, sol = example_4(*x)
    if bits == 5:
        p, sol = example_5(*x)
    elif bits == 6:
        p, sol = example_6(*x)
    elif bits == 8:
        p, sol = example_8(*x)
    elif bits == 10:
        p, sol = example_10(*x)
    elif bits == 12:
        p, sol = example_12(*x)
    elif bits == 14:
        p, sol = example_14(*x)
    elif bits == 16:
        p, sol = example_16(*x)
    elif bits == 18:
        p, sol = example_18(*x)

    print(f'Actual solution(s) for the MQ problem: {sol}\n')

    if tru:
        truncated = truncated_embedding(bits, p, x)
        tru_H = truncated.create_hamiltonian()
        tru_sym = truncated.get_symbols()

        tru_dwave = dwave_annealing(tru_H, bits, tru_sym, numruns=1000, chainstrength=None)

        tru_result, tru_best_energy, tru_it = tru_dwave.dwave_iterative(iterations=it, threshold=th, inspect=ins)
        tru_en_hist = tru_dwave.energies_histogram

        savedata(tru_en_hist, title=f'truncated_{bits}_bits_energies')

    if pen:
        penalization = penalization_embedding(bits, p, x, p_cnot=3, p_sol=2)
        pen_H = penalization.create_hamiltonian()
        pen_sym = penalization.get_symbols()
        

        pen_dwave = dwave_annealing(pen_H, bits, pen_sym, numruns=1000, chainstrength=None)

        pen_result, pen_best_energy, pen_it = pen_dwave.dwave_iterative(iterations=it, threshold=th, inspect=ins)
        pen_en_hist = pen_dwave.energies_histogram

        savedata(pen_en_hist, title=f'penalty_{bits}_bits_energies')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", default=5, type=int)
    parser.add_argument("--tru", action="store_true")
    parser.add_argument("--pen", action="store_true")
    parser.add_argument("--it", default=5, type=int)
    parser.add_argument("--th", default=10, type=int)
    parser.add_argument("--ins", action="store_true")
    args = vars(parser.parse_args())
    main(**args)