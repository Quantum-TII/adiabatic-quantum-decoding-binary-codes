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


def main(qmax, tru, pen):
    """Main program that collects the number of qubits needed to embed a Hamiltonian
    for solving the MQ problem into a D-Wave device.

    Args:
        qmax (int): maximum number of bits to collect data from (even numbers from 4 to 18)
        dire (Bool): set to True to use the direct embedding.
        pen (Bool): set to True to use the penalization embedding.

    Returns:
        csv files with the number of logical and physical qubits needed.

    """
    ps = []
    xs = []
    bit_list = [i for i in range(4, qmax+1, 2)]

    for bits in bit_list:
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
        ps.append(p)
        xs.append(x)

    if tru:
        tru_qubits = [[], [], []]
        for i in range(len(bit_list)):
            truncated = truncated_embedding(bit_list[i], ps[i], xs[i])
            tru_H = truncated.create_hamiltonian()
            tru_sym = truncated.get_symbols()

            tru_dwave = dwave_annealing(tru_H, bit_list[i], tru_sym, numruns=1)

            log_q, phys_q = tru_dwave.get_phys_qubits(tru_H, tru_sym, average=4)

            tru_qubits[0].append(bit_list[i])
            tru_qubits[1].append(log_q)
            tru_qubits[2].append(phys_q)

        savedata(tru_qubits, title=f'truncated_num_qubits')

    if pen:
        pen_qubits = [[], [], []]
        for i in range(len(bit_list)):
            penalization = penalization_embedding(bit_list[i], ps[i], xs[i])
            pen_H = penalization.create_hamiltonian()
            pen_sym = penalization.get_symbols()

            pen_dwave = dwave_annealing(pen_H, bit_list[i], pen_sym, numruns=1)

            log_q, phys_q = pen_dwave.get_phys_qubits(pen_H, pen_sym, average=4)

            pen_qubits[0].append(bit_list[i])
            pen_qubits[1].append(log_q)
            pen_qubits[2].append(phys_q)

        savedata(pen_qubits, title=f'penalty_num_qubits')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--qmax", default=12, type=int)
    parser.add_argument("--tru", action="store_true")
    parser.add_argument("--pen", action="store_true")
    args = vars(parser.parse_args())
    main(**args)