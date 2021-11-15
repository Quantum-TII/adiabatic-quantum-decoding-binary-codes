# Solving systems of Boolean multivariate equations with quantum annealing - Code

This repository contains the supporting code to solve Boolean Multivariate Quadratic (MQ) problems using D-Wave devices as part of the paper **Solving systems of Boolean multivariate equations with quantum annealing**. This code is intended to run with an existing D-Wave account or through the IDE Workspaces provided by D-Wave Leap (https://cloud.dwavesys.com/leap/).

The code is prepared to solve various instances using the different encoding methods presented in the paper, that is the direct embedding, the truncated embedding and the penalty embedding.

The examples for the truncated and penalization embedding can be loaded through `examples.py` while the examples for the direct embedding (after the ANF to NNF conversion) are stored in `examples_nnf.py`. Both are loaded by their corresponding codes.

Classes that implement the different embedding from a general MQ problem into the ground state of a two-body Hamiltonian are contained in the `embedding_classes.py` file. The required functions to deploy a given Hamiltonian to a D-Wave annealer are in the `dwave_classes.py` file.


## Running the code

There are three main code files for different returns.

The file `main_direct.py` solves the instance for the **direct embedding**. It can be modified by the following arguments:

- `--bits` (int): sets up the example to use in the program. It supports `5` and `9` bits. Default = `5`
- `--it` (int): number of iterations to perform for the iterative method. Default = `5`
- `--th` (int): number of low energy solution to check for the iterative method. Default = `10`
- `--ins` (Bool): Show the problem inspector after the run. Default = `False`.

The code returns a detailed analysis of the steps the algorithm takes in order to solve the MQ problem. At the end the code returns the solution that the quantum device has found.

To get the energy and solution values for the **truncated embedding** and the **penalty embedding** the file `main_energy.py` is to be used. The returns can be customized through the following arguments:

- `--bits` (int): sets up the example to use in the program. It supports `(4, 5, 6, 8, 10, 12, 14, 16, 18)` bits. Default = `5`
- `--tru` (Bool): Run the instance with the truncated embedding. Default = `False`.
- `--pen` (Bool): Run the instance with the penalty embedding. Default = `False`.
- `--it` (int): number of iterations to perform for the iterative method. Default = `5`
- `--th` (int): number of low energy solution to check for the iterative method. Default = `10`
- `--ins` (Bool): Show the problem inspector after the run. Default = `False`.

The code returns a detailed analysis of the steps the algorithm takes in order to solve the MQ problem. At the end the code returns the solution that the quantum device has found for both embeddings if selected.

To recover information regarding the number of physical qubits needed to map the **truncated embedding** and the **penalty embedding** the file `main_qubits.py` is to be used. In this case, the following arguments are provided:

- `--qmax` (int): Up to how many qubits to check for the quantum resources needed to solve the problem. Default = `12`
- `--tru` (Bool): Run the instance with the truncated embedding. Default = `False`.
- `--pen` (Bool): Run the instance with the penalty embedding. Default = `False`.




