# Adiabatic quantum decoding of binary codes - Code

This repository contains code to solve instances of 3-7 and 8-16 bit nonlinear codes using D-Wave devices as part of the paper **Adiabatic quantum decoding of binary codes**. This code is intended to run with an existing D-Wave account or through the IDE Workspaces provided by D-Wave Leap (https://cloud.dwavesys.com/leap/).

The code is prepared to decode two examples. An instance of a 3-7 nonlinear code, and an instance of a 8-16 nonlinear code. The specifications of both codes are detailen underneath:

### 3-8 nonlinear code example:

Defining polynomials (ANF):

*x* = (*x*<sub>1</sub>, *x*<sub>3</sub>, *x*<sub>3</sub>) → *F* = (*f*<sub>1</sub>, ..., *f*<sub>8</sub>)

*f*<sub>1</sub> : *x*<sub>1</sub>*x*<sub>2</sub>*x*<sub>3</sub> + *x*<sub>1</sub>*x*<sub>2</sub> + *x*<sub>1</sub>*x*<sub>3</sub> + *x*<sub>1</sub> + *x*<sub>2</sub>*x*<sub>3</sub> + *x*<sub>2</sub>

### 8-16 nonlinear code example:

Defining polynomials (ANF):

*x* = (*x*<sub>1</sub>, ..., *x*<sub>8</sub>) → *F* = (*f*<sub>1</sub>, ..., *f*<sub>16</sub>)

*f*<sub>1</sub> : *x*<sub>1</sub> + *x*<sub>2</sub>

*f*<sub>2</sub> : *x*<sub>2</sub> + *x*<sub>3</sub>

*f*<sub>3</sub> : *x*<sub>3</sub> + *x*<sub>4</sub>

*f*<sub>4</sub> : *x*<sub>4</sub> + *x*<sub>5</sub>

*f*<sub>5</sub> : *x*<sub>5</sub> + *x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> + 1

*f*<sub>6</sub> : *x*<sub>2</sub> + *x*<sub>4</sub> + *x*<sub>5</sub> + *x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>6</sub>*x*<sub>8</sub> + *x*<sub>7</sub> + *x*<sub>8</sub>

*f*<sub>7</sub> : *x*<sub>1</sub> + *x*<sub>2</sub> + *x*<sub>4</sub> + *x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>6</sub>*x*<sub>8</sub> + *x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>7</sub> + *x*<sub>8</sub> + 1

*f*<sub>8</sub> : *x*<sub>1</sub> + *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>4</sub> + *x*<sub>5</sub> + *x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>6</sub>*x*<sub>7</sub> + *x*<sub>6</sub>*x*<sub>8</sub> + *x*<sub>7</sub>*x*<sub>8</sub>

*f*<sub>9</sub> : *x*<sub>2</sub> + *x*<sub>5</sub> + *x*<sub>6</sub>*x*<sub>7</sub> + *x*<sub>8</sub>

*f*<sub>10</sub> : *x*<sub>1</sub> + *x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>6</sub>*x*<sub>8</sub> + *x*<sub>6</sub> + *x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>7</sub> + *x*<sub>8</sub>

*f*<sub>11</sub> : *x*<sub>1</sub> + *x*<sub>4</sub> + *x*<sub>6</sub>*x*<sub>7</sub> + *x*<sub>6</sub>*x*<sub>8</sub> + *x*<sub>6</sub> + *x*<sub>7</sub>

*f*<sub>12</sub> : *x*<sub>1</sub> + *x*<sub>3</sub> + *x*<sub>5</sub> + *x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>6</sub>*x*<sub>8</sub> + *x*<sub>6</sub> + *x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>8</sub> + 1

*f*<sub>13</sub> : *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>4</sub> + *x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>6</sub>*x*<sub>8</sub> + *x*<sub>6</sub> + *x*<sub>7</sub>*x*<sub>8</sub>

*f*<sub>14</sub> : *x*<sub>1</sub> + *x*<sub>3</sub> + *x*<sub>4</sub> + *x*<sub>5</sub> + *x*<sub>6</sub>*x*<sub>8</sub> + *x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>7</sub> + 1

*f*<sub>15</sub> : *x*<sub>1</sub> + *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>5</sub> + *x*<sub>6</sub> + *x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>8</sub> + 1

*f*<sub>16</sub> : *x*<sub>3</sub> + *x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>6</sub>*x*<sub>7</sub> + *x*<sub>6</sub> + *x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>7</sub> + 1

*x* = (0, 1, 0, 0, 0, 1, 0, 1)

*F*(*x*) = (1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0)

*e* = (0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0)

*w*(*e*) = 2

*y* = *F*(*x*) + *e* = (1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0)

This examples have to be modified in order to be run on an adiabatic quantum computer. Detailed in the original paper are the procedures in more detail, and we present an outline below.

First the change <img src="https://render.githubusercontent.com/render/math?math=(x_i + x_j)\mod{2} \longrightarrow x_i + x_j - 2 x_i\cdot x_j"> is applied to the original ANF
function, to transform it into NNF. This is needed as Hamiltonians do not work on the binary basis. Then the problem Hamiltonian
*H*<sub>*p*</sub> is
constructed via the penalty <img src="https://render.githubusercontent.com/render/math?math=H_p = \sum_{i=1}^m\left(f_i(x)-r_i\right)^2">.

The transformation into NNF introduces several multi-qubit terms. All possible combination of products of
*x*<sub>*n*</sub> can appear in each function *f*<sub>*n*</sub>.

In order to reduce *H*<sub>*p*</sub> into 2-body interactions using a
minimal amount of ancillary qubits, a special mapping inspired by the following graph
is proposed.

<img src="https://github.com/Quantum-TII/hamming_codes/blob/master/ancillagraph.png" width="400">

Qubits are divided equally into two groups, within each they are fully
connected. Each vertex represents a qubit and each edge represents an
ancilla.

If *H*<sub>*p*</sub> has interaction up to *n* qubits, it would need up
to <img src="https://render.githubusercontent.com/render/math?math=2^{\frac{n %2B 2}{2}}-2"> total qubits to represent all possible
combination in each graph if *n* is even
(<img src="https://render.githubusercontent.com/render/math?math=3\times 2^{\frac{n-1}{2}}-2"> if *n* is odd).

The ancilla substitution scheme proposed for the 8-16 nonlinear code example is:

*x*<sub>1</sub>*x*<sub>2</sub> → *x*<sub>12</sub>

*x*<sub>3</sub>*x*<sub>4</sub> → *x*<sub>34</sub>

*x*<sub>1</sub>*x*<sub>3</sub> → *x*<sub>13</sub>

*x*<sub>2</sub>*x*<sub>4</sub> → *x*<sub>24</sub>

*x*<sub>1</sub>*x*<sub>4</sub> → *x*<sub>14</sub>

*x*<sub>2</sub>*x*<sub>3</sub> → *x*<sub>23</sub>

*x*<sub>5</sub>*x*<sub>6</sub> → *x*<sub>56</sub>

*x*<sub>7</sub>*x*<sub>8</sub> → *x*<sub>78</sub>

*x*<sub>5</sub>*x*<sub>7</sub> → *x*<sub>57</sub>

*x*<sub>6</sub>*x*<sub>8</sub> → *x*<sub>68</sub>

*x*<sub>5</sub>*x*<sub>8</sub> → *x*<sub>58</sub>

*x*<sub>6</sub>*x*<sub>7</sub> → *x*<sub>67</sub>

*x*<sub>1</sub>*x*<sub>34</sub> → *x*<sub>134</sub>

*x*<sub>2</sub>*x*<sub>34</sub> → *x*<sub>234</sub>

*x*<sub>3</sub>*x*<sub>12</sub> → *x*<sub>123</sub>

*x*<sub>4</sub>*x*<sub>12</sub> → *x*<sub>124</sub>

*x*<sub>5</sub>*x*<sub>78</sub> → *x*<sub>578</sub>

*x*<sub>6</sub>*x*<sub>78</sub> → *x*<sub>678</sub>

*x*<sub>7</sub>*x*<sub>56</sub> → *x*<sub>567</sub>

*x*<sub>8</sub>*x*<sub>56</sub> → *x*<sub>568</sub>

*x*<sub>12</sub>*x*<sub>34</sub> → *x*<sub>1234</sub>

*x*<sub>56</sub>*x*<sub>78</sub> → *x*<sub>5678</sub>

The maximum number of ancillary qubits is needed only if we have one
term in the ANF function that includes all members
*x*<sub>1</sub>, ..., *x*<sub>*n*</sub>, which is not always the case.

## Running the code

The code can be run with different options via the `main.py` file. The different arguments include:

- `--bits` (int): sets up the example to use in the program. It supports `3` and `8` bits of the examples detailed above.
- `--T` (float): determines the total annealing time. It defaults to 20 microseconds.
- `--chainstrength` (float): set a fixed *chainstrength* paramter for the D-Wave implementation. Leave blank for the automatic calculation based on the Hamiltonian.
- `--numruns` (int): number of samples for the quantum device. Default = `500`.
- `--control` (Bool): call in order to build the Hamiltonian minimizing the control issues.
- `--dwave` (Bool): call to run the example on dwave. Leave blank to only contruct the Hamiltonian.
- `--inspect` (Bool): open the D-Wave inspector after finishing the run.
- `--repeat` (Bool): perform the iterative approach outlined in the paper.
- `--iterations` (int): maximum number of iterations for the iterative approach. Deafault = `5`.
- `--rep` (int): number of repetitions of the entire algorithm.

The code returns a detailed analysis of the steps the algorithm takes in order to decode the binary code. At the end the code returns the solution that the quantum device has found to be correct.


