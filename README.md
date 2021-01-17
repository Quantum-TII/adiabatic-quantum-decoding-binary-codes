# Solving linear and nonlinear codes using adiabatic evolution

(currently in progress)

This repository contains code to solve instances of nonlinear codes using the simulator Qibo and allows it to be run on dwave devices.

In this section we use a 8-16 code example to explain the procedure from
the ANF polynomials to the result provided by D-Wave.

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

First the change <img src="https://render.githubusercontent.com/render/math?math=(x_i + x_j)\mod{2} \longrightarrow x_i + x_j - 2 x_i\cdot x_j"> is applied to the original ANF
function, to transform it into NNF. Then the problem Hamiltonian
*H*<sub>*p*</sub> = *H*<sub>1</sub> + ... + *H*<sub>16</sub> is
constructed via the penalty <img src="https://render.githubusercontent.com/render/math?math=H_p = \sum_{i=1}^m\left(f_i(x)-r_i\right)^2">.

For example,

<div class="center">

*f*<sub>10</sub>(*NNF*) =  − 2*x*<sub>1</sub>*x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> + 4*x*<sub>1</sub>*x*<sub>6</sub>*x*<sub>7</sub> + 2*x*<sub>1</sub>*x*<sub>6</sub>*x*<sub>8</sub> − 2*x*<sub>1</sub>*x*<sub>6</sub> + 2*x*<sub>1</sub>*x*<sub>7</sub>*x*<sub>8</sub> − 2*x*<sub>1</sub>*x*<sub>7</sub> − 2*x*<sub>1</sub>*x*<sub>8</sub> + *x*<sub>1</sub> + *x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> − 2*x*<sub>6</sub>*x*<sub>7</sub> − *x*<sub>6</sub>*x*<sub>8</sub> + *x*<sub>6</sub> − *x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>7</sub> + *x*<sub>8</sub>

</div>

<div class="center">

*H*<sub>10</sub> = 2*x*<sub>1</sub>*x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> − 4*x*<sub>1</sub>*x*<sub>6</sub>*x*<sub>7</sub> − 2*x*<sub>1</sub>*x*<sub>6</sub>*x*<sub>8</sub> + 2*x*<sub>1</sub>*x*<sub>6</sub> − 2*x*<sub>1</sub>*x*<sub>7</sub>*x*<sub>8</sub> + 2*x*<sub>1</sub>*x*<sub>7</sub> + 2*x*<sub>1</sub>*x*<sub>8</sub> − *x*1 − *x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> + 2*x*<sub>6</sub>*x*<sub>7</sub> + *x*<sub>6</sub>*x*<sub>8</sub> − *x*<sub>6</sub> + *x*<sub>7</sub>*x*<sub>8</sub> − *x*<sub>7</sub> − *x*<sub>8</sub> + 1

When the ANF function is developed into NNF, the recursive action of <img src="https://render.githubusercontent.com/render/math?math=(x_i + x_j)\mod{2} \longrightarrow x_i + x_j - 2 x_i\cdot x_j"> can introduce all possible combination of products of
*x*<sub>*n*</sub> in each function *f*<sub>*n*</sub>. Every possible
interaction of qubits *x*<sub>1</sub>, ..., *x*<sub>8</sub> is appeared
in *H*<sub>*p*</sub> since *f*<sub>8</sub> includes every terms.

In order to reduce *H*<sub>*p*</sub> into 2-body interaction using a
minimal amount of ancillary qubits, a special mapping inspired by graph
is proposed.

<img src="https://github.com/Quantum-TII/hamming_codes/blob/master/ancillagraph.png" width="400">

Qubits are divided equally into two groups, within each they are fully
connected. Each vertex represents a qubit and each edge represents an
ancilla.

For *n* qubits, if *H*<sub>*p*</sub> includes all possibilities for up
to 3-qubit interaction. A total number of <img src="https://render.githubusercontent.com/render/math?math=\frac{n(n %2B 2)}{4}"> qubits is
needed to represent the problem Hamiltonian if *n* is even
(<img src="https://render.githubusercontent.com/render/math?math=\frac{(n %2B 1)^2}{4}"> if *n* is odd).

If *H*<sub>*p*</sub> has interaction up to *n* qubits, it would need up
to <img src="https://render.githubusercontent.com/render/math?math=2^{\frac{n %2B 2}{2}}-2"> total qubits to represent all possible
combination in each graph if *n* is even
(<img src="https://render.githubusercontent.com/render/math?math=3\times 2^{\frac{n-1}{2}}-2"> if *n* is odd).

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

Here is a 16-32 linear code example:

Defining polynomials (ANF):

*x* = (*x*<sub>1</sub>, ..., *x*<sub>16</sub>) → *F* = (*f*<sub>1</sub>, ..., *f*<sub>32</sub>)

*f*<sub>1</sub> = *x*<sub>1</sub> + *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>6</sub> + *x*<sub>7</sub> + *x*<sub>10</sub> + *x*<sub>12</sub> + *x*<sub>16</sub>

*f*<sub>2</sub> = *x*<sub>1</sub> + *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>5</sub> + *x*<sub>7</sub> + *x*<sub>8</sub> + *x*<sub>10</sub> + *x*<sub>11</sub> + *x*<sub>16</sub>

*f*<sub>3</sub> = *x*<sub>3</sub> + *x*<sub>4</sub> + *x*<sub>5</sub> + *x*<sub>7</sub> + *x*<sub>8</sub> + *x*<sub>14</sub>

*f*<sub>4</sub> = *x*<sub>4</sub> + *x*<sub>5</sub> + *x*<sub>8</sub> + *x*<sub>11</sub> + *x*<sub>13</sub> + *x*<sub>15</sub> + *x*<sub>16</sub>

*f*<sub>5</sub> = *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>4</sub> + *x*<sub>7</sub> + *x*<sub>8</sub> + *x*<sub>12</sub> + *x*<sub>14</sub> + *x*<sub>16</sub>

*f*<sub>6</sub> = *x*<sub>5</sub> + *x*<sub>6</sub> + *x*<sub>7</sub> + *x*<sub>8</sub> + *x*<sub>11</sub> + *x*<sub>12</sub> + *x*<sub>14</sub> + *x*<sub>16</sub>

*f*<sub>7</sub> = *x*<sub>2</sub> + *x*<sub>4</sub> + *x*<sub>6</sub> + *x*<sub>7</sub> + *x*<sub>9</sub> + *x*<sub>10</sub> + *x*<sub>11</sub> + *x*<sub>12</sub> + *x*<sub>15</sub> + *x*<sub>16</sub>

*f*<sub>8</sub> = *x*<sub>5</sub> + *x*<sub>7</sub> + *x*<sub>8</sub> + *x*<sub>9</sub> + *x*<sub>11</sub> + *x*<sub>12</sub> + *x*<sub>15</sub>

*f*<sub>9</sub> = *x*<sub>1</sub> + *x*<sub>2</sub> + *x*<sub>5</sub> + *x*<sub>6</sub> + *x*<sub>9</sub> + *x*<sub>12</sub> + *x*<sub>14</sub>

*f*<sub>10</sub> = *x*<sub>1</sub> + *x*<sub>2</sub> + *x*<sub>4</sub> + *x*<sub>5</sub> + *x*<sub>7</sub> + *x*<sub>8</sub> + *x*<sub>10</sub> + *x*<sub>11</sub> + *x*<sub>12</sub> + *x*<sub>14</sub> + *x*<sub>15</sub> + *x*<sub>16</sub>

*f*<sub>11</sub> = *x*<sub>3</sub> + *x*<sub>4</sub> + *x*<sub>6</sub> + *x*<sub>8</sub> + *x*<sub>10</sub> + *x*<sub>12</sub> + *x*<sub>13</sub> + *x*<sub>15</sub> + *x*<sub>16</sub>

*f*<sub>12</sub> = *x*<sub>1</sub> + *x*<sub>5</sub> + *x*<sub>8</sub> + *x*<sub>11</sub> + *x*<sub>13</sub> + *x*<sub>16</sub>

*f*<sub>13</sub> = *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>6</sub> + *x*<sub>9</sub> + *x*<sub>12</sub> + *x*<sub>14</sub> + *x*<sub>15</sub>

*f*<sub>14</sub> = *x*<sub>1</sub> + *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>4</sub> + *x*<sub>5</sub> + *x*<sub>7</sub>

*f*<sub>15</sub> = *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>6</sub> + *x*<sub>8</sub> + *x*<sub>9</sub> + *x*<sub>12</sub> + *x*<sub>14</sub> + *x*<sub>16</sub>

*f*<sub>16</sub> = *x*<sub>2</sub> + *x*<sub>4</sub> + *x*<sub>5</sub> + *x*<sub>6</sub> + *x*<sub>7</sub> + *x*<sub>13</sub> + *x*<sub>15</sub> + *x*<sub>16</sub>

*f*<sub>17</sub> = *x*<sub>4</sub> + *x*<sub>8</sub> + *x*<sub>9</sub> + *x*<sub>11</sub> + *x*<sub>12</sub> + *x*<sub>13</sub> + *x*<sub>15</sub> + *x*<sub>16</sub>

*f*<sub>18</sub> = *x*<sub>5</sub> + *x*<sub>6</sub> + *x*<sub>7</sub> + *x*<sub>10</sub> + *x*<sub>12</sub> + *x*<sub>13</sub> + *x*<sub>16</sub>

*f*<sub>19</sub> = *x*<sub>2</sub> + *x*<sub>4</sub> + *x*<sub>5</sub> + *x*<sub>6</sub> + *x*<sub>7</sub> + *x*<sub>8</sub> + *x*<sub>11</sub> + *x*<sub>12</sub> + *x*<sub>14</sub> + *x*<sub>16</sub>

*f*<sub>20</sub> = *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>7</sub> + *x*<sub>8</sub> + *x*<sub>12</sub> + *x*<sub>13</sub> + *x*<sub>14</sub> + *x*<sub>16</sub>

*f*<sub>21</sub> = *x*<sub>3</sub> + *x*<sub>5</sub> + *x*<sub>7</sub> + *x*<sub>8</sub> + *x*<sub>10</sub> + *x*<sub>11</sub>

*f*<sub>22</sub> = *x*<sub>2</sub> + *x*<sub>4</sub> + *x*<sub>5</sub> + *x*<sub>6</sub> + *x*<sub>8</sub> + *x*<sub>11</sub> + *x*<sub>12</sub> + *x*<sub>14</sub> + *x*<sub>15</sub> + *x*<sub>16</sub>

*f*<sub>23</sub> = *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>4</sub> + *x*<sub>9</sub> + *x*<sub>10</sub> + *x*<sub>11</sub> + *x*<sub>13</sub> + *x*<sub>15</sub>

*f*<sub>24</sub> = *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>4</sub> + *x*<sub>6</sub> + *x*<sub>8</sub> + *x*<sub>12</sub> + *x*<sub>13</sub> + *x*<sub>14</sub>

*f*<sub>25</sub> = *x*<sub>1</sub> + *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>7</sub> + *x*<sub>8</sub> + *x*<sub>9</sub> + *x*<sub>10</sub> + *x*<sub>11</sub> + *x*<sub>12</sub> + *x*<sub>13</sub> + *x*<sub>15</sub> + *x*<sub>16</sub>

*f*<sub>26</sub> = *x*<sub>1</sub> + *x*<sub>2</sub> + *x*<sub>5</sub> + *x*<sub>6</sub> + *x*<sub>7</sub> + *x*<sub>8</sub> + *x*<sub>9</sub> + *x*<sub>10</sub> + *x*<sub>11</sub> + *x*<sub>14</sub> + *x*<sub>15</sub>

*f*<sub>27</sub> = *x*<sub>2</sub> + *x*<sub>4</sub> + *x*<sub>5</sub> + *x*<sub>6</sub> + *x*<sub>8</sub> + *x*<sub>9</sub> + *x*<sub>10</sub> + *x*<sub>11</sub> + *x*<sub>16</sub>

*f*<sub>28</sub> = *x*<sub>5</sub> + *x*<sub>6</sub> + *x*<sub>7</sub> + *x*<sub>10</sub> + *x*<sub>11</sub> + *x*<sub>12</sub> + *x*<sub>13</sub> + *x*<sub>14</sub>

*f*<sub>29</sub> = *x*<sub>1</sub> + *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>4</sub> + *x*<sub>5</sub> + *x*<sub>8</sub> + *x*<sub>11</sub> + *x*<sub>13</sub> + *x*<sub>14</sub> + *x*<sub>16</sub>

*f*<sub>30</sub> = *x*<sub>2</sub> + *x*<sub>5</sub> + *x*<sub>6</sub> + *x*<sub>9</sub> + *x*<sub>10</sub> + *x*<sub>11</sub> + *x*<sub>13</sub> + *x*<sub>15</sub>

*f*<sub>31</sub> = *x*<sub>1</sub> + *x*<sub>2</sub> + *x*<sub>3</sub> + *x*<sub>13</sub> + *x*<sub>14</sub> + *x*<sub>15</sub> + *x*<sub>16</sub>

*f*<sub>32</sub> = *x*<sub>1</sub> + *x*<sub>3</sub> + *x*<sub>6</sub> + *x*<sub>8</sub> + *x*<sub>9</sub> + *x*<sub>10</sub> + *x*<sub>12</sub> + *x*<sub>13</sub> + *x*<sub>14</sub> + *x*<sub>15</sub> + *x*<sub>16</sub>

*x* = (1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1)

*F*(*x*) = (0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0,
0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1)

*e* = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

*w*(*e*) = 2

*y* = *F*(*x*) + *e* = (0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1,
1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1)

In this example, the total number of *x* is 16, but the maximum number
of terms included in a function is *m* = 12 (for *f*<sub>10</sub> and
*f*<sub>25</sub>).

In order to solve any 16-32 code, <img src="https://render.githubusercontent.com/render/math?math=2^{\frac{16 %2B 2}{2}}-2=510"> qubits is
needed. But for this case, it can be reduced down to 401 qubits.

The total number of qubits scales exponentially with *m*, with an upper
limit <img src="https://render.githubusercontent.com/render/math?math=2^{\frac{n %2B 2}{2}}-2">.
