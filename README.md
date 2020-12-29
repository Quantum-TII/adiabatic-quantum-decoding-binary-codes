# Solving linear and nonlinear codes using adiabatic evolution

(currently in progress)

This repository contains code to solve instances of nonlinear codes using the simulator Qibo and allows it to be run on dwave devices.

As of right now, only the example with 3 bits can be properly run, as the one with 8 has a problem when scaling down the interaction to two body.

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

First the change Eq. [\[eq:ANFNNF\]][1] is applied to the original ANF
function, to transform it into NNF. Then the problem Hamiltonian
*H*<sub>*p*</sub> = *H*<sub>1</sub> + ... + *H*<sub>16</sub> is
constructed via the penalty Eq. [\[eq:penalty\]][2].

For example,

<div class="center">

*f*<sub>10</sub>(*N**N**F*) =  − 2*x*<sub>1</sub>*x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> + 4*x*<sub>1</sub>*x*<sub>6</sub>*x*<sub>7</sub> + 2*x*<sub>1</sub>*x*<sub>6</sub>*x*<sub>8</sub> − 2*x*<sub>1</sub>*x*<sub>6</sub> + 2*x*<sub>1</sub>*x*<sub>7</sub>*x*<sub>8</sub> − 2*x*<sub>1</sub>*x*<sub>7</sub> − 2*x*<sub>1</sub>*x*<sub>8</sub> + *x*<sub>1</sub> + *x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> − 2*x*<sub>6</sub>*x*<sub>7</sub> − *x*<sub>6</sub>*x*<sub>8</sub> + *x*<sub>6</sub> − *x*<sub>7</sub>*x*<sub>8</sub> + *x*<sub>7</sub> + *x*<sub>8</sub>

</div>

<div class="center">

*H*<sub>10</sub> = 2*x*<sub>1</sub>*x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> − 4*x*<sub>1</sub>*x*<sub>6</sub>*x*<sub>7</sub> − 2*x*<sub>1</sub>*x*<sub>6</sub>*x*<sub>8</sub> + 2*x*<sub>1</sub>*x*<sub>6</sub> − 2*x*<sub>1</sub>*x*<sub>7</sub>*x*<sub>8</sub> + 2*x*<sub>1</sub>*x*<sub>7</sub> + 2*x*<sub>1</sub>*x*<sub>8</sub> − *x*1 − *x*<sub>6</sub>*x*<sub>7</sub>*x*<sub>8</sub> + 2*x*<sub>6</sub>*x*<sub>7</sub> + *x*<sub>6</sub>*x*<sub>8</sub> − *x*<sub>6</sub> + *x*<sub>7</sub>*x*<sub>8</sub> − *x*<sub>7</sub> − *x*<sub>8</sub> + 1

</div>

  [1]: #eq:ANFNNF
  [2]: #eq:penalty
