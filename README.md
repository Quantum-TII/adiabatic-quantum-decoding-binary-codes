# Solving linear and nonlinear codes using adiabatic evolution

(currently in progress)

This repository contains code to solve instances of nonlinear codes using the simulator Qibo and allows it to be run on dwave devices.

As of right now, only the example with 3 bits can be properly run, as the one with 8 has a problem when scaling down the interaction to two body.

In this section we use a 8-16 code example to explain the procedure from the ANF polynomials to the result provided by D-Wave.

Defining polynomials (ANF): 



x = (x_1,...,x_8) \longrightarrow F =(f_1,...,f_{16})

f_{1}: x_1 + x_2

f_{2}: x_2 + x_3

f_{3}: x_3 + x_4

$`f_{4}: x_4 + x_5`$

$f_{5}: x_5 + x_{6} x_{7} x_{8} + 1$

$f_{6}: x_2 + x_4 + x_5 + x_{6} x_{7} x_{8} + x_{6} x_{8} + x_7 + x_8$

$f_{7}: x_1 + x_2 + x_4 + x_{6} x_{7} x_{8} + x_{6} x_{8} + x_{7} x_{8} + x_{7} + x_8 + 1$

$f_{8}: x_1 + x_2 + x_3 + x_4 + x_5 + x_{6} x_{7} x_{8} + x_{6} x_{7}+ x_{6} x_{8} + x_{7} x_{8}$

$f_{9}: x_2 + x_5 + x_6 x_7 + x_8$

$f_{10}: x_1 + x_6 x_7 x_8 + x_6 x_8 + x_6 + x_7 x_8 + x_7 + x_8$

$f_{11}: x_1 + x_4 + x_{6} x_{7} + x_{6} x_{8} + x_{6} + x_{7}$

$f_{12}: x_1 + x_3 + x_5 + x_{6} x_{7} x_{8} + x_{6} x_{8} + x_6 + x_{7} x_{8} + x_8 + 1$

$f_{13}: x_2 + x_3 + x_4 + x_{6} x_{7} x_{8} + x_{6} x_{8} + x_6 + x_{7} x_{8}$

$f_{14}: x_1 + x_3 + x_4 + x_5 + x_{6} x_{8} + x_{7} x_{8} + x_{7} + 1$

$f_{15}: x_1 + x_2 + x_3 + x_5 + x_6 + x_{7} x_{8} + x_8 + 1$

$f_{16}: x_3 + x_{6} x_{7} x_{8} + x_{6} x_{7} + x_6 + x_{7} x_{8} + x_7 + 1$

$x  =  (0,1,0,0,0,1,0,1)$

$F(x)=( 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 )$

$e   =  (0,0,0,0,0,0,0,0,0,1,0,1,0,0,0)$

$w(e)=2$

$y=F(x)+e=(1,1,0,0,1,1,0,0,0,1,1,0,0,0,0,0)$



First the change is applied to the original ANF function, to transform it into NNF. Then the problem Hamiltonian H_p = H_1 + ... + H_16 is constructed via the penalty.

For example,

```
f_10 (NNF)= -2 x_1 x_6 x_7 x_8 + 4 x_1 x_6 x_7 + 2 x_1 x_6 x_8 - 2 x_1 x_6 + 2 x_1 x_7 x_8 - 2 x_1 x_7 - 2 x_1 x_8 + x_1 + x_6 x_7 x_8 - 2 x_6 x_7 - x_6 x_8 + x_6 - x_7 x_8 + x_7 + x_8
```

```
H_10=2 x_1 x_6 x_7 x_8 - 4 x_1 x_6 x_7 - 2 x_1 x_6 x_8 + 2 x_1 x_6 - 2 x_1 x_7 x_8 + 2 x_1 x_7 + 2 x_1 x_8 - x1 - x_6 x_7 x_8 + 2 x_6 x_7 + x_6 x_8 - x_6 + x_7 x_8 - x_7 - x_8 + 1
```

When the ANF function is developed into NNF, the recursive action of Eq can introduce all possible combination of products of x_n in each function f_n. Every possible interaction of qubits x_{1},...,x_{8} is appeared in H_p since f_8 includes every terms.

In order to reduce H_p into 2-body interaction using a minimal amount of ancillary qubits, a special mapping inspired by graph is proposed.

Qubits are divided equally into two groups, within each they are fully connected. Each vertex represents a qubit and each edge represents an ancilla.

For n qubits, if H_p includes all possibilities for up to 3-qubit interaction. A total number of n(n+2)/4 qubits is needed to represent the problem Hamiltonian if n is even ((n+1)^2/4 if n is odd).

If $H_p$ has interaction up to $n$ qubits, it would need up to 2^((n+2)/2)-2 total qubits to represent all possible combination in each graph if n is even (3* 2^((n-1)/2)-2 if n is odd).

```
    $x_1x_2\longrightarrow x_{12}$\;
    $x_3x_4\longrightarrow x_{34}$\;
    $x_1x_3\longrightarrow x_{13}$\;
    $x_2x_4\longrightarrow x_{24}$\;
    $x_1x_4\longrightarrow x_{14}$\;
    $x_2x_3\longrightarrow x_{23}$\;
    \;
    $x_5x_6\longrightarrow x_{56}$\;
    $x_7x_8\longrightarrow x_{78}$\;
    $x_5x_7\longrightarrow x_{57}$\;
    $x_6x_8\longrightarrow x_{68}$\;
    $x_5x_8\longrightarrow x_{58}$\;
    $x_6x_7\longrightarrow x_{67}$\;
    \;
    $x_1x_{34}\longrightarrow x_{134}$\;
    $x_2x_{34}\longrightarrow x_{234}$\;
    $x_3x_{12}\longrightarrow x_{123}$\;
    $x_4x_{12}\longrightarrow x_{124}$\;
    \;
    $x_5x_{78}\longrightarrow x_{578}$\;
    $x_6x_{78}\longrightarrow x_{678}$\;
    $x_7x_{56}\longrightarrow x_{567}$\;
    $x_8x_{56}\longrightarrow x_{568}$\;
    \;
    $x_{12}x_{34}\longrightarrow x_{1234}$\;
    $x_{56}x_{78}\longrightarrow x_{5678}$\;
```

The maximum number of ancillary qubits is needed only if we have one term in the ANF function that includes all members x_1,...,x_n, which is not always the case.

Here is a 16-32 linear code example:

Defining polynomials (ANF): 

```

$x = (x_1,...,x_{16}) \longrightarrow F =(f_1,...,f_{32})$

$f_1 = x_1 + x_2 + x_3 + x_6 + x_7 + x_{10} + x_{12} + x_{16}$

$f_2 = x_1 + x_2 + x_3 + x_5 + x_7 + x_8 + x_{10} + x_{11} + x_{16}$

$f_3 = x_3 + x_4 + x_5 + x_7 + x_8 + x_{14}$

$f_4 = x_4 + x_5 + x_8 + x_{11} + x_{13} + x_{15} + x_{16}$

$f_5 = x_2 + x_3 + x_4 + x_7 + x_8 + x_{12} + x_{14} + x_{16}$

$f_6 = x_5 + x_6 + x_7 + x_8 + x_{11} + x_{12} + x_{14} + x_{16}$

$f_7 = x_2 + x_4 + x_6 + x_7 + x_9 + x_{10} + x_{11} + x_{12} + x_{15} + x_{16}$

$f_8 = x_5 + x_7 + x_8 + x_9 + x_{11} + x_{12} + x_{15}$

$f_9 = x_1 + x_2 + x_5 + x_6 + x_9 + x_{12} + x_{14}$

$f_{10} = x_1 + x_2 + x_4 + x_5 + x_7 + x_8 + x_{10} + x_{11} + x_{12} + x_{14} + x_{15} + x_{16}$

$f_{11} = x_3 + x_4 + x_6 + x_8 + x_{10} + x_{12} + x_{13} + x_{15} + x_{16}$

$f_{12} = x_1 + x_5 + x_8 + x_{11} + x_{13} + x_{16}$

$f_{13} = x_2 + x_3 + x_6 + x_9 + x_{12} + x_{14} + x_{15}$

$f_{14} = x_1 + x_2 + x_3 + x_4 + x_5 + x_7$

$f_{15} = x_2 + x_3 + x_6 + x_8 + x_9 + x_{12} + x_{14} + x_{16}$

$f_{16} = x_2 + x_4 + x_5 + x_6 + x_7 + x_{13} + x_{15} + x_{16}$

$f_{17} = x_4 + x_8 + x_9 + x_{11} + x_{12} + x_{13} + x_{15} + x_{16}$

$f_{18} = x_5 + x_6 + x_7 + x_{10} + x_{12} + x_{13} + x_{16}$

$f_{19} = x_2 + x_4 + x_5 + x_6 + x_7 + x_8 + x_{11} + x_{12} + x_{14} + x_{16}$

$f_{20} = x_2 + x_3 + x_7 + x_8 + x_{12} + x_{13} + x_{14} + x_{16}$

$f_{21} = x_3 + x_5 + x_7 + x_8 + x_{10} + x_{11}$

$f_{22} = x_2 + x_4 + x_5 + x_6 + x_8 + x_{11} + x_{12} + x_{14} + x_{15} + x_{16}$

$f_{23} = x_2 + x_3 + x_4 + x_9 + x_{10} + x_{11} + x_{13} + x_{15}$

$f_{24} = x_2 + x_3 + x_4 + x_6 + x_8 + x_{12} + x_{13} + x_{14}$

$f_{25} = x_1 + x_2 + x_3 + x_7 + x_8 + x_9 + x_{10} + x_{11} + x_{12} + x_{13} + x_{15} + x_{16}$

$f_{26} = x_1 + x_2 + x_5 + x_6 + x_7 + x_8 + x_9 + x_{10} + x_{11} + x_{14} + x_{15}$

$f_{27} = x_2 + x_4 + x_5 + x_6 + x_8 + x_9 + x_{10} + x_{11} + x_{16}$

$f_{28} = x_5 + x_6 + x_7 + x_{10} + x_{11} + x_{12} + x_{13} + x_{14}$

$f_{29} = x_1 + x_2 + x_3 + x_4 + x_5 + x_8 + x_{11} + x_{13} + x_{14} + x_{16}$

$f_{30} = x_2 + x_5 + x_6 + x_9 + x_{10} + x_{11} + x_{13} + x_{15}$

$f_{31} = x_1 + x_2 + x_3 + x_{13} + x_{14} + x_{15} + x_{16}$

$f_{32} = x_1 + x_3 + x_6 + x_8 + x_9 + x_{10} + x_{12} + x_{13} + x_{14} + x_{15} + x_{16}$

```
```

x  =  (1,1,0,1,0,0,0,1,0,0,1,1,1,1,1,1)

F(x) = (0,1,1,0,0,1,0,0,0,1,0,1,0,1,1,1,1,1,1,0,0,$
$0,1,0,0,0,1,0,0,0,0,1)

e   =  (0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,$
$0,0,0,0,0,0,0,0,0,0)

w(e)=2

y=F(x)+e=(0,1,1,0,0,1,0,0,0,1,1,1,0,1,1,0,1,$
$1,1,0,0,0,1,0,0,0,1,0,0,0,0,1)

```

In this example, the total number of x is 16, but the maximum number of terms included in a function is m=12 (for f_{10} and f_{25}).

In order to solve any 16-32 code, 2^((16+2)/2)-2=510 qubits is needed. But for this case, it can be reduced down to 401 qubits.

The total number of qubits scales exponentially with $m$, with an upper limit 2^((n+2)/2)-2.
