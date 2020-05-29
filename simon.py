import numpy as np
from qiskit import *
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt
from qiskit import Aer
from qiskit.quantum_info.operators import Operator

def get_bitstring_permutations(index, lst, n, args):
    """
    This function populates a list with all the combinations of bit strings of length n
    """
    if (index == n):
        # append combination to list
        lst.append(list.copy(args))
    else:
        # handle the case where the preceding bit is 0
        args[index] = 0
        get_bitstring_permutations(index + 1, lst, n, args)

        # handle the case where the preceding bit is 1
        args[index] = 1
        get_bitstring_permutations(index + 1, lst, n, args)
   
   
def generate_uf_simons(f, n):
    """
    Parameters: f is an anonymous function and n is the number of bits in input: f:{0,1}^n -> {0,1}^n
    This function returns an oracle gate representing the function f for all x in {0,1}^n and y in {0,1}^n,
    the desired result is of the oracle is mapping the input x,y> to |x, y + f(x)> where + is addition modulo 2.
    The function first finds the list of bitstring permutations of n bits, it then establishes a mapping which is
    representative of the decimal number of the bitstring represents. For each |x,y>, it calculates
    |x, y + f(x)>. Finally it constructs a permutation gate which treats each permutation as a different basis vector
    in the 2^(n+1) dimensional complex Hilbert space that represents a system of 2*n qubits.
    Returns: the permutation gate
    """

    # generate list of all bitstrings of size n
    bitstrings = []
    get_bitstring_permutations(0, bitstrings, 2 * n, [0] * (2 * n))
    
    # initialize mapping and permutation list
    perm_dict = dict()
    perm_list = []
    
    # populate mapping
    for permutation, bitstring in enumerate(bitstrings):
        values = [0]*(len(bitstrings))
        values[permutation] = 1
        perm_dict["".join(str(bit) for bit in bitstring)] = values

    # Send each |xy> to |x, f(x) + y>
    for bitstring in bitstrings:
        params = bitstring[:n]
        params2 = bitstring[n:2 * n]
        f_values = f(bitstring[:n])
        for i in range(n):
            params.append((params2[i] + f_values[i]) % 2)
        perm_list.append(perm_dict["".join(str(bit) for bit in params)])
    return Operator(np.array(perm_list))
    

def simons_solver(Y, n):
    """
    Inputs: Y is a linear system of n-1 equations in matrix form. n is the dimension of the input into f.
    This function acts as a binary linear matrix solver.
    Returns: the key string s, if found, or the zero bitstring
    """
    # Create all possible bit strings to test for s
    bitstrings = []
    get_bitstring_permutations(0, bitstrings, n, [0] * n)

    # For each possible s, test to see if it's a candidate
    for s in bitstrings:
        if s == [0] * n:
            continue
        candidate = True

        # For each equation in Y, bit by bit test that y*s = [0]*n
        for y in Y:
            value = 0
            for i in np.arange(n):
                value = value + s[i] * y[i]

            # If a bit doesn't evaluate to 0...
            if (value % 2 == 1):
                candidate = False
        if (candidate):
            return s

    return [0] * n


def simons_algorithm(f, n):
    oracle = generate_uf_simons(f, n)
    circuit = QuantumCircuit(2*n, 2*n)
    simulator = Aer.get_backend("qasm_simulator")
    indices = list(range(2*n))
    indices.reverse()
    for i in range(n):
        circuit.h(i)
    circuit.unitary(oracle, indices, label="oracle")
    for i in range(n):
        circuit.h(i)
    indices = list(range(n))
    circuit.measure(indices, indices)
#    circuit.draw('mpl')
#    plt.show()
    for i in range(20):
        s = set()
        s_trials = []
        while(len(s) < n-1):
            job = execute(circuit, simulator, shots=1)
            result = job.result()
            counts = result.get_counts()
            for count in counts:
                s.add(count[2*n:n-1:-1])
        for bitstring in s:
            s_trials.append([int(bit) for bit in bitstring])
        s_trials = np.array(s_trials)
        val = simons_solver(s_trials, n)
        if val == [0] * n:
            continue
        f_val = f(val)
        if f_val == f([0]*n):
            return val
    return [0] * n
  #  print(counts)
        
    
def f(args):
    return [0]

def f3(x):
    if x == [0, 0, 0]:
        return [1, 0, 1]
    elif x == [0, 0, 1]:
        return [0, 1, 0]
    elif x == [0, 1, 0]:
        return [0, 0, 0]
    elif x == [0, 1, 1]:
        return [1, 1, 0]
    elif x == [1, 0, 0]:
        return [0, 0, 0]
    elif x == [1, 0, 1]:
        return [1, 1, 0]
    elif x == [1, 1, 0]:
        return [1, 0, 1]
    elif x == [1, 1, 1]:
        return [0, 1, 0]
    else:
        return NULL
        
def f2(x):
    if x == [0, 0]:
        return [0, 1]
    elif x == [0, 1]:
        return [1, 1]
    elif x == [1, 0]:
        return [0, 1]
    elif x == [1, 1]:
        return [1, 1]
        
print(simons_algorithm(f3,3))
