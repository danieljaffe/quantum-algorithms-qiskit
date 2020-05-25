from qiskit import *
from qiskit.quantum_info.operators import Operator
import numpy as np
import matplotlib.pyplot as plt


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


def generate_uf(f, n):
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
    get_bitstring_permutations(0, bitstrings, n + 1, [0] * (n + 1))

    # initialize mapping and permutation list
    perm_dict = dict()
    perm_list = []

    # populate mapping
    for permutation, bitstring in enumerate(bitstrings):
        values = [0] * (len(bitstrings))
        values[permutation] = 1
        perm_dict["".join(str(bit) for bit in bitstring)] = values

    # Send each |xy> to |x, f(x) + y>
    for bitstring in bitstrings:
        params = bitstring[:n]
        params.append((f(params) + bitstring[-1]) % 2)
        perm_list.append(perm_dict["".join(str(bit) for bit in params)])
    return Operator(np.array(perm_list))


def initialize(n):
    """
    This function sets initial states and applies Hadamards to each qubit
    Note: apply_H isn't called because it is actually more efficient to initialize in one loop as opposed to 2.
    """
    # apply H to first n qubits and X H to the last qubit (ancilla qubit)
    quantum_register = QuantumRegister(n + 1)
    classical_register = ClassicalRegister(n)
    quantum_circuit = QuantumCircuit(quantum_register, classical_register)
    # In qiskit, all quantum registers start in the low energy |0> state so we must apply an x gate to our helper bit
    # in order to put it in the state |1>
    quantum_circuit.x(quantum_register[-1])
    for index in range(n + 1):
        quantum_circuit.h(quantum_register[index])
    quantum_circuit.barrier()
    return quantum_circuit, quantum_register, classical_register


def deutsch_jozsa_algorithm(f, n, shots=1024, threshold=0.9):
    """
    This function is intended to determine if f is constant or balanced for a given function f s.t.
        f:{0,1}^n -> {0,1}. The algorithm initializes the qubits with H for the first n qubits and X and H for the last
        qubit. The algorithm then constructs a Uf oracle gate based on the function input f. It then applies Uf to all
        the qubits and applies H to the first n qubits. Finally, the simulator is run on the circuit and measures the
        results. If upon measurement, the first n qubits are all 0, 1 is returned and the function is constant,
        otherwise 0 is returned and the function is balanced.
    This function has an anonymous function and integer n as parameters.
    This function uses 9q-squared-qvm, so it assumes that n <= 9.
    """
    # apply H to first n qubits and X H to the last qubit (ancilla qubit)
    quantum_circuit, quantum_register, classical_register = initialize(n)

    # Generate Uf oracle from f (anonymous function)
    uf_gate = generate_uf(f, n)
    # Applying the uf_gate must be done with the qubits in the reverse order due to the implementation of qiskit
    rv_qr = list()
    for index in range(n + 1):
        rv_qr.append(index)
    rv_qr.reverse()
    quantum_circuit.unitary(uf_gate, rv_qr)
    quantum_circuit.barrier()

    # Apply Hadamards to first n qubits
    for index in range(n):
        quantum_circuit.h(quantum_register[index])

    # Run simulator
    quantum_simulator = Aer.get_backend('qasm_simulator')
    quantum_circuit.measure(quantum_register[0:n], classical_register)
    quantum_circuit.draw('mpl')
    plt.show()

    # Evaluate the job results
    job = execute(quantum_circuit, quantum_simulator, shots=shots)
    results = job.result()
    counts = results.get_counts(quantum_circuit)

    # NOTE: 1 = constant, 0 = balanced
    # To compensate for the error rates of actual quantum machines, the threshold for being consider balanced is
    # set to be threshold * shots, or a percentage of shots given at runtime.

    # Function is constant
    key = '0' * n
    if key in counts:
        if counts[key] >= threshold * shots:
            return 1
    key = '1' * n
    if key in counts:
        if counts[key] >= threshold * shots:
            return 1

    # Function is balanced
    return 0


def f(args):
    # return 1
    return (args[0] + args[1]) % 2
    # return args[0]

print(deutsch_jozsa_algorithm(f, 2))