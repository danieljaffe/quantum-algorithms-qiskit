from qiskit import *
from qiskit.quantum_info.operators import Operator
from math import floor, pi, sqrt
import matplotlib.pyplot as plt
import numpy as np


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


def initialize(n):
    """
    This function sets initial states and applies Hadamards to each qubit
    Note: apply_H isn't called because it is actually more efficient to initialize in one loop as opposed to 2.
    """
    # apply H to first n qubits and X H to the last qubit (ancilla qubit)
    quantum_register = QuantumRegister(n)
    classical_register = ClassicalRegister(n)
    quantum_circuit = QuantumCircuit(quantum_register, classical_register)
    for index in range(n):
        quantum_circuit.h(quantum_register[index])
    quantum_circuit.barrier()
    return quantum_circuit, quantum_register, classical_register


def get_Z0(n):
    """
    This function generates the Z0 gate satisfying the conditions for x in {0,1}^n Z0|x> -> -|x> iff x = 0^n
        otherwise Z0|x> -> |x>
    The parameter to this function is only the size n, a 2^n x 2^n dimensional matrix is created satisfying the
        conditions above.
    This function has one dependency, the DefGate function defined in pyquil.quil
    This function is designed to absorb the negative in G, so the returned gate is actually -Z0
    Returns -Z0
    """
    # Create a 2^n x 2^n matrix with all 0's
    gate = np.zeros((2 ** n, 2 ** n), dtype=int)
    # since this is -Z0, set first element to 1 not -1
    gate[0][0] = 1
    # set all other elements on the diagonal to -1, again not 1 because this is -Z0
    for i in range(1, 2 ** n):
        gate[i][i] = -1
    # Return gate
    return Operator(gate)


def get_Zf(f, n):
    """
    This function generates the Zf gate satisfying the condition for x in {0,1}^n where Zf|x> -> (-1)^f(X)|x>
    This function requires that f(x) be calculated for all x, so f is passed as an anonymous function, the other
        parameter is n.
    The function has one dependency, the DefGate function defined in pyquil.quil
    This function finds all permutations of bitstrings of length n, then initializes a 2^n x 2^n matrix of all 0's,
        and sets all elements along the diagonal to either 1 or -1 depending on f(x)
    Finally a gate representation of this matrix is returned.
    """
    # generate bitstring permutations
    bitstrings = list()
    get_bitstring_permutations(0, bitstrings, n, [0] * n)
    # initialize a 2^n x 2^n matrix of all 0's
    gate = np.zeros((2 ** n, 2 ** n), dtype=int)
    # set diagonals of matrix based on f(x)
    for i in range(2 ** n):
        gate[i][i] = -1 if f(bitstrings[i]) == 1 else 1
    # create and return gate
    return Operator(gate)


def grovers_algorithm(f, n, shots=1024):
    """
    This function is intended to determine if there exists an x in {0,1}^n s.t. f(x) = 1 for a given function f s.t.
        f:{0,1}^n -> {0,1}. The algorithm first constructs Zf, -Z0 gates, initializes with Hanamard matrices, and
        applies G = -H^n o Z0 o H^n o Zf. This algorithm is not deterministic, so G is applied multiple times. More
        specifically, G is run (pi / 4 * sqrt(n)) times. Furthermore, there are 10 trials to minimize the chance of a
        false negative.
    This function has an anonymous function and integer n as parameters.
    This function runs the algorithm as described for each 10 trials, and then checks if for any of the outputted states
        x, if f(x) = 1. If this is true, then 1 is returned, otherwise 0 is returned. The function returns 0 if there
        is an issue with the simulator.
    This function uses 9q-squared-qvm, so it assumes that n <= 9
    """
    # Initialize the circuit and apply Hadamards to all qubits
    quantum_circuit, quantum_register, classical_register = initialize(n)

    # Needed for application of custom gates since Operator works in reverse
    rv_qr = list()
    for index in range(n):
        rv_qr.append(index)
    rv_qr.reverse()

    # Define and generate Z0 gate (really -Z0)
    z0_gate = get_Z0(n)

    # Define and generate Zf gate
    zf_gate = get_Zf(f, n)

    # Determine the number of times to apply G
    iteration_count = floor(pi / 4 * sqrt(2 ** n))
    # Apply G iteration_count times
    for i in range(iteration_count):
        # Apply Zf
        quantum_circuit.unitary(zf_gate, rv_qr)
        # Apply H to all qubits
        for index in range(n):
            quantum_circuit.h(quantum_register[index])
        # Apply -Z0
        quantum_circuit.unitary(z0_gate, rv_qr)
        # Apply H to all qubits
        for index in range(n):
            quantum_circuit.h(quantum_register[index])
        quantum_circuit.barrier()

    # Run simulator
    quantum_simulator = Aer.get_backend('qasm_simulator')
    quantum_circuit.measure(quantum_register[0:n], classical_register)

    # Display circuit diagram
    quantum_circuit.draw('mpl')
    plt.show()

    # Execute and evaluate the job results
    job = execute(quantum_circuit, quantum_simulator, shots=shots)
    results = job.result()
    counts = results.get_counts(quantum_circuit)

    # Parse results and return 1 or 0 accordingly
    dict = {}
    for key in counts:
        if counts[key] >= (shots/(2**n)):
           dict[key] = counts[key]
    for key in dict:
        element = list(key)
        element = [int(i) for i in element]
        element.reverse()
        if f(element) == 1:
            return 1
    return 0


def f(args):
    return 0
    # n = len(args)
    # marked = [0]*(n-1)
    # marked.append(1)
    # if(args == marked):
    #     return 1
    # return 0
    # return x[0]


grovers_algorithm(f, 4, 1024)
