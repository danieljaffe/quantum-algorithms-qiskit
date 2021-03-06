from qiskit import *
from qiskit.visualization import plot_histogram
from qiskit.quantum_info.operators import Operator
from qiskit.tools.monitor import job_monitor
import numpy as np
import matplotlib.pyplot as plt


def apply_H(circuit, apply_to_list):
    """
    Apply Hadamards to all specified qubits (if apply_to_list[index] == 1).
    Designed for a large amount of Hadamards being applied at once.
    """
    for index, qubit in enumerate(apply_to_list):
        if qubit == 1:
            circuit.h(index)
    return circuit


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
    Parameters: f is an anonymous function and n is the number of bits in input: f:{0,1}^n -> {0,1}
    This function returns an oracle gate representing the function f
        for all x in {0,1}^n and y in {0,1}, the desired result of the oracle is mapping the input
        |x,y> to |x, y + f(x)> where + is addition modulo 2. The function first finds the list of bitstring
        permutations of n bits, it then establishes a mapping which is representative of the decimal number of the
        bitstring represents. It then determines for each |x,y>, it calculates |x, f(x) + y>. Finally it constructs a
        permutation gate which treats each permutation as a different basis vector in the 2^(n+1) dimensional complex
        hilbert space that represents a system of n + 1 qubits. The permutation gate is returned.
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


def initialize(states):
    """
    This function sets initial states and applies Hadamards to each qubit
    Note: apply_H isn't called because it is actually more efficient to initialize in one loop as opposed to 2.
    """
    n = len(states)
    circuit = QuantumCircuit(n, n - 1)
    for index, state in enumerate(states):
        if state == 1:
            circuit.x(index)
        circuit.h(index)
    return circuit


def bv_algorithm(f, n, shots=1024, token=""):
    """
    This function is intended to determine to find a,b s.t. f(x) = a*x + b where a*x is bitwise inner product and +
        is addition modulo 2, for a given function f s.t. f:{0,1}^n -> {0,1}. It is assumed that f can be represented as
        f(x) = a*x + b. The algorithm first calculates b by solving for f(0^n). It then initializes the qubits with H
        for the first n qubits and X and H for the last qubit. The algorithm then constructs a Uf oracle gate based on
        the function input f. It then applies Uf to all the qubits and applies H to the first n qubits. Finally, the
        simulator is run on the circuit and measures the results. The function then returns the first n measured qubits.
    
    UPDATE: This program will now handle running on IBM's quantum computers. A new parameter token is provided which
            can be used to pass in a user token. If no token is specified, the running machine will attempt to use
            a previously saved token if one can be found. If no token is found, the program will default to running on
            the simulator.
    
    This function has an anonymous function and integer n as parameters.
    """
    # Account and backend setup
    using_simulator = False
    if token != "":
        # Sets the IBMQ token
        IBMQ.save_account(token)
    if comp != "":
        # Attempts to load IBMQ based on a previously stored token
        IBMQ.load_account()
        provider = IBMQ.get_provider('ibm-q')
        backend = provider.get_backend("ibmq_16_melbourne")
    else:
        # Failure loading an IBMQ account will default to simulator usage
        print("Error in loading IBMQ account. Running simulation instead.")
        backend = Aer.get_backend('qasm_simulator')
        using_simulator = True

    initialize_list = [0] * n

    # calculate b by f(0^n) = b
    b = f(initialize_list)
    print("b is: ", b)

    # Initialize circuit by applying H to first n qubits and X H to last qubit (ancilla qubit)
    initialize_list.append(1)
    qubits = list(range(len(initialize_list)))
    circuit = initialize(initialize_list)

    # Generate Uf oracle from f (anonymous function)
    uf_gate = generate_uf(f, n)
    # Applying the uf_gate must be done with the qubits in the reverse order due to the implementation of qiskit
    rv_qubits = qubits[::-1]
    circuit.unitary(uf_gate, rv_qubits)

    # Apply H to all qubits except for the last qubit
    apply_to_list = [1] * n
    apply_to_list.append(0)
    circuit = apply_H(circuit, apply_to_list)
    circuit.measure(range(n), range(n))

    # Draw the circle
    circuit.draw('mpl')
    plt.show()

    # run circuit and measure qubits
    job = execute(circuit, backend, shots=shots, optimization_level=3)
    if not using_simulator:
        job_monitor(job)

    try:
        result = job.result()
    except:
        print(job.error_message())
        return
    print("Time taken:", result.time_taken)
    counts = result.get_counts(circuit)

    plot_histogram(counts)

    # Will return the key with the highest number of hits. While the simulation will be deterministic,
    # the any real quantum computer will have errors. This helps to mitigate those errors
    max = 0
    for count in counts:
        if max < counts[count]:
            max = counts[count]
            a = count
    return a, b