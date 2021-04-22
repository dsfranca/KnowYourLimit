
import cirq
import recirq
import networkx as nx
from cirq.contrib.svg import SVGCircuit
import numpy as np
import matplotlib.pyplot as plt

from pytket.predicates import CompilationUnit, ConnectivityPredicate
from pytket.passes import SequencePass, RoutingPass, DecomposeSwapsToCXs
from pytket.routing import GraphPlacement



def rzz(rads):
    """Returns a gate with the matrix exp(-i Z⊗Z rads)."""
    return cirq.ZZPowGate(exponent=2 * rads / np.pi, global_shift=-0.5)


def qaoa_max_cut_unitary(qubits, betas, gammas, graph):  # Nodes should be integers
    for beta, gamma in zip(betas, gammas):
        yield (rzz(-0.5 * gamma*problem_graph[i][j]['weight']).on(qubits[i], qubits[j]) for i, j in graph.edges)
        yield cirq.rx(2 * beta).on_each(*qubits)


def qaoa_max_cut_circuit(qubits, betas, gammas, graph):  # Nodes should be integers
    return cirq.Circuit(
        # Prepare uniform superposition
        cirq.H.on_each(*qubits),
        # Apply QAOA unitary
        qaoa_max_cut_unitary(qubits, betas, gammas, graph),
        # Measure
        cirq.measure(*qubits, key='m'),
    )


def cut_values(bitstrings, graph):
    mat = networkx.adjacency_matrix(graph, nodelist=sorted(graph.nodes))
    vecs = (-1) ** bitstrings
    vals = 0.5 * np.sum(vecs * (mat @ vecs.T).T, axis=-1)
    vals = 0.5 * (graph.size() - vals)
    return vals
## define the length
#length = 3
#
## define the qubits
#qubits = [cirq.GridQubit(i, j) for i in range(length) for j in range(length)]
#
## print the qubits
#print(qubits)
#
#
## define a circuit
#circuit = cirq.Circuit()
#circuit.append(cirq.H(q) for q in qubits if (q.row + q.col) % 2 == 0)
#print(circuit)
#
#circuit.append(cirq.X(q) for q in qubits if (q.row + q.col) % 2 == 1)
#
#circuit.append(cirq.X(q) for q in qubits if (q.row + q.col) % 2 == 0)
#
#
#
## Set problem parameters
#n = 9
#p = 2
#
## Generate a random 3-regular graph on n nodes
#graph = nx.random_regular_graph(2, n)
#
## Make qubits
#qubits = cirq.LineQubit.range(n)
#qubits2 = [cirq.GridQubit(i, j) for i in range(length) for j in range(length)]
## Print an example circuit
#betas = np.random.uniform(-np.pi, np.pi, size=p)
#gammas = np.random.uniform(-np.pi, np.pi, size=p)
#circuit = qaoa_max_cut_circuit(qubits, betas, gammas, graph)
#print('Example QAOA circuit:')
#print(circuit)
#depth=0
#for moments in circuit:
#    print("new layer\n")
#    print(moments)
#    print("\n")
#    depth+=1
#    
#depth2=0 
#circuit2 = qaoa_max_cut_circuit(qubits2, betas, gammas, graph)
#print('Example QAOA circuit:')
#print(circuit2)
#
#for moments in circuit2:
#    print("new layer\n")
#    print(moments)
#    print("\n")
#    depth2+=1
    


from recirq.qaoa.problem_circuits import get_generic_qaoa_circuit
from recirq.qaoa.gates_and_compilation import compile_problem_unitary_to_arbitrary_zz, \
    compile_driver_unitary_to_rx
all_to_all=[]
square=[]
for k in range(22,23):
#k=50
#problem_graph = nx.random_regular_graph(d=4, n=k**2)
problem_graph=nx.complete_graph(k**2)
nx.set_edge_attributes(problem_graph, values=1, name='weight')
circuit_qubits = cirq.LineQubit.range(1, k**2+1)
gammas = [0.1,0.3,0.5]
betas = [0.5,0.3,0.1]
circuit = get_generic_qaoa_circuit(
    problem_graph=problem_graph,
    qubits=circuit_qubits,
    gammas=gammas,
    betas=betas)
circuit = compile_problem_unitary_to_arbitrary_zz(circuit)
circuit = compile_driver_unitary_to_rx(circuit)
SVGCircuit(circuit)
initial_depth=0
for moment in circuit:
    initial_depth+=1
all_to_all.append(initial_depth)
print(initial_depth)
import cirq.contrib.routing as ccr

uncompiled_c_graph = ccr.get_circuit_connectivity(circuit)
#nx.draw_networkx(uncompiled_c_graph)
plt.show()
plt.close()
import cirq.google as cg

dev_graph = ccr.xmon_device_to_graph(cg.Sycamore23)
#nx.draw_networkx(dev_graph)
plt.show()
device = cg.Sycamore23

#device=nx.grid_2d_graph(4, 5)
import pytket
from pytket.circuit import Node

def _qubit_index_edges():
    dev_graph = ccr.xmon_device_to_graph(device)
    for n1, n2 in dev_graph.edges:
        yield Node('grid', n1.row, n1.col), Node('grid', n2.row, n2.col)

def _device_to_tket_device():
    arc = pytket.routing.Architecture(
        list(_qubit_index_edges())
    )
    return pytket.device.Device({}, {}, arc)

tk_circuit = pytket.extensions.cirq.cirq_to_tk(circuit)
# pytket.device.Device({}, {}, pytket.routing.SquareGrid(7,7))
tk_device =_device_to_tket_device()
tk_device=pytket.device.Device({}, {}, pytket.routing.SquareGrid(k,k))


from pytket.predicates import CompilationUnit, ConnectivityPredicate
from pytket.passes import SequencePass, RoutingPass, DecomposeSwapsToCXs, PlacementPass
from pytket.routing import GraphPlacement



unit = CompilationUnit(tk_circuit, [ConnectivityPredicate(tk_device)])
passes = SequencePass([
    PlacementPass(GraphPlacement(tk_device)),
    RoutingPass(tk_device)])
passes.apply(unit)
valid = unit.check_all_predicates()
assert valid
unit.initial_map
def tk_to_cirq_qubit(tk):
    ind = tk.index
    print(ind)
    #print(ind[0])
    #print(*ind,len(ind))
    return cirq.LineQubit(ind[0]) if len(ind) == 1 else cirq.GridQubit(*ind)


#initial_map2 = {(n1): (n2) for n1, n2 in unit.initial_map.items()}
#initial_map = {tk_to_cirq_qubit(n1): tk_to_cirq_qubit(n2) for n1, n2 in unit.initial_map.items()}
#initial_map
#unit.final_map
#final_map = {tk_to_cirq_qubit(n1): tk_to_cirq_qubit(n2)
#             for n1, n2 in unit.final_map.items()}
#final_map
unit.circuit.qubits
routed_circuit = pytket.extensions.cirq.tk_to_cirq(unit.circuit)
SVGCircuit(routed_circuit)
routed_c_graph = ccr.get_circuit_connectivity(routed_circuit)
plt.close()
#nx.draw_networkx(routed_c_graph)
initial_depth=0
for moment in routed_circuit:
    initial_depth+=1
print(initial_depth)
square.append(initial_depth)
    
    
    
moments=[] 
for moment in circuit:
    moments.append(moment)
for k in moments[1]:
    C=k