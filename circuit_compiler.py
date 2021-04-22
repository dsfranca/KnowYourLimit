

import cirq
import recirq
import networkx as nx
from cirq.contrib.svg import SVGCircuit
import numpy as np
import matplotlib.pyplot as plt

from pytket.predicates import CompilationUnit, ConnectivityPredicate
from pytket.passes import SequencePass, RoutingPass, DecomposeSwapsToCXs
from pytket.routing import GraphPlacement
from recirq.qaoa.problem_circuits import get_generic_qaoa_circuit
from recirq.qaoa.gates_and_compilation import compile_problem_unitary_to_arbitrary_zz, \
    compile_driver_unitary_to_rx
import pytket
from pytket.circuit import Node

from pytket.predicates import CompilationUnit, ConnectivityPredicate
from pytket.passes import SequencePass, RoutingPass, DecomposeSwapsToCXs, PlacementPass
from pytket.routing import GraphPlacement
import cirq.google as cg
import cirq.contrib.routing as ccr


#some of these functions were copy-pasted from QAOA from Google. How to deal with this?
def rzz(rads):
    """Returns a gate with the matrix exp(-i Z⊗Z rads)."""
    return cirq.ZZPowGate(exponent=2 * rads / np.pi, global_shift=-0.5)


def qaoa_weighted_max_cut_unitary(qubits, betas, gammas, graph):  # Nodes should be integers
    for beta, gamma in zip(betas, gammas):
        yield (rzz(-0.5 * gamma*problem_graph[i][j]['weight']).on(qubits[i], qubits[j]) for i, j in graph.edges)
        yield cirq.rx(2 * beta).on_each(*qubits)

def qaoa_max_cut_unitary(qubits, betas, gammas, graph):  # Nodes should be integers
    for beta, gamma in zip(betas, gammas):
        yield (rzz(-0.5 * gamma).on(qubits[i], qubits[j]) for i, j in graph.edges)
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
    




def tk_to_cirq_qubit(tk):
    ind = tk.index
    print(ind)
    #print(ind[0])
    #print(*ind,len(ind))
    return cirq.LineQubit(ind[0]) if len(ind) == 1 else cirq.GridQubit(*ind)
    
def compiled_routed_qaoa(problem_graph,gammas,betas,tk_device):
    k=problem_graph.order()
    #have to adjust later
    nx.set_edge_attributes(problem_graph, values=1, name='weight')
    circuit_qubits = cirq.LineQubit.range(1, k+1)

    circuit = get_generic_qaoa_circuit(
        problem_graph=problem_graph,
        qubits=circuit_qubits,
        gammas=gammas,
        betas=betas)
    circuit = compile_problem_unitary_to_arbitrary_zz(circuit)
    circuit = compile_driver_unitary_to_rx(circuit)
    SVGCircuit(circuit)

    
#    uncompiled_c_graph = ccr.get_circuit_connectivity(circuit)
#    #nx.draw_networkx(uncompiled_c_graph)
#    plt.show()
#    plt.close()
    
    
#    dev_graph = ccr.xmon_device_to_graph(cg.Sycamore23)
#    #nx.draw_networkx(dev_graph)
#    plt.show()
#    device = cg.Sycamore23
    
    #device=nx.grid_2d_graph(4, 5)

    

    
    tk_circuit = pytket.extensions.cirq.cirq_to_tk(circuit)
    # pytket.device.Device({}, {}, pytket.routing.SquareGrid(7,7))
    
    
    

    
    
    
    unit = CompilationUnit(tk_circuit, [ConnectivityPredicate(tk_device)])
    passes = SequencePass([
        PlacementPass(GraphPlacement(tk_device)),
        RoutingPass(tk_device)])
    passes.apply(unit)
    valid = unit.check_all_predicates()
    assert valid
    unit.initial_map

    
    
    #initial_map2 = {(n1): (n2) for n1, n2 in unit.initial_map.items()}
    #initial_map = {tk_to_cirq_qubit(n1): tk_to_cirq_qubit(n2) for n1, n2 in unit.initial_map.items()}
    #initial_map
    #unit.final_map
    #final_map = {tk_to_cirq_qubit(n1): tk_to_cirq_qubit(n2)
    #             for n1, n2 in unit.final_map.items()}
    #final_map
    unit.circuit.qubits
    routed_circuit = pytket.extensions.cirq.tk_to_cirq(unit.circuit)

    #nx.draw_networkx(routed_c_graph)

    return routed_circuit


def compiled_routed_weighted_qaoa(problem_graph,gammas,betas,tk_device):
    k=problem_graph.order()
    circuit_qubits = cirq.LineQubit.range(1, k+1)
    try:
        circuit = get_generic_qaoa_circuit(
            problem_graph=problem_graph,
            qubits=circuit_qubits,
            gammas=gammas,
            betas=betas)
    except ValueError:
        print("Setting edge weights to 1")
        nx.set_edge_attributes(problem_graph, values=1, name='weight')
        circuit = get_generic_qaoa_circuit(
            problem_graph=problem_graph,
            qubits=circuit_qubits,
            gammas=gammas,
            betas=betas)
        
        
    circuit = compile_problem_unitary_to_arbitrary_zz(circuit)
    circuit = compile_driver_unitary_to_rx(circuit)
    SVGCircuit(circuit)

    
    
#    uncompiled_c_graph = ccr.get_circuit_connectivity(circuit)
#    #nx.draw_networkx(uncompiled_c_graph)
#    plt.show()
#    plt.close()
    
    
#    dev_graph = ccr.xmon_device_to_graph(cg.Sycamore23)
#    #nx.draw_networkx(dev_graph)
#    plt.show()
#    device = cg.Sycamore23
    
    #device=nx.grid_2d_graph(4, 5)

    

    
    tk_circuit = pytket.extensions.cirq.cirq_to_tk(circuit)
    # pytket.device.Device({}, {}, pytket.routing.SquareGrid(7,7))
    
    
    

    
    
    
    unit = CompilationUnit(tk_circuit, [ConnectivityPredicate(tk_device)])
    passes = SequencePass([
        PlacementPass(GraphPlacement(tk_device)),
        RoutingPass(tk_device)])
    passes.apply(unit)
    valid = unit.check_all_predicates()
    assert valid
    unit.initial_map

    
    
    #initial_map2 = {(n1): (n2) for n1, n2 in unit.initial_map.items()}
    #initial_map = {tk_to_cirq_qubit(n1): tk_to_cirq_qubit(n2) for n1, n2 in unit.initial_map.items()}
    #initial_map
    #unit.final_map
    #final_map = {tk_to_cirq_qubit(n1): tk_to_cirq_qubit(n2)
    #             for n1, n2 in unit.final_map.items()}
    #final_map
    unit.circuit.qubits
    routed_circuit = pytket.extensions.cirq.tk_to_cirq(unit.circuit)

    #nx.draw_networkx(routed_c_graph)

    return routed_circuit





