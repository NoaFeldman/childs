import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

outPhase = 1
n = 2
circuitNodes = [1 for i in range(2**n)]
circuitEdges = []
daggerNodes = [1 for i in range(2**n)]
daggerEdges = []
fullPhase = 8


def addNode(wireNum, option='circuit'):
    if option == 'circuit':
        nodes = circuitNodes
        edges = circuitEdges
    else:
        nodes = daggerNodes
        edges = daggerEdges
    nodes[wireNum] += 1
    edges.append([wireNum, nodes[wireNum] - 1, wireNum, nodes[wireNum] - 2])

def breakGates():
    for wire in range(2**n):
        addNode(wire)
        addNode(wire, option='dagger')
    global outPhase
    outPhase += 2


def addPhase(qubit: int):
    for wire in range(2**n):
        if wire & 2**(n - qubit - 1) != 0:
            addNode(wire)
            for i in range(fullPhase - 1):
                addNode(wire, option='dagger')

def addUcBasicGadget(qubit: int, option='circuit'):
    if option == 'circuit':
        nodes = circuitNodes
        edges = circuitEdges
    else:
        nodes = daggerNodes
        edges = daggerEdges
    for wire in range(2**n):
        if wire & 2**qubit > 0:
            wire2 = wire ^ 2**qubit
            currNode = nodes[wire]
            addNode(wire, option)
            nodes[wire] += 3
            edges.append([wire, currNode, wire, currNode + 1])
            edges.append([wire, currNode, wire, currNode + 3])
            edges.append([wire, currNode + 2, wire, currNode + 3])
            currNode2 = nodes[wire2]
            addNode(wire2, option)
            nodes[wire2] += 1
            edges.append([wire, currNode + 1, wire2, currNode2])
            edges.append([wire2, currNode2, wire2, currNode2 + 1])
            edges.append([wire2, currNode2 + 1, wire, currNode + 2])

def addUc(qubit: int):
    # circuit Uc
    addUcBasicGadget(qubit)
    # Uc^dagger = -ZUcZ
    for wire in range(2**n):
        if wire & 2**qubit > 0:
            for i in range(int(fullPhase / 2)):
                addNode(wire, option='dagger')
    addUcBasicGadget(qubit, option='dagger')
    for wire in range(2 ** n):
        if wire & 2 ** qubit == 0:
            for i in range(int(fullPhase / 2)):
                addNode(wire, option='dagger')
    global outPhase
    outPhase += 1

def addCnot(control: int, target: int):
    for wire in range(2**n):
        if wire & 2**control == 0:
            addNode(wire, 'circuit')
            addNode(wire, 'dagger')
        else:
            if wire < wire ^ 2**target:
                circuitNodes[wire] += 1
                circuitNodes[wire ^ 2**target] += 1
                circuitEdges.append([wire, circuitNodes[wire] - 1, wire ^ 2**target, circuitNodes[wire ^ 2**target] - 2])
                circuitEdges.append([wire, circuitNodes[wire] - 2, wire ^ 2**target, circuitNodes[wire ^ 2**target] - 1])
                daggerNodes[wire] += 1
                daggerNodes[wire ^ 2**target] += 1
                daggerEdges.append([wire, daggerNodes[wire] - 1, wire ^ 2**target, daggerNodes[wire ^ 2**target] - 2])
                daggerEdges.append([wire, daggerNodes[wire] - 2, wire ^ 2**target, daggerNodes[wire ^ 2**target] - 1])
    global outPhase
    outPhase += 2

def findEdge():
    for i in range(len(H)):
        if np.sum(H[i, :]) == 1:
            return i
    return -1

def removeTrails():
    while True:
        edge = findEdge()
        if edge == -1:
            return
        else:
            global H
            H = np.delete(np.delete(H, edge, 0), edge, 1)

addUc(0)
addCnot(0, 1)
midWiresLength = fullPhase - outPhase % fullPhase
circleWireLength = 6
circuitSize = sum(circuitNodes)
daggerSize = sum(daggerNodes)
nodesNum = circuitSize + daggerSize + midWiresLength * 2**n + circleWireLength
H = np.zeros((nodesNum, nodesNum))
# first nodes are the circuit nodes
for edge in circuitEdges:
    a = sum(circuitNodes[:edge[0]]) + edge[1]
    b = sum(circuitNodes[:edge[2]]) + edge[3]
    H[a, b] = 1

# second nodes are the midWires
for wire in range(2**n):
    start = circuitSize + wire * midWiresLength
    H[start, sum(circuitNodes[:wire+1]) - 1] = 1
    for i in range(midWiresLength - 1):
        H[start + i, start + i + 1] = 1
    end = start + midWiresLength - 1
    # connect to last site of dagger wire
    H[end, circuitSize + midWiresLength * 2**n + sum(daggerNodes[:wire+1]) - 1] = 1

# dagger wires
for edge in daggerEdges:
    a = sum(daggerNodes[:edge[0]]) + edge[1] + circuitSize + midWiresLength * 2**n
    b = sum(daggerNodes[:edge[2]]) + edge[3] + circuitSize + midWiresLength * 2**n
    # print([a, b, edge])
    H[a, b] = 1

# Add cirle wire to the 0th wire
start = circuitSize + midWiresLength * 2**n + daggerSize
H[start, circuitSize + midWiresLength * 2**n] = 1
for i in range(circleWireLength - 1):
    H[start + i, start + i + 1] = 1
H[start + circleWireLength - 1, 0] = 1

H = H + np.transpose(H)
removeTrails()
[w, v] = np.linalg.eigh(H)
G = nx.from_numpy_matrix(H)
nx.draw(G, with_labels=True)

tst = np.round(np.array([1/(np.arccos(val / 2) / np.pi) for val in w[np.abs(w) < 2]]), 4) - 4
f = np.where(tst == 0)
inds = np.array(list(range(8)) + list(range(16, 21)) + list(range(29, 45)) + list(range(57, 64)))
b = 1
plt.show()