import numba
import numpy as np
from numba import types
from numba.typed import Dict


@numba.njit
def non_recursive_dfs(graph, start, visited):
    stack = []
    stack.append(start)
    component = []

    while len(stack) > 0:
        node = stack.pop()

        if not visited[node]:
            visited[node] = True
            component.append(node)
            # Push neighbors to the stack
            for neighbor in graph[node]:
                if not visited[neighbor]:
                    stack.append(neighbor)

    return component


@numba.njit
def find_connected_components(graph, num_nodes):
    visited = np.zeros(num_nodes, dtype=np.bool_)
    components = []

    for node in range(num_nodes):
        if not visited[node]:
            component = non_recursive_dfs(graph, node, visited)
            components.append(component)

    return components


# Example usage with numba.typed.Dict
def create_graph():
    graph = Dict.empty(key_type=types.int32, value_type=types.ListType(types.int32))

    graph[types.int32[0]] = [1]
    graph[1] = [0, 2]
    graph[2] = [1]
    graph[3] = [4]
    graph[4] = [3]

    return graph


graph = Dict.empty(key_type=types.int32, value_type=types.ListType(types.int32))


# Create the graph
graph = create_graph()
num_nodes = len(graph)

# Find connected components
components = find_connected_components(graph, num_nodes)
print("Connected Components:", components)
