import copy
import decors as dec


@dec.timer
def kmer_dict2graph(kmers: dict) -> dict:
    graph = dict()

    for kmer, count in kmers.items():
        node_from = kmer[:-1]
        node_to = kmer[1:]
        if count > 0:
            graph.setdefault(node_from, list())
            graph.setdefault(node_to, list())
            graph[node_from].extend([node_to] * count)

    return graph


@dec.timer
def eulerian_path(graph: dict) -> list:
    path = list()

    # Make a mutable copy of the graph
    g = {node: [neighbor for neighbor in neighbors] for node, neighbors in graph.items()}

    # Calculate in-degrees and out-degrees
    out_deg = {node: len(neighbors) for node, neighbors in graph.items()}
    in_deg = {node: 0 for node in graph}
    for neighbors in graph.values():
        for neighbor in neighbors:
            in_deg[neighbor] += 1

    # Find start
    start = next((node for node in graph if out_deg[node] - in_deg[node] == 1),
                 next((node for node in graph if out_deg[node] > 0), None))

    stack = [start]

    while stack:
        curr_node = stack[-1]
        if g[curr_node]:
            next_node = g[curr_node].pop()
            stack.append(next_node)
        else:
            path.append(stack.pop())

    return path[::-1]


def all_eulerian_paths(graph: dict) -> list:
    paths = set()
    stack = list()

    # Make a mutable copy of the graph
    g = {node: [neighbor for neighbor in neighbors] for node, neighbors in graph.items()}

    # Calculate in-degrees and out-degrees
    out_deg = {node: len(neighbors) for node, neighbors in graph.items()}
    in_deg = {node: 0 for node in graph}
    for neighbors in graph.values():
        for neighbor in neighbors:
            in_deg[neighbor] += 1

    # Find start
    start = next((node for node in graph if out_deg[node] - in_deg[node] == 1),
                 next((node for node in graph if out_deg[node] > 0), None))

    node_stack = [start]
    path = list()

    stack.append((g, node_stack, path))

    while stack:
        curr_g, curr_node_stack, curr_path = stack[-1]

        # doing eulpath on curr_g
        while curr_node_stack:
            curr_node = curr_node_stack[-1]
            if curr_g[curr_node]:
                neighbors_left = set(curr_g[curr_node])  # ignore duplicates
                next_node = curr_g[curr_node].pop()
                curr_node_stack.append(next_node)
                neighbors_left.remove(next_node)

                # for every neighbor that might have been picked as next node, a new graph is added to the stack where that neighbor is picked
                if neighbors_left:
                    for neighbor in neighbors_left:
                        new_g = copy.deepcopy(curr_g)
                        new_g[curr_node].append(next_node)
                        new_g[curr_node].remove(neighbor)
                        new_node_stack = curr_node_stack.copy()
                        new_node_stack.pop()
                        new_node_stack.append(neighbor)
                        stack.append((new_g, new_node_stack, curr_path.copy()))
            else:
                curr_path.append(curr_node_stack.pop())

        # once a graph is exhausted, add the path to paths and remove it from the stack
        paths.add(condense(curr_path[::-1]))
        stack.pop()

    return list(paths)


@dec.timer
def condense(path: list) -> str:
    sequence = list()

    for node in path:
        sequence.append(node[0])
    sequence.append(path[-1][1:])

    return "".join(sequence)
