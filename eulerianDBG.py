import networkx as nx
import matplotlib.pyplot as plt


def get_k_minus_one_mers(kmer_list):
    k_size = len(kmer_list[0])
    new_k = k_size - 1

    k_minus_one_list = []
    for kmer in kmer_list:
        k_minus_one_list.append(kmer[0:new_k])
        k_minus_one_list.append(kmer[k_size - new_k:])

    return k_minus_one_list


def get_edges(k_minus_one_list):
    edge_list = []
    for i in range(0, len(k_minus_one_list) - 1, 2):
        edge = (k_minus_one_list[i], k_minus_one_list[i + 1])
        edge_list.append(edge)

    return edge_list


def create_deBruijn(kmer_list):
    graph = nx.DiGraph()

    k_minus_one_nodes = get_k_minus_one_mers(kmer_list)
    k_minus_one_edges = get_edges(k_minus_one_nodes)

    graph.add_nodes_from(k_minus_one_nodes)
    graph.add_edges_from(k_minus_one_edges)

    return graph


def order_slice_visited(edge_list):
    positions = []
    for i in range(0, len(edge_list), 2):
        if edge_list[i][1] is not edge_list[i + 1][0]:
            pos = i
            positions.append(pos)

    if len(positions) == 0:
        return edge_list
    else:
        # Order the edges in initial eulerian walk fashion
        temp = edge_list[positions[0]]
        edge_list.remove(temp)
        edge_list.insert(positions[1], temp)

        # Slice the array for the other eulerian walk ordering
        before = edge_list[:positions[0] - 1]
        end = edge_list[positions[1]:]

        # Find where end and middle connect
        connection = 0
        for i in range(positions[0], positions[1]):
            if end[-1][1] == edge_list[i][0]:
                connection = i

        middle = edge_list[connection:positions[1]] + edge_list[positions[0] - 1:connection]

        sliced_edge_list = before + end + middle

        return edge_list, sliced_edge_list


def get_eulerian_walks(graph):
    # Adjacency matrix representation
    nodes_list = list(graph.nodes())
    adj_mat = nx.convert_matrix.to_numpy_array(graph, dtype=int)

    visited = []
    for i in range(adj_mat.shape[0]):
        for j in range(adj_mat.shape[1]):
            edge = (nodes_list[i], nodes_list[j])
            if adj_mat[i][j] == 1 and edge not in visited:
                visited.append(edge)

    ordered_edges = order_slice_visited(visited)
    walk1 = ordered_edges[0]
    walk2 = ordered_edges[1]
    return walk1, walk2


def get_sequence(walk):
    seq = walk[0][0] + walk[0][1][-1]

    for i in range(1, len(walk)):
        seq += walk[i][1][-1]
    return seq


def main():
    kmer_list = ['ATGT', 'TGTC', 'GTCT', 'TCTA', 'CTAG', 'TAGT', 'AGTG', 'GTGA', 'TGAA', 'GAAC',
                 'AACG', 'ACGT', 'CGTA', 'GTAG', 'TAGG', 'AGGC', 'GGCC', 'GCCT', 'CCTG', 'CTGA']

    graph = create_deBruijn(kmer_list)
    walk1, walk2 = get_eulerian_walks(graph)
    seq1 = get_sequence(walk1)
    seq2 = get_sequence(walk2)

    print('Assembly 1: {}'.format(seq1) + '\n' + 'Assembly 2: {}'.format(seq2))


if __name__ == '__main__':
    main()
