from utils import rev_comp, pcov
import networkx as nx
with open("gasm_test_2", "r") as f:
    line = f.readline().strip()
    dnas = []
    while line:
        dnas.append(line)
        dnas.append(rev_comp(line))
        line = f.readline().strip()


    k = 2
    min1, min2, max1, max2 = 0, 0, 0, 0
    H = nx.Graph()

    while nx.number_connected_components(H) != 2 or min1 != 1 or max1 != 1 or min2 != 1 or max2 != 1:
        kmers = []
        G = nx.DiGraph()
        for dna in dnas:
            for i in range(len(dnas[0])-k+1):
                szo = dna[i:i+k]
                if szo not in kmers:
                    kmers.append(szo)
                    G.add_edge(szo[:-1], szo[1:])
        min1 = min([s[1] for s in G.in_degree])
        max1 = max([s[1] for s in G.in_degree])
        min2 = min([s[1] for s in G.out_degree])
        max2 = max([s[1] for s in G.out_degree])

        H = G.to_undirected()
        k += 1

    l = list(list(nx.connected_components(H))[0])
    print(pcov(l))




