import sys
import networkx as nx
import matplotlib.pyplot as plt

def compute_energy_matrix(sequence, alpha):
    length = len(sequence)
    energy_matrix = [[0] * length for _ in range(length)] 

    for gap in range(1, length):
        for i in range(length - gap):
            j = i + gap
            min_energy = float('inf')
            min_energy = min(min_energy, energy_matrix[i + 1][j])
            min_energy = min(min_energy, energy_matrix[i][j - 1])
            if sequence[i] + sequence[j] in ['CG', 'GC', 'AU', 'UA', 'GU', 'UG']:
                min_energy = min(min_energy, energy_matrix[i + 1][j - 1] + alpha[sequence[i] + sequence[j]])
            for k in range(i + 1, j):
                min_energy = min(min_energy, energy_matrix[i][k - 1] + energy_matrix[k][j])

            energy_matrix[i][j] = min_energy

    return energy_matrix

def predict_secondary_structure(sequence, alpha):
    length = len(sequence)
    energy_matrix = compute_energy_matrix(sequence, alpha)
    score = energy_matrix[0][length - 1]

    return energy_matrix, score

def traceback(sequence, energy_matrix, alpha):
    length = len(sequence)
    traceback_pairs = ["."] * length 
    paired_positions = []
    structures = []

    def traceback_helper(i, j):
        nonlocal traceback_pairs, paired_positions, structures

        if i >= j:
            return
        
        if sequence[i] + sequence[j] in ['CG', 'GC', 'AU', 'UA', 'GU', 'UG'] and energy_matrix[i][j] == energy_matrix[i + 1][j - 1] + alpha[sequence[i] + sequence[j]]:
            traceback_pairs[i] = "("
            traceback_pairs[j] = ")"
            paired_positions.append((i + 1, j + 1))  
            traceback_helper(i + 1, j - 1)
            return

        if energy_matrix[i][j] == energy_matrix[i + 1][j]:
            traceback_helper(i + 1, j)
            return
    
        if energy_matrix[i][j] == energy_matrix[i][j - 1]:
            traceback_helper(i, j - 1)
            return
        
        for k in range(i + 1, j):
            if energy_matrix[i][j] == energy_matrix[i][k - 1] + energy_matrix[k][j]:
                traceback_helper(i, k - 1)
                traceback_helper(k, j)
                return
    
    traceback_helper(0, length - 1)
    
    # Identificar estructuras: bucles internos, tallos, protuberancias
    visited = [False] * length
    for i, pair in enumerate(traceback_pairs):
        if pair == "(" and not visited[i]:
            j = traceback_pairs.index(")", i + 1)
            visited[i] = visited[j] = True
            if j - i > 1:
                if all(traceback_pairs[x] == "." for x in range(i + 1, j)):
                    structures.append(f"Bucle interno: {i + 1}-{j + 1}")
                else:
                    structures.append(f"Tallo: {i + 1}-{j + 1}")
            else:
                if i == 0 or j == length - 1:
                    structures.append(f"Base no emparejada: {i + 1}-{j + 1}")
                else:
                    structures.append(f"Bulbo: {i + 1}-{j + 1}")

    return traceback_pairs, paired_positions, structures

def plot_structure(sequence, pairs, filename, node_color="#808080", edge_color="#808080", highlight_color="#FF5733"): 
    plt.figure(figsize=(max(4, len(sequence) / 8), max(4, len(sequence) / 8)))

    G = nx.Graph()
    for idx, base in enumerate(sequence):
        G.add_node(idx + 1, base=base, color=node_color)

    for idx in range(1, len(sequence)):
        G.add_edge(idx, idx + 1, color=edge_color, style="-")

    for pair in pairs:
        G.add_edge(pair[0], pair[1], color=highlight_color, style="--")

    pos = nx.kamada_kawai_layout(G)

    nx.draw_networkx_nodes(G, pos, node_color=node_color, node_size=200, edgecolors=node_color, linewidths=1.5)
    nx.draw_networkx_labels(G, pos, labels=nx.get_node_attributes(G, "base"), font_color="black")

    edges = G.edges()
    edge_colors = [G[u][v]['color'] for u, v in edges]
    edge_styles = [G[u][v]['style'] for u, v in edges]
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=edge_colors, style=edge_styles, width=3)

    plt.axis('off')
    plt.savefig(filename)  # Guarda la figura con el nombre especificado
    plt.close()

sequence_input = "ucaagcguuagagaagucauuaugugauaaaaaaauucaacuugguaucaacuuaacuaa"

sequence = sequence_input.upper()
alpha_dict = {"CG": -1, "GC": -1, "AU": -1, "UA": -1, "GU": -1, "UG": -1}

energy_matrix, min_energy_score = predict_secondary_structure(sequence, alpha_dict)
traceback_structure, paired_positions, structures = traceback(sequence, energy_matrix, alpha_dict)

# Guardar la estructura en un archivo con el nombre 'structure.png'
plot_structure(sequence, paired_positions, filename='structure3.png', node_color="#808080", edge_color="#808080", highlight_color="#FF5733")

# Guardar el puntaje de energía, las estructuras y los pares emparejados en un archivo de salida
with open('output3.txt', 'w') as f:
    f.write(f"Puntaje mínimo de energía: {min_energy_score}\n\n")
    
    f.write("Matriz de energía:\n")
    for row in energy_matrix:
        f.write("\t".join(map(str, row)) + "\n")
    
    f.write("\nPosiciones emparejadas:\n")
    for pair in paired_positions:
        f.write(f"{pair[0]}-{pair[1]}\n")
    
    f.write("\nEstructuras:\n")
    for struct in structures:
        f.write(f"{struct}\n")
