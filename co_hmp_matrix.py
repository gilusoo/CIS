import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

'''First create a dictionary in the format {run: {<protein>:<count> ...} ...} '''

runs_dict = {}  # Large dictionary with format { <run>, { <protein> : <count> ...} ...}
pc = {}
with open('HMPhits_redoscan.csv', 'r') as infile:
    next(infile)
    for line in infile:
        line = line.strip().split(',')
        currun = line[0]
        pc['181.1'] = int(line[2])
        pc['198.1'] = int(line[3])
        pc['199.1'] = int(line[4])
        pc['200.1'] = int(line[5])
        pc['201.1'] = int(line[6])
        pc['202.1'] = int(line[7])
        pc['203.1'] = int(line[8])
        pc['205.1'] = int(line[9])
        pc['206.1'] = int(line[10])
        pc['207.1'] = int(line[11])
        pc['208.1'] = int(line[12])
        pc['394.1'] = int(line[13])
        pc['393.1'] = int(line[14])
        pc['392.1'] = int(line[15])
        pc['209.1'] = int(line[16])
        pc['210.1'] = int(line[17])
        pc['211.1'] = int(line[18])
        pc['212.1'] = int(line[19])
        runs_dict[currun] = pc
        pc = {}

'''Next, create a pandas dataframe using the dictionary created above.
    Proteins should be the rows and columns and where they intersect should be CO count.'''

proteins = ['181.1', '198.1', '199.1', '200.1', '201.1', '202.1', '203.1',
            '205.1', '206.1', '207.1', '208.1', '394.1', '393.1', '392.1', '209.1',
            '210.1', '211.1', '212.1']
df = pd.DataFrame(columns=proteins, index=proteins)
df[:] = int(0)

'''The following code creates a list of total number of hits for each protein across all runs.'''

hits = []
for protein in proteins:
    count = 0
    for value in runs_dict.values():
        for prot in value.keys():
            if prot == protein and value[prot] != 0:
                count += value[prot]
    hits.append(count)

'''Adjusting total hits.'''

hits = [x/5 for x in hits]

for key, value in runs_dict.items():                    # Edited for weighted CO. Original in CIS/coMatrix.py
    present = {}
    for prot, count in value.items():
        if count != 0:
            present[prot] = count
    # print(present)
    for protein1 in proteins:
        for protein2 in proteins:
            if protein1 in present and protein2 in present:
                if present[protein1] < present[protein2]:
                    df[protein1][protein2] += present[protein1]
                    df[protein2][protein1] += present[protein1]
                else:
                    df[protein1][protein2] += present[protein2]
                    df[protein2][protein1] += present[protein2]

edge_list = []
for index, row in df.iterrows():
    i = 0
    for col in row:
        weight = float(col)/500         # Change 100 to any max; line weight will be relative to this number
        edge_list.append((index, df.columns[i], weight))
        i += 1
# Remove lines if there is no CO
updated_edge_list = [x for x in edge_list if not x[2] == 0.0]

node_list = []
for i in proteins:
    for e in updated_edge_list:
        if i == e[0] and i == e[1]:
            node_list.append(i)
for i in node_list:
    if i[1] == 0.0:
        node_list.remove(i)

for i in updated_edge_list:
    if i[0] == i[1]:
        updated_edge_list.remove(i)

# Set canvas size
plt.subplots(figsize=(10,10))

# Networkx graph
G = nx.Graph()
for i in sorted(node_list):
    G.add_node(i[0], size = i[1])
G.remove_node('1')
G.remove_node('2')
G.remove_node('3')
G.add_weighted_edges_from(updated_edge_list)

node_order = ['181.1', '198.1', '199.1', '200.1', '201.1', '202.1', '203.1',
            '205.1', '206.1', '207.1', '208.1', '394.1', '393.1', '392.1', '209.1',
            '210.1', '211.1', '212.1']

# Reorder node list
updated_node_order = []
for i in node_order:
    for x in node_list:
        if x[0] == i:
            updated_node_order.append(x)

# Reorder edge list
test = nx.get_edge_attributes(G, 'weight')
updated_again_edges = []
for i in nx.edges(G):
    for x in test.keys():
        if i[0] == x[0] and i[1] == x[1]:
            updated_again_edges.append(test[x])
# print(updated_again_edges)


# Drawing customization
node_scalar = 80
edge_scalar = 100
sizes = [x[1] for x in node_list]
widths = [x for x in updated_again_edges]

# Draw the graph!!!
pos = nx.spring_layout(G)

nx.draw(G, pos, with_labels=True, font_size=14, font_weight='bold', width = widths, node_size=hits)
plt.savefig("COmatrix_weighted.png")
plt.show()