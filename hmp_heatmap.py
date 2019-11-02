import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns


hitsfile = 'HMPhits_redoscan.csv'

prot_to_idx = {'181.1':0, '198.1':1, '199.1':2, '200.1':3, '201.1':4, '202.1':5, '203.1':6, '205.1':7, '206.1':8,
               '207.1':9, '208.1':10, '394.1':11, '393.1':12, '392.1':13, '209.1':14, '210.1':15, '211.1':16, '212.1':17}

group_lines = {'G_DNA_Stool':[], 'mouth':[], 'G_DNA_Anterior nares':[], 'retroauricular':[], 'G_DNA_Posterior fornix':[]}
curr_counts = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

mouth_sites = ['G_DNA_Supragingival plaque', 'G_DNA_Tongue dorsum', 'G_DNA_Buccal mucosa', 'G_DNA_Subgingival plaque',
               'G_DNA_Palatine Tonsils', 'G_DNA_Saliva']
retroauricular_sites = ['G_DNA_L_Retroauricular crease', 'G_DNA_R_Retroauricular crease']
skip_sites = ['G_DNA_R_Antecubital fossa', 'G_DNA_Throat', 'G_DNA_Attached/Keratinized gingiva', 'G_DNA_Mid vagina',
              'G_DNA_Vaginal introitus', 'G_DNA_Hard palate', 'TEST', 'G_DNA_Posterior fornix']

with open(hitsfile, 'r') as infile:
    next(infile)
    for line in infile:
        line = line.strip().split(',')
        run = line[0]
        site = line[1]
        hits = [int(x) for x in line[2:]]

        if site in skip_sites:
            continue
        elif site in mouth_sites:
            group_lines['mouth'] = group_lines.get('mouth', []) + [hits]
        elif site in retroauricular_sites:
            continue
        else:
            group_lines[site] = group_lines.get(site, []) + [hits]

figure = []
for k,v in group_lines.items():
    if k == "Unlabeled":
        continue
    for run in sorted(v, key=sum, reverse=True)[0:501]:
        run = [x if x < 3 else 3 for x in run]
        figure.append(run)

colors = ['#DA5526', '#ED7D4A', '#ED934A', '#EFAD61', '#FFD093', '#F4F1E6']
colors = [x for x in colors[::-1]]
color_palette = sns.color_palette(colors)
color_map = ListedColormap(sns.color_palette(colors))

matplotlib.use("TkAgg")
ax = sns.heatmap(data=figure, cmap=color_map, xticklabels=prot_to_idx.keys())
plt.yticks(rotation=0)
plt.show(block=True)
