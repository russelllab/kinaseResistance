#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A script to create a circular bar plot.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import gzip

sys.path.append('../ML/')
import fetchData

def get_label_rotation(angle, offset):
    # Rotation must be specified in degrees :(
    rotation = np.rad2deg(angle + offset)
    if angle <= np.pi/2 or angle >= 3*np.pi/2:
        alignment = "left"
        # rotation = rotation + 180
    else: 
        alignment = "right"
        rotation = rotation + 180
    return rotation, alignment

def add_labels(angles, values, labels, offset, ax):
    # This is the space between the end of the bar and the label
    padding = 0.5
    
    # Iterate over angles, values, and labels, to add all of them.
    for angle, value, label, in zip(angles, values, labels):
        angle = angle
        
        # Obtain text rotation and alignment
        rotation, alignment = get_label_rotation(angle, offset)

        # And finally add the text
        ax.text(
            x=angle, 
            y=value + padding, 
            s=label, 
            ha=alignment, 
            va="center", 
            size=30,
            rotation=rotation, 
            rotation_mode="anchor"
        ) 


# Ensures reproducibility of random numbers
rng = np.random.default_rng(123)

df = pd.DataFrame({
    "name": [f"item {i}" for i in range(1, 51)],
    "value": rng.integers(low=10, high=50, size=50),
    "group": ["A"] * 10 + ["B"] * 20 + ["C"] * 12 + ["D"] * 8
})

# Define a dictionary to store a region's start and end positions
# dic_region -> {region_name: [start, end, order_num]}
# where order_num is the order in which the region appears in the file
dic_region = {}
order_num = 0
for line in gzip.open('../alignments/humanKinasesHitsSplitTrimmed_ss.tsv.gz', 'rt', encoding='utf-8'):
    if line.startswith('#'): continue
    if line.split()[0] in ['DFG-motif', 'HrD-motif', 'APE-motif',
                            'Catalytic-Lys', 'Gly-rich-loop']: continue
    order_num += 1
    name = line.split()[0]
    start, end = line.split()[1].strip().split('-')
    dic_region[name] = [int(start), int(end), order_num]

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

# Define a dictionary to store the pfam positions and the corresponding amino acids
# dif_pfam2aa -> {pfampos: pfamaa}
dif_pfam2aa = {}
mycursor.execute("SELECT pfampos, pfamaa FROM hmm")
hits = mycursor.fetchall()
for hit in hits:
    pfampos, pfamaa = hit
    if pfampos == '-': continue
    pfampos = int(pfampos)
    dif_pfam2aa[pfampos] = pfamaa

# Define a dictionary to store the pfam positions and number of corresponding mutations
# at that position
# dic_pfam -> {pfampos: {mut_type: num}}

# select all mutations from the database
mycursor.execute("SELECT mutation, mut_type, pfampos, acc, gene FROM mutations\
                where pfampos!='%s'" % ('-', ))
hits = mycursor.fetchall()
dic_pfam = {}
for hit in hits:
    mutation, mut_type, pfampos, acc, gene = hit
    pfampos = int(pfampos)
    if pfampos not in dic_pfam:
        dic_pfam[pfampos] = {'constitutive-activation':0, 'increase':0, 'loss':0, 'decrease':0, 'resistance':0, 'neutral': 0}
    if mut_type in ['activating']:
        dic_pfam[pfampos]['constitutive-activation'] += 1
    else:
        # if str(pfampos) == '127' and mut_type == 'loss':
        #     dic_pfam[pfampos][mut_type] = 100
        #     continue
        dic_pfam[pfampos][mut_type] += 1
    
# for pfampos in dic_pfam:
#     for mut_type in dic_pfam[pfampos]:
#         dic_pfam[pfampos][mut_type] = np.log2(dic_pfam[pfampos][mut_type] + 1)*2

# Define a dictionary to store the pfam positions and PTM information
# dic_pfam2psite -> {pfampos: {ptmsite: num}}
dic_pfam2ptmsite = {}
mycursor.execute("SELECT ptmtype, pfampos FROM ptms\
                where pfampos!='%s'" % ('-', ))
hits = mycursor.fetchall()
for hit in hits:
    ptmtype, pfampos = hit
    pfampos = int(pfampos)
    if pfampos not in dic_pfam2ptmsite:
        dic_pfam2ptmsite[pfampos] = {}
    if ptmtype not in dic_pfam2ptmsite[pfampos]:
        dic_pfam2ptmsite[pfampos][ptmtype] = 0
    dic_pfam2ptmsite[pfampos][ptmtype] += 1

def generate_shades_of_cyan(number):
    number = int(number)
    number_log2 = np.log2(number)
    number_log2 = int(number_log2)
    # number = number / 10
    # if number < 1 or number > 10:
    #     raise ValueError("Number must be between 1 and 10.")

    # Determine the step size to create shades
    step = 255 // 10

    # Calculate the RGB values for the given number
    red = 0
    green = 255 - (number_log2 * step)
    blue = 255

    # print (red, green, blue)

    # Convert RGB values to hexadecimal color code
    hex_color = f"#{red:02x}{green:02x}{blue:02x}"
    return hex_color

dic_ptm_colors = {}
for pfampos in dic_pfam2ptmsite:
    for ptmtype in dic_pfam2ptmsite[pfampos]:
        if dic_pfam2ptmsite[pfampos][ptmtype] >= 20:
            print (pfampos, ptmtype, dic_pfam2ptmsite[pfampos][ptmtype])
            # dic_ptm_colors[pfampos] = generate_shades_of_cyan(dic_pfam2ptmsite[pfampos][ptmtype])
            dic_ptm_colors[pfampos] = 'cyan' if ptmtype == 'p' else 'grey'

print (dic_ptm_colors)
# sys.exit()

def make_fractions(mut_dic, total):
    total_log2 = np.log2(total + 1)*2
    for mut_type in mut_dic:
        mut_type_frac = float(mut_dic[mut_type])/total
        mut_type_frac_new = mut_type_frac*total_log2
        mut_dic[mut_type] = mut_type_frac_new
    return total_log2

num = 0
unique_pfam = {}
name = []
value = []
group = []
order = []
const_act = []
inc = []
loss = []
dec = []
resistance = []
neutral = []
ptm = []
all_var = []
# Now we construct the dataframe. But before that we save the data in lists
# so that we can sort them according to the order in which the regions appear
# in the file. We will create lists for the following:
# name, value, group, order, const_act, inc, loss, dec, resistance, neutral, all_var

# for pfampos in dic_pfam:
for pfampos in range(1, 882):
    total = 0
    # Calculate the total number of mutations at that pfam position
    if pfampos in dic_pfam:
        for mut_type in dic_pfam[pfampos]:
            total += dic_pfam[pfampos][mut_type]
    # else:
    #     dic_pfam[pfampos] = {'constitutive-activation':0, 'increase':0, 'loss':0, 'decrease':0, 'resistance':0, 'neutral': 0}
    if total > 0:
        total_log2 = make_fractions(dic_pfam[pfampos], total)
        num += 1
        for region in dic_region:
            if dic_region[region][0] <= int(pfampos) <= dic_region[region][1]:
                if region not in unique_pfam:
                    unique_pfam[region] = []
                unique_pfam[region].append(pfampos)
                name.append(pfampos)
                value.append(total)
                group.append(region)
                order.append(dic_region[region][2])
                const_act.append(dic_pfam[pfampos]['constitutive-activation'])
                inc.append(dic_pfam[pfampos]['increase'])
                loss.append(dic_pfam[pfampos]['loss'])
                dec.append(dic_pfam[pfampos]['decrease'])
                resistance.append(dic_pfam[pfampos]['resistance'])
                neutral.append(dic_pfam[pfampos]['neutral'])
                ptm.append(1)
                all_var.append(total_log2)
                break
        print (pfampos, dic_pfam[pfampos])
# print(df)
# print (len(dec))
# sys.exit()
print (num)
# unique_pfam = list(set(unique_pfam))
# print (unique_pfam)

# Now we create the dataframe
df = pd.DataFrame({
    "name": name,
    "value": value,
    "constitutive-activation": const_act,
    "increase": inc,
    "loss": loss,
    "decrease": dec,
    "resistance": resistance,
    "neutral": neutral,
    "ptm": ptm,
    'all_var': all_var,
    "group": group,
    "order": order
})
df = df.sort_values(by=['order', 'name'])
GROUP_NAMES = []
for region in df["group"].values:
    if region not in GROUP_NAMES:
        GROUP_NAMES.append(region)

print (df)
# sys.exit()


# Determines where to place the first bar. 
# By default, matplotlib starts at 0 (the first bar is horizontal)
# but here we say we want to start at pi/2 (90 deg)
OFFSET = np.pi / 2
OFFSET = 0

# All this part is like the code above
VALUES = df["constitutive-activation"].values
TOTAL_LOG2 = df['all_var'].values
LABELS = df["name"].values
NEW_LABELS = []
for i in range(len(LABELS)):
    pfampos = LABELS[i]
    # NEW_LABELS += [str(dif_pfam2aa[pfampos]) + ' (' + str(pfampos) + ')']
    NEW_LABELS += [str(dif_pfam2aa[pfampos])]
    # print (LABELS[i], dif_pfam2aa[LABELS[i]])
    # NEW_LABELS += dif_pfam2aa[LABELS[i]]
LABELS = NEW_LABELS
GROUP = df["group"].values

PAD = 1
ANGLES_N = 1 + len(VALUES) + PAD * len(np.unique(GROUP))
ANGLES = np.linspace(0, (2 * np.pi), num=ANGLES_N, endpoint=False)
WIDTH = (2 * np.pi) / len(ANGLES)

offset = 0
IDXS = []
# GROUPS_SIZE = [10, 20, 12, 8]
# GROUPS_SIZE = [len(group_name) for group_name in GROUP_NAMES]
GROUPS_SIZE = [len(unique_pfam[group_name]) for group_name in GROUP_NAMES]
print (GROUPS_SIZE)
for size in GROUPS_SIZE:
    IDXS += list(range(offset + PAD, offset + size + PAD))
    offset += size + PAD

print(ANGLES, len(ANGLES))
print(IDXS, len(IDXS))
print (VALUES)
bottom = [0 for i in range(0, len(VALUES))]
bottom = np.array(bottom)
# sys.exit()

fig, ax = plt.subplots(figsize=(100, 100), subplot_kw={"projection": "polar"})
ax.set_theta_offset(OFFSET)
ax.set_ylim(-50, 50)
ax.set_frame_on(False)
ax.xaxis.grid(False)
ax.yaxis.grid(False)
ax.set_xticks([])
ax.set_yticks([])

# GROUPS_SIZE = [10, 20, 12, 8]
COLORS = [f"C{i}" for i, size in enumerate(GROUPS_SIZE) for _ in range(size)]
COLORS = ['red' for i in range(0, len(VALUES))]

# for i in range(0, 2):
for col, arr in zip(['-', 'green', 'lightgreen', 'red', 'coral', 'blue', '#F2E34C'],
                    [
                    df["ptm"].values,
                    df["constitutive-activation"].values,
                    df["increase"].values,
                    df["loss"].values,
                    df["decrease"].values,
                    df["resistance"].values,
                    df["neutral"].values
                    ]):
    if col == '-':
        COLORS = [dic_ptm_colors[pfampos] if pfampos in dic_ptm_colors else 'white' for pfampos in df['name'].values]
    else:
        COLORS = [col for i in range(0, len(VALUES))]
    ax.bar(
        ANGLES[IDXS], arr, width=WIDTH, color=COLORS, 
        edgecolor="white", linewidth=2, bottom=bottom
    )
    bottom = bottom + arr

# add_labels(ANGLES[IDXS], VALUES, LABELS, OFFSET, ax)
add_labels(ANGLES[IDXS], bottom, LABELS, OFFSET, ax)

# Extra customization below here --------------------

# This iterates over the sizes of the groups adding reference
# lines and annotations.

offset = 0
N = 10
# for group, size in zip(["A", "B", "C", "D"], GROUPS_SIZE):
for group, size in zip(GROUP_NAMES, GROUPS_SIZE):
    # Add line below bars
    x1 = np.linspace(ANGLES[offset + PAD], ANGLES[offset + size + PAD - 1], num=N)
    ax.plot(x1, [-2] * N, color="#333333")
    
    # Add text to indicate group
    ax.text(
        np.mean(x1), -12, group, color="#333333", fontsize=50, 
        fontweight="bold", ha="center", va="center"
    )
    
    # Add reference lines at 20, 40, 60, and 80
    x2 = np.linspace(ANGLES[offset], ANGLES[offset + PAD - 1], num=N)
    # ax.plot(x2, [3] * N, color="#bebebe", lw=0.8)
    # ax.plot(x2, [6] * N, color="#bebebe", lw=0.8)
    # ax.plot(x2, [9] * N, color="#bebebe", lw=0.8)
    # ax.plot(x2, [8] * N, color="#bebebe", lw=0.8)
    # ax.plot(x2, [10] * N, color="#bebebe", lw=0.8)
    
    offset += size + PAD

# plt.show()
plt.savefig('circular_bar_plot.png', format='png', dpi=400)
plt.savefig('circular_bar_plot.svg', format='svg', dpi=400)
# plt.savefig('circular_bar_plot_all_positions.png', format='png', dpi=400)
# plt.savefig('circular_bar_plot_all_positions.svg', format='svg', dpi=400)