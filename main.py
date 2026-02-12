import kagglehub
import matplotlib.pyplot as plt


# Download latest version
path = kagglehub.dataset_download("paultimothymooney/coronavirus-genome-sequence")
print("Path to dataset files:", path) #code to download dataset

import os
files = os.listdir(path)   #allows pycharm to read yo file
print(files)
def read_fasta(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    sequence = "".join([line.strip() for line in lines if not line.startswith(">")])
    return sequence
genome_path = path + r"\MN908947.fna"   #r avoids Windows path issues
genome = read_fasta(genome_path)

print("Genome length:", len(genome))
print("First 100 bases:", genome[:100])

def gc_content(genome):
    g = genome.count("G")
    c = genome.count("C")
    total = len(genome)
    return (g+c)/total * 100   #calcing GC% here
gc = gc_content(genome)
print("GC content of entire genome of covid 19:", gc)

window_size = 1000
windows = []
no_windows=0
for i in range(0,len(genome),window_size):
    window = genome[i:i+window_size]
    windows.append(window)
    no_windows+=1

gc_values = []
for w in windows:
    gc_values.append(gc_content(w))
print(gc_values)
print('no. of windows:',no_windows)

gc_rich = [i for i, x in enumerate(gc_values) if x > 40]
print("especially gc rich windows of COVID-19:",gc_rich)
mean_gcvalues = [x for x, x in enumerate(gc_values) if 35<x<39]
print("mean gc values of COVID-19:", mean_gcvalues)

# moving on to M. smegmatis genome

with open("data/GCA_000283295.1_ASM28329v1_genomic.fna", "r") as file:
    lines = file.readlines()
actino = "".join(line.strip() for line in lines if not line.startswith(">"))
def gc_content(actino):
    g = actino.count("G")
    c = actino.count("C")
    total = len(actino)
    return (g+c)/total * 100   #calcing GC% here
GC = gc_content(actino)
print("GC content of entire genome of M. smegmatis:", GC)

window_size = 1000 #trynna do window sliding calculation
windows = []
no_windows=0
for i in range(0,len(actino),window_size): #allows us to make each window like 0-1000 of genome, and next 1k-2k of genome etc
    window = actino[i:i+window_size]
    windows.append(window)
    no_windows+=1

GC_values = []  #making a list to store the gc vals of each window
for w in windows:
    GC_values.append(gc_content(w))
print('no. of windows:',no_windows)

GC_rich = [i for i, x in enumerate(GC_values) if x > 74]
print("especially gc rich windows of M. smegmatis:",GC_rich)

mean_GCvalues = [x for x, x in enumerate(GC_values) if 65<x<70]
print("mean gc vals of mycobacterium:", mean_GCvalues)

#for plasmodium falciparum
with open("data/GCA_000002765.3_GCA_000002765_genomic.fna", "r") as file:
    lines = file.readlines()
plasmodium = "".join(line.strip() for line in lines if not line.startswith(">"))
def gc_content(plasmodium):
    g = plasmodium.count("G")
    c = plasmodium.count("C")
    total = len(plasmodium)
    return (g+c)/total * 100   #calcing GC% here
gcc = gc_content(plasmodium)
print("GC content of entire genome of plasmodium falciparum:", gcc)

window_size = 1000 #trynna do window sliding calculation
windows = []
no_windows=0
for i in range(0,len(plasmodium),window_size): #allows us to make each window like 0-1000 of genome, and next 1k-2k of genome etc
    window = plasmodium[i:i+window_size]
    windows.append(window)
    no_windows+=1

gcc_values = []  #making a list to store the gc vals of each window
for w in windows:
    gcc_values.append(gc_content(w))
print('no. of windows:',no_windows)

mean_gccvalues = [x for x, x in enumerate(gcc_values) if 11<x<14]
print("mean gc values of plasmodium:", mean_gccvalues)

#graphing a line graph of gc_content in COVID-19 genome
plt.figure(figsize=(10,10))
plt.plot(gc_values, marker='o', linestyle='-', color='blue')
plt.title("GC content of entire genome of COVID-19")
plt.xlabel("window index")
plt.ylabel("GC Content(%)")
plt.grid(True)
plt.show()

#graph of a line graph of gc_content for Mycobacterium
plt.figure(figsize=(10,10))
plt.plot(GC_values, marker='o', linestyle='-', color='blue')
plt.title("GC content of entire genome of M.smegmatis")
plt.xlabel("window index")
plt.ylabel("GC Content(%)")
plt.grid(True)
plt.show()


plt.boxplot([mean_gcvalues, mean_GCvalues, mean_gccvalues])
plt.xticks([1,2,3], ["SARS-Cov-2", "Mycobacterium", "Plasmodium"])
plt.ylabel("GC content (%)")
plt.title("Comparing mean GC content values of the 3 selected genomes")
plt.grid(True)
plt.show()