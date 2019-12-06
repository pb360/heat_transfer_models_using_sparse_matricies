# command line usage: 
# $ python3 postprocess.py input1.txt solution<#>.txt
import sys 
import math
import matplotlib.pyplot as plt
import numpy as np
import time

if len(sys.argv) != 3:
    print('Usage:')
    print('  $ python3 postprocess.py input1.txt solution<#>.txt')
    sys.exit(0)

# read in input file 
input = sys.argv[1]
f = open(input, 'r')
line1 = f.readline()
line2 = f.readline()
f.close()
line1_list = line1.split(' ')
line2_list = line2.split(' ')

# get values from input file. 
length = float(line1_list[0])
width  = float(line1_list[1])
h      = float(line1_list[2])
T_c    = float(line2_list[0])
T_h    = float(line2_list[1])

# read in solution 
solution = sys.argv[2]
f = open(solution, 'r')
line = f.readline()
f.close()
soln_list = line.split(' ')
del soln_list[-1] # adjustment for last space being counted as an item
soln_list = [float(i) for i in soln_list]

# function to generate bottom row of heat map 
def get_cold_temp(j): # left out function parameters: , int ncols, float T_c
    x = j/ncols*length
    exp_term = x-length/2
    
    cold = -T_c*(math.exp(-10*(exp_term)**2)-2)
    return cold

# make heatmap matrix 
nrows = int(width/h  + 1)
ncols = int(length/h + 1)
heat_matrix = np.zeros((nrows, ncols), dtype = float)
sol_index = 0 # for walking through solution list

for j in range(ncols): 
    for i in range(nrows): 
        # top row set to T_h 
        if i == 0: 
            heat_matrix[i][j] = T_h

        # set bottom row to T_c function 
        elif i == nrows - 1:
            heat_matrix[i][j] = get_cold_temp(j)

        # fill inside from solution vector 
        elif j != ncols-1: 
            heat_matrix[i][j] = soln_list[sol_index]
            sol_index += 1 
            # print("i: " + str(i) + ",  j: " + str(j), + ",  j: " + str(j))

        else: # last column is periodic boundary (set equal to first column)
            heat_matrix[i][j] = heat_matrix[i][0]

def heatmap2d(arr: np.ndarray):
    # fig, ax = plt.imshow(arr, cmap='viridis')
    plt.imshow(arr, cmap='viridis')
    plt.colorbar()
    plt.title(solution)
    xlocs, labels = plt.xticks() 
    xlabs_len = len(xlocs)-2
    ylocs, labels = plt.yticks() 
    ylabs_len = len(ylocs)-2

    # x axis labels 
    xlabs = list(np.linspace(0,length, xlabs_len))
    xlabs = [round(i, 3) for i in xlabs]
    xtics = list(np.linspace(0,ncols, xlabs_len))

    # y axis labels
    ylabs = list(np.linspace(width, 0, ylabs_len))
    ylabs = [round(i, 2) for i in ylabs]
    ytics = list(np.linspace(0,nrows, ylabs_len))
    plt.xticks(xtics, xlabs)
    plt.yticks(ytics, ylabs)
    plt.show()

heatmap2d(heat_matrix)

print("post process ran fully")



