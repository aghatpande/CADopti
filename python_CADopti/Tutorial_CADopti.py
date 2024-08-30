import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from CADopti import CADopti
from assembly_assignment_matrix import assembly_assignment_matrix
from assembly_activity_function import assembly_activity_function

# Load data
data = loadmat('Data.mat')
spM = data['spM']
nneu = spM.shape[0]  # number of recorded units

# Detection
BinSizes = [0.015, 0.025, 0.04, 0.06, 0.085, 0.15, 0.25, 0.4, 0.6, 0.85, 1.5]
MaxLags = [10] * len(BinSizes)

# Assembly detection
As_across_bins, As_across_bins_index, assembly = CADopti(spM, MaxLags, BinSizes)

# Visualization
display = 'raw'  # or 'clustered'
Amatrix, Binvector, Unit_order, As_order = assembly_assignment_matrix(As_across_bins, nneu, BinSizes, display)

plt.figure()
plt.imshow(Amatrix)
plt.show()

# Assembly Activation
lagChoice = 'duration'
act_count = 'full'
assembly_activity = assembly_activity_function(As_across_bins, assembly, spM, BinSizes, lagChoice, act_count)

plt.figure(figsize=(10, 10))
for i, activity in enumerate(assembly_activity):
    plt.subplot(len(assembly_activity), 1, i+1)
    plt.plot(activity[:, 0], activity[:, 1])
plt.tight_layout()
plt.show()