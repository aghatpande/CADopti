import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from CADopti import CADopti
from assembly_assignment_matrix import assembly_assignment_matrix
from assembly_activity_function import assembly_activity_function

if __name__ == '__main__':
    # Load data
    mat_data = sio.loadmat('Data.mat')

    # Print the keys in the mat file for troubleshooting if needed
    # print("Keys in Data.mat:", mat_data.keys()) # uncomment for troubleshooting

    # Print information about spM if needed for troubleshooting loading the data
    # print("spM shape:", mat_data['spM'].shape)
    # print("spM data type:", mat_data['spM'].dtype)
    # print("Sample of spM (first 5 rows, first 10 columns):")
    # print(mat_data['spM'][:5, :10])

    # Assign spM correctly
    spM = mat_data['spM']

    # Convert spM to a list of spike times, removing NaNs
    spike_times = [row[~np.isnan(row)] for row in spM]

    # Print some statistics to check if the data is correct
    print("Number of neurons:", len(spike_times))
    print("Number of spikes per neuron:")
    print([len(spikes) for spikes in spike_times])

    try:
        # Define parameters
        MaxLags = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10] # 10 is the maximum lag for each bin size
        BinSizes = [0.015, 0.025, 0.04, 0.06, 0.085, 0.15, 0.25, 0.4, 0.6, 0.85, 1.5]

        # Call CADopti function
        As_across_bins, As_across_bins_index, assembly = CADopti(spM, MaxLags, BinSizes)

        if As_across_bins is None:
            raise ValueError("No valid assemblies found. Check the input data and parameters.")

        # Visualization
        nneu = len(spike_times)  # number of recorded units
        display = 'raw'  # or 'clustered'
        Amatrix, Binvector, Unit_order, As_order = assembly_assignment_matrix(As_across_bins, nneu, BinSizes, display)

        plt.figure()
        plt.imshow(Amatrix)
        plt.title('Assembly Assignment Matrix')
        plt.xlabel('Assembly')
        plt.ylabel('Neuron')
        plt.colorbar(label='Assignment')
        plt.savefig('assembly_assignment_matrix.png')
        plt.close()

        # Assembly Activation
        lagChoice = 'duration'
        act_count = 'full'
        assembly_activity = assembly_activity_function(As_across_bins, assembly, spike_times, BinSizes, lagChoice, act_count)

        plt.figure(figsize=(10, 10))
        for i, activity in enumerate(assembly_activity):
            plt.subplot(len(assembly_activity), 1, i+1)
            plt.plot(activity[:, 0], activity[:, 1])
            plt.title(f'Assembly {i+1} Activity')
            plt.xlabel('Time')
            plt.ylabel('Activity')
        plt.tight_layout()
        plt.savefig('assembly_activity.png')
        plt.close()

        # Print summary of detected assemblies
        print("\nDetected Assemblies:")
        for i, assembly in enumerate(As_across_bins):
            print(f"Assembly {i+1}:")
            print(f"  Elements: {assembly['elements']}")
            print(f"  Bin size: {assembly['bin']}")
            print(f"  Lags: {assembly['lag']}")
            print(f"  p-values: {assembly['pr']}")
            print(f"  Occurrences: {assembly['Noccurrences']}")
            print()
    except Exception as e:
        print(f"An error occurred: {e}")