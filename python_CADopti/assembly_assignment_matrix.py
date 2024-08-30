import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform

def assembly_assignment_matrix(As_across_bins, nneu, BinSizes, display='raw'):
    """
    Display the Assembly assignment matrix.

    Args:
    As_across_bins (list): List of dictionaries containing assembly information
    nneu (int): Number of units
    BinSizes (list): List of investigated bin sizes
    display (str): Display options: 'raw', 'ordunit', or 'clustered'

    Returns:
    tuple: Amatrix, Binvector, Unit_order, As_order
    """
    if display not in ['raw', 'ordunit', 'clustered']:
        print("'raw', 'ordunit', 'clustered' are the only acceptable style specifications")
        return

    nAss_final = len(As_across_bins)
    AAT = np.full((nneu, nAss_final), np.nan)
    Binvector = np.full(nAss_final, np.nan)

    for i in range(nAss_final):
        AAT[As_across_bins[i]['elements'], i] = As_across_bins[i]['lag']
        Binvector[i] = As_across_bins[i]['bin']

    if display == 'raw':
        Unit_order = np.arange(1, nneu + 1)
        As_order = np.arange(1, len(As_across_bins) + 1)

    elif display == 'ordunit':
        aus = np.all(np.isnan(AAT), axis=1)
        idx_nan = np.where(aus)[0]
        idx_activ = np.where(~aus)[0]
        AAT = np.vstack([AAT[np.sum(~np.isnan(AAT), axis=1) > 0], 
                         AAT[np.sum(~np.isnan(AAT), axis=1) == 0]])
        Unit_order = np.concatenate([idx_activ + 1, idx_nan + 1])
        As_order = np.arange(1, len(As_across_bins) + 1)

    elif display == 'clustered':
        aus = np.all(np.isnan(AAT), axis=1)
        idx_nan = np.where(aus)[0]
        idx_activ = np.where(~aus)[0]
        aus = ~np.all(np.isnan(AAT), axis=1)
        AAT_activ = AAT[aus]

        # Order based on units co-occurrence
        A01 = ~np.isnan(AAT_activ)
        M_assemb = np.zeros((AAT_activ.shape[0], AAT_activ.shape[0]))

        for n in range(A01.shape[1]):
            aus = np.where(A01[:, n])[0]
            M_assemb[np.ix_(aus, aus)] += 1

        M_assemb[M_assemb == 0] = 0.0001
        d_assemb = 1 / M_assemb
        np.fill_diagonal(d_assemb, 0)
        Q = linkage(d_assemb, 'average')
        _, _, perm = dendrogram(Q, no_plot=True)
        perm1 = idx_activ[perm]

        # Order based on assemblies distance
        D = pdist(~np.isnan(AAT_activ).T.astype(float))
        Z = squareform(D)
        Q2 = linkage(Z, 'average')
        _, _, perm2 = dendrogram(Q2, no_plot=True)

        AAT = AAT_activ[perm][:, perm2]
        AAT = np.vstack([AAT, np.full((len(idx_nan), AAT.shape[1]), np.nan)])
        Unit_order = np.concatenate([perm1 + 1, idx_nan + 1])
        Binvector = Binvector[perm2]
        As_order = perm2 + 1

    Amatrix = AAT
    binmat = np.log10(Binvector)
    AAT = np.pad(AAT, ((0, 1), (0, 1)), mode='constant', constant_values=np.nan)
    AAT[-1, -1] = 0
    binmat = np.pad(binmat, (0, 1), mode='constant', constant_values=0)

    # Plotting
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

    im1 = ax1.pcolormesh(AAT, cmap='jet', edgecolors='lightgray', linewidth=0.5)
    ax1.set_ylim(ax1.get_ylim()[::-1])  # Invert y-axis
    ax1.set_ylabel('Unit #')
    ax1.set_xticks([])
    cbar1 = fig.colorbar(im1, ax=ax1)
    cbar1.set_label('Time lag l (# bins)')
    ax1.set_yticks(np.arange(0.5, len(Unit_order) + 0.5))
    ax1.set_yticklabels(Unit_order)

    im2 = ax2.pcolormesh(binmat.reshape(1, -1), cmap='jet')
    ax2.set_xlabel('Assembly #')
    ax2.set_yticks([])
    cbar2 = fig.colorbar(im2, ax=ax2, orientation='horizontal')
    cbar2.set_label('Temporal precision Δ (sec)')
    
    if len(BinSizes) > 1:
        cbar2.set_clim(np.log10(min(BinSizes)), np.log10(max(BinSizes)))
    else:
        cbar2.set_clim(np.log10(BinSizes - 0.001), np.log10(BinSizes + 0.001))
    
    L = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03,
         0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
         1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    cbar2.set_ticks(np.log10(L))
    cbar2.set_ticklabels(L)

    ax2.set_xticks(np.arange(0.5, len(As_order) + 0.5, 2))
    ax2.set_xticklabels(As_order[::2])

    plt.tight_layout()
    plt.show()

    return Amatrix, Binvector, Unit_order, As_order