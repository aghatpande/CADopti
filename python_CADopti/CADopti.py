import numpy as np
from scipy import stats
import multiprocessing as mp

def CADopti(spM, MaxLags, BinSizes, ref_lag=None, alph=None, No_th=None, O_th=None, bytelimit=None):
    """
    This function returns cell assemblies detected in spM spike matrix binned 
    at a temporal resolution specified in 'BinSizes' vector and testing for all 
    lags between '-MaxLags(i)' and 'MaxLags(i)'

    USAGE: [assembly] = Main_assemblies_detection(spM, MaxLags, BinSizes, ref_lag, alph, Dc, No_th, O_th, bytelimit)

    ARGUMENTS:
    spM     := matrix with population spike trains; each row is the spike train (time stamps, not binned) relative to a unit. 
    BinSizes:= vector of bin sizes to be tested;
    MaxLags := vector of maximal lags to be tested. For a binning dimension of BinSizes(i) the program will test all pairs configurations with a time shift between -MaxLags(i) and MaxLags(i);
    (optional) ref_lag      := reference lag. Default value 2
    (optional) alph      := alpha level. Default value 0.05
    (optional) No_th     := minimal number of occurrences required for an assembly (all assemblies, even if significant, with fewer occurrences than No_th are discarded). Default value 0.
    (optional) O_th      := maximal assembly order (the algorithm will return assemblies of composed by maximum O_th elements).
    (optional) bytelimit := maximal size (in bytes) allocated for all assembly structures detected with a bin dimension. When the size limit is reached the algorithm stops adding new units. 

    RETURNS:
    assembly - structure containing assembly information:
        assembly['parameters']       - parameters used to run Main_assemblies_detection
        assembly['bin'][i] contains information about assemblies detected with 
                        'BinSizes[i]' bin size tested for all lags between 
                        '-MaxLags[i]' and 'MaxLags[i]'
    
           assembly['bin'][i]['bin_edges'] - bin edges (common to all assemblies in assembly['bin'][i])
           assembly['bin'][i]['n'][j] information about the j-th assembly detected with BinSizes[i] bin size 
                      elements: vector of units taking part to the assembly (unit order correspond to the agglomeration order)
                           lag: vector of time lags. '.lag[z]' is the activation delay between .elements[0] and .elements[z+1]
                            pr: vector of pvalues. '.pr[z]' is the pvalue of the statistical test between performed adding .elements[z+1] to the structure .elements[0:z]
                          Time: assembly activation time. It reports how many times the complete assembly activates in that bin. .Time always refers to the activation of the first listed assembly element (.elements[0]), that doesn't necessarily corresponds to the first unit firing.
                  Noccurrences: number of assembly occurrence. '.Noccurrences[z]' is the occurrence number of the structure composed by the units .elements[0:z+1] 

    As_across_bins         - structure containing assembly information (exactly same information contained in "assembly" but collected across different temporal resolutions)
    As_across_bins_index   - information to link assemblies in "As_across_bins"
                             back to the structure "assembly":
                             assembly As_across_bins[i] is contained in assembly['bin'][As_across_bins_index[i][0]]['n'][As_across_bins_index[i][1]].

    © 2020 Russo
    for information please contact eleonora.russo@zi-mannheim.de
    """

    if ref_lag is None:
        ref_lag = 2
    if alph is None:
        alph = 0.05
    if No_th is None:
        No_th = 0  # no limitation on the number of assembly occurrences
    if O_th is None:
        O_th = float('inf')  # no limitation on the assembly order (=number of elements in the assembly)
    if bytelimit is None:
        bytelimit = float('inf')  # no limitation on assembly dimension

    nneu = spM.shape[0]  # number of units
    testit = np.ones(len(BinSizes))
    binM = [None] * len(BinSizes)
    number_tests = 0

    # matrix binning at all bins
    for gg in range(len(BinSizes)):
        int_val = BinSizes[gg]
        tb = np.arange(np.min(spM), np.max(spM) + int_val, int_val)
        
        binM[gg] = np.zeros((nneu, len(tb) - 1), dtype=np.uint8)
        number_tests += nneu * (nneu - 1) * (2 * MaxLags[gg] + 1) // 2
        
        for n in range(nneu):
            binM[gg][n, :], _ = np.histogram(spM[n, :], tb)
        
        assembly = {'bin': [{'n': [], 'bin_edges': tb} for _ in range(len(BinSizes))]}
        
        if binM[gg].shape[1] - MaxLags[gg] < 100:
            print(f'Warning: testing bin size={int_val}. The time series is too short, consider taking a longer portion of spike train or diminish the bin size to be tested')
            testit[gg] = 0

    print('order 1')
    Assemblies_all_orders = []
    O = 1
    Dc = 100  # length (in # bins) of the segments in which the spike train is divided to compute #abba variance (parameter k).

    assembly_selected_xy = []
    p_values = []

    import multiprocessing as mp
from itertools import combinations

def process_pair(args):
    w1, w2, binM, MaxLags, BinSizes, Dc, ref_lag = args
    assemblybin = [None] * len(BinSizes)
    p_by_bin = []
    for gg in range(len(BinSizes)):
        assemblybin[gg] = FindAssemblies_recursive_prepruned(
            np.vstack((binM[gg][w1, :], binM[gg][w2, :])),
            w1, w2, MaxLags[gg], Dc, ref_lag
        )
        p_by_bin.append(assemblybin[gg]['pr'][-1])
        assemblybin[gg]['bin'] = BinSizes[gg]
    
    b = np.argmin(p_by_bin)
    return assemblybin[b], p_by_bin[b]

# First order assembly
print('order 1')
assembly_selected_xy = []
p_values = []

# Prepare arguments for parallel processing
pair_args = [
    (w1, w2, binM, MaxLags, BinSizes, Dc, ref_lag)
    for w1, w2 in combinations(range(nneu), 2)
]

# Use multiprocessing to parallelize the computation
with mp.Pool() as pool:
    results = pool.map(process_pair, pair_args)

# Process the results
for result, p_value in results:
    assembly_selected_xy.append(result)
    p_values.append(p_value)

assembly_selected = assembly_selected_xy
    # ... (rest of the function implementation)

    return As_across_bins, As_across_bins_index, assembly, Assemblies_all_orders

# Helper functions (to be implemented)
def FindAssemblies_recursive_prepruned(binM, w1, w2, MaxLag, Dc, ref_lag):
    # Implementation needed
    pass

def TestPair_ref(assembly, spikeTrain2, w2, MaxLag, Dc, ref_lag):
    # Implementation needed
    pass

def assemblies_across_bins(assembly, BinSizes):
    # Implementation needed
    pass