import numpy as np

def assembly_activity_function(As_across_bins, assembly, spM, BinSizes, lagChoice='duration', act_count='full'):
    """
    This function returns assembly activity.

    Args:
    As_across_bins (list): List containing assembly information output of ASSEMBLIES_ACROSS_BINS
    assembly (dict): Dictionary containing assembly information output of Main_assemblies_detection.m
    spM (np.array): Matrix with population spike trains; each row is the spike train (time stamps, not binned) relative to a unit.
    BinSizes (list): List of bin sizes to be tested
    lagChoice (str, optional): 'beginning' for only assembly activation, 'duration' for activity during whole assembly duration
    act_count (str, optional): 'full' for full assembly occurrences, 'partial' for partial activations, 'combined' for weighted partial activations

    Returns:
    list: Assembly activation profile across recording time
    """

    assembly_activity = [None] * len(As_across_bins)
    
    for n in range(len(As_across_bins)):
        tb = assembly['bin'][next(i for i, size in enumerate(BinSizes) if size == As_across_bins[n]['bin'])]['bin_edges']
        
        if act_count == 'full':
            activity = detect_assembly(spM, As_across_bins[n], tb)
        elif act_count == 'partial':
            activity = detect_assembly_additive(spM, As_across_bins[n], tb)
        elif act_count == 'combined':
            activity = detect_assembly_additive_weighted(spM, As_across_bins[n], tb)
        
        ta = tb[:-1] + As_across_bins[n]['bin'] / 2
        WOO = np.column_stack((ta, activity))
        
        if lagChoice == 'beginning':
            assembly_activity[n] = WOO
        elif lagChoice == 'duration':
            alag = As_across_bins[n]['lag'][-1]
            activity_lag = np.vstack([np.pad(activity[:-i], (i, 0), 'constant') for i in range(alag + 1)])
            VOO = np.column_stack((ta, np.sum(activity_lag, axis=0)))
            assembly_activity[n] = VOO
    
    return assembly_activity

def detect_assembly(spM, Assembly, tb):
    """Detect full assembly occurrences"""
    elements = Assembly['elements']
    alag = Assembly['lag']
    Sbin_assembly = np.array([np.histogram(spM[e], tb)[0] for e in elements])
    
    Sbin_shifted = np.full_like(Sbin_assembly, np.nan)
    for e1, lag in enumerate(alag):
        if lag == 0:
            Sbin_shifted[e1] = Sbin_assembly[e1]
        elif lag > 0:
            Sbin_shifted[e1, :-lag] = Sbin_assembly[e1, lag:]
        else:
            Sbin_shifted[e1, -lag:] = Sbin_assembly[e1, :lag]
    
    aus = np.sum(Sbin_shifted, axis=0)
    Sbin_shifted[:, np.isnan(aus)] = 0
    return np.min(Sbin_shifted, axis=0)

def detect_assembly_additive(spM, Assembly, tb):
    """Detect partial assembly activations"""
    elements = Assembly['elements']
    alag = Assembly['lag']
    Sbin_assembly = np.array([np.histogram(spM[e], tb)[0] for e in elements])
    
    maxrate = np.max(Sbin_assembly)
    Sbin_p = [(Sbin_assembly >= i).astype(int) for i in range(1, maxrate + 1)]
    
    activity = np.zeros((maxrate, Sbin_assembly.shape[1]))
    for i, Sbin in enumerate(Sbin_p):
        Sbin_shifted = np.full_like(Sbin, np.nan)
        for e1, lag in enumerate(alag):
            if lag == 0:
                Sbin_shifted[e1] = Sbin[e1]
            elif lag > 0:
                Sbin_shifted[e1, :-lag] = Sbin[e1, lag:]
            else:
                Sbin_shifted[e1, -lag:] = Sbin[e1, :lag]
        
        activity[i] = np.sum(Sbin_shifted, axis=0)
        activity[i, activity[i] == 1] = 0
    
    activationtot = np.sum(activity, axis=0)
    activationtot[np.isnan(activationtot)] = 0
    return activationtot / len(elements)

def detect_assembly_additive_weighted(spM, Assembly, tb):
    """Detect weighted partial assembly activations"""
    elements = Assembly['elements']
    alag = Assembly['lag']
    Sbin_assembly = np.array([np.histogram(spM[e], tb)[0] for e in elements])
    
    maxrate = np.max(Sbin_assembly)
    Sbin_p = [(Sbin_assembly >= i).astype(int) for i in range(1, maxrate + 1)]
    
    activity = np.zeros((maxrate, Sbin_assembly.shape[1]))
    for i, Sbin in enumerate(Sbin_p):
        Sbin_shifted = np.full_like(Sbin, np.nan)
        for e1, lag in enumerate(alag):
            if lag == 0:
                Sbin_shifted[e1] = Sbin[e1]
            elif lag > 0:
                Sbin_shifted[e1, :-lag] = Sbin[e1, lag:]
            else:
                Sbin_shifted[e1, -lag:] = Sbin[e1, :lag]
        
        activity[i] = np.sum(Sbin_shifted, axis=0)
        activity[i] = activity[i] ** 5
        activity[i, activity[i] == 1] = 0
    
    activationtot = np.sum(activity, axis=0)
    activationtot[np.isnan(activationtot)] = 0
    return activationtot / (len(elements) ** 5)