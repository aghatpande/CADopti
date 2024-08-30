import numpy as np

def restyle_assembly_lags_time(A, A_index):
    """
    REORDER LAGS AND SHIFT ASSEMBLY'S OCCURRENCE
    Lags and elements are reordered so that all lags are expressed with respect
    to the first firing unit. Assembly activation times are shifted
    accordingly to match the firing of the first assembly unit.

    © 2020 Russo
    For information please contact eleonora.russo@zi-mannheim.de
    """
    As_restyled = [None] * len(A)
    As_restyled_index = A_index.copy()

    num_ass = len(As_restyled)
    for i in range(num_ass):
        llag = np.concatenate(([0], A[i]['lag']))
        sorted_indices = np.argsort(llag)
        a = llag[sorted_indices]
        min_lag = np.min(a)

        As_restyled[i] = {
            'elements': [A[i]['elements'][j] for j in sorted_indices],
            'lag': a - min_lag,
            'Time': np.zeros_like(A[i]['Time']),
            'pr': A[i]['pr'][-1],
            'Noccurrences': A[i]['Noccurrences'][-1],
            'bin': A[i]['bin']
        }

        aus = A[i]['Time'].astype(float)
        act_bins = np.nonzero(aus)[0]
        As_restyled[i]['Time'][act_bins + min_lag] = A[i]['Time'][act_bins]

    return As_restyled, As_restyled_index