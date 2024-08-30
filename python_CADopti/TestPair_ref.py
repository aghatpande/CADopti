import numpy as np
from scipy import stats

def test_pair_ref(ensemble, spike_train2, n2, max_lag, dc, reference_lag):
    """
    Test if two spike trains have repetitive patterns occurring more frequently than chance.

    Parameters:
    - ensemble: structure with the previously formed assembly and its spike train
    - spike_train2: spike train of the new unit to be tested for significance
    - n2: new unit tested
    - max_lag: maximum lag to be tested
    - dc: size (in bins) of chunks to divide spike trains for variance computation
    - reference_lag: lag of reference; if zero or negative, reference_lag = -l

    Returns:
    - assem_d: dictionary containing assembly information
    """
    couple = np.column_stack([
        ensemble['Time'] - np.min(ensemble['Time']),
        spike_train2 - np.min(spike_train2)
    ])
    nu = 2
    ntp = couple.shape[0]

    # Divide into parallel trials with 0/1 elements
    max_rate = int(np.max(couple))
    zaa = [np.zeros(couple.shape, dtype=np.uint8) for _ in range(max_rate)]
    exp_abi = np.zeros(max_rate)
    
    for i in range(max_rate):
        zaa[i][couple >= i + 1] = 1
        exp_abi[i] = np.prod(np.sum(zaa[i], axis=0)) / couple.shape[0]

    # ... (skipping the lag calculation part for brevity)

    # Implement the lag calculation logic here
    # This part needs careful translation and testing

    # ... (skipping to the main calculation part)

    exp_ab = np.sum(exp_abi)

    if a == 0 or exp_ab <= 5 or exp_ab >= (min(np.sum(couple)) - 5):
        assem_d = {
            'elements': ensemble['elements'] + [n2],
            'lag': ensemble['lag'] + [99],
            'pr': ensemble['pr'] + [1],
            'Time': [],
            'Noccurrences': ensemble['Noccurrences'] + [0]
        }
    else:
        # Construct the activation time series for the couple
        len_couple = couple.shape[0]
        time = np.zeros(len_couple, dtype=np.uint8)

        # Implement the time series construction logic here
        # This part needs careful translation and testing

        # ... (skipping the chunking and variance calculation for brevity)

        # Implement the chunking and variance calculation logic here
        # This part needs careful translation and testing

        # ... (skipping to the final calculations)

        x = tp_r_m_tot - tp_r_m_tot.T
        if np.abs(x[0, 1]) > 0:
            x = np.abs(tp_r_m_tot - tp_r_m_tot.T) - 0.5  # Yates correction

        if var_x_tot[0, 1] == 0:
            pr_f = 1
        else:
            f = x ** 2 / var_x_tot
            pr_f = stats.f.sf(f[0, 1], 1, n)

        assem_d = {
            'elements': ensemble['elements'] + [n2],
            'lag': ensemble['lag'] + [l_],
            'pr': ensemble['pr'] + [pr_f],
            'Time': time,
            'Noccurrences': ensemble['Noccurrences'] + [np.sum(time)]
        }

    return assem_d