def assemblies_across_bins(assembly, bin_sizes):
    """
    Rearrange assemblies of all bins in a different structure (no assembly modification, only formal).

    Args:
    assembly (dict): Dictionary containing assembly information
    bin_sizes (list): List of investigated bin sizes

    Returns:
    tuple: (as_across_bins, as_across_bins_index)
        as_across_bins (list): List of dictionaries containing assembly information
        as_across_bins_index (list): Information to link assemblies in "as_across_bins"
                                     back to the structure "assembly"
    """
    as_across_bins = []
    as_across_bins_index = []

    # Find the first non-empty bin
    for e, bin_data in enumerate(assembly['bin']):
        if 'n' in bin_data:
            n_ass = len(bin_data['n'])
            as_across_bins = bin_data['n']
            for i in range(n_ass):
                as_across_bins[i]['bin'] = bin_sizes[e]
                as_across_bins_index.append([e, i])
            j = n_ass
            identity = list(range(n_ass))
            id_count = n_ass
            break
    else:
        raise ValueError("No non-empty bins found")

    # Process remaining bins
    for gg in range(e + 1, len(assembly['bin'])):
        bin_data = assembly['bin'][gg]
        if 'n' in bin_data:
            n_ass = len(bin_data['n'])
            for i in range(n_ass):
                new_assembly = {
                    'elements': bin_data['n'][i]['elements'],
                    'pr': bin_data['n'][i]['pr'],
                    'lag': bin_data['n'][i]['lag'],
                    'Time': bin_data['n'][i]['Time'],
                    'Noccurrences': bin_data['n'][i]['Noccurrences'],
                    'bin': bin_sizes[gg]
                }
                as_across_bins.append(new_assembly)
                as_across_bins_index.append([gg, i])
                id_count += 1
                identity.append(id_count)
                j += 1

    # Reorder lags and shift assembly's occurrence
    as_across_bins, as_across_bins_index = restyle_assembly_lags_time(as_across_bins, as_across_bins_index)

    return as_across_bins, as_across_bins_index

# Note: The restyle_assembly_lags_time function is not provided in the original code,
# so you'll need to implement it separately or provide its implementation.