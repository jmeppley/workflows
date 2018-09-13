"""
Helper functions for the calculation of normalization factors from spiked
sequences
"""
import numpy
import pandas

def calculate_factors(counts_table, spiked_amounts, table_format=None):
    """ calculate a scaling factor from two input tables """
    # load counts of recovered standards from this sample
    pandas_args = {'header': None, 'index_col': None, 'skiprows':1,
                   'names':['Ref', 'Counts']}
    if table_format == 'bbduk':
        pandas_args.update({'skiprows':4,
                            'names':['Ref', 'Counts', 'Pct']})

    count_table = pandas.read_table(counts_table, **pandas_args)
    count_table['Ref'] = [r.split()[0] for r in count_table['Ref']]
    count_table.set_index('Ref', inplace=True)

    # load spiked in amounts
    spike_table = pandas.read_table(spiked_amounts, header=None, index_col=0,
                                   names=['Ref', 'Spiked'])

    # get data as lists in same order
    standard_list = sorted(list(spike_table.index))
    counts = [count_table.Counts.get(s,0) for s in standard_list]
    spikes = [spike_table.Spiked[s] for s in standard_list]

    # calculate the scale factor and save
    scale_factor = get_best_fit(counts, spikes, force_intercept=True)[0]
    return scale_factor

def get_best_fit(xd, yd, force_intercept=False, force_slope=False):
    """Return coeefs for a line of best fit"""
    #Calculate trendline
    if force_intercept:
        # intercept of 0
        x = numpy.array(xd)[:,numpy.newaxis]
        slope, _, _, _ = numpy.linalg.lstsq(x, yd)
        coeffs = [slope[0], 0]
        if force_slope:
            # We shouldn't get here, but let's just return the fixed values
            coeffs = (1, 0)
    elif force_slope:
        # slope of 1: intercept is average of difference
        intercept = numpy.mean(yd-xd)
        coeffs = [1,intercept]
    else:
        coeffs = numpy.polyfit(xd, yd, 1)

    return coeffs


