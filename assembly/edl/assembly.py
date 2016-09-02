import os, sys, numpy, logging, re, pandas
from Bio import SeqIO, SeqUtils
logger=logging.getLogger(__name__)

###
# hook for scriptifying
def main():
    """
    Simple hook for running some of the functions below as a script. Only works with positional arguments that are strings.

    Examples:
     python assembly.py contig_read_counts file_1
       will run countig_read_counts("file_1")
       """

    log_level = logging.WARN
    while sys.argv[1]=='-v':
        sys.argv.pop(1)
        log_level-=10
    logger.setLevel(log_level)
    logger.debug("Log level: {}".format(log_level))

    function = eval(sys.argv[1])
    args=[]
    kwargs={}
    for arg in sys.argv[2:]:
        try:
            param,value = arg.split("=",1)
            try:
                value=eval(value)
            except NameError:
                pass
            kwargs[param]=value
        except ValueError:
            args.append(arg)

    logger.debug("Function: {}\nArguments: {}\nKWArgs: {}".format(function,
                                                                  args,
                                                                  kwargs))
    function(*args,**kwargs)

###
# Code for getting contig stats contigs file
#
def get_contig_stats(contigs_fasta,
                     contig_depth_file=None,
                     contig_read_counts_file=None,
                     contig_stats_file=None,
                     contig_histogram_file=None,
                     **kwargs):
    """
    Extracts GC and length from contigs fasta

    CAn optionally merge with read counts and mapped coverage if
    samtools output files given.

    provide a contig_stats_file location to write data to disk instead of just returning a pandas DataFrame.

    provide contig_histogram_file to produce a file with summary stats and histgrams for each metric. See contig_length_stats() and numpy.histogram() for additional kwargs that can be passed when using this option.
    """
    # parse contigs fasta
    logger.info("Parsing contig fasta file: {}".format(contigs_fasta))
    contig_stats = get_stats_from_contigs(contigs_fasta)
    
    # add other files if requested
    if contig_read_counts_file is not None:
        # read counts
        logger.info("Parsing read count file: {}"\
                            .format(contig_read_counts_file))
        read_count_table = pandas.read_table(contig_read_counts_file,delim_whitespace=True,names=['read count','contig']).set_index('contig')
        contig_stats=contig_stats.join(read_count_table,how='left')

    if contig_depth_file is not None:
        # convert base by base depth data into coverage
        logger.info("Parsing read depth file: {}"\
                            .format(contig_depth_file))
        mapping_depth_table = get_samtool_depth_table(contig_depth_file)
        contig_stats=contig_stats.join(mapping_depth_table, how='left')
    
    # sort and get cumulative length
    contig_stats.fillna(0,inplace=True)
    contig_stats.sort_values(by='length',ascending=False,inplace=True)
    contig_stats['cumul length']=contig_stats.length.cumsum()
    for col in ['length', 'read count','mx cov','mn cov','cumul length']:
        if col in contig_stats.columns:
            contig_stats[col]=contig_stats[col].astype(int)

    if contig_stats_file is not None:
        logger.info("Writing stats table to: {}".format(contig_stats_file))
        contig_stats.to_csv(contig_stats_file,sep='\t',float_format="%0.2f")

    if contig_histogram_file is not None:
        with open(contig_histogram_file,'w') as OUTF:
            if 'min_lengths' in kwargs:
                min_lengths = kwargs.pop('min_lengths')
            elif 'min_length' in kwargs:
                min_lengths=[kwargs.pop('min_length'),]
            else:
                min_lengths = [0,500,2000]

            for i,min_length in enumerate(min_lengths):
                logger.info("Making report for contigs >= {}"\
                                .format(min_length))
                if i>0:
                    OUTF.write("\n==============================="\
                                +"================\n\n")
                OUTF.write("CONTIGS longer or equal to {}bp:\n\n"\
                              .format(min_length))
                OUTF.write(contig_length_stats(contig_stats,
                                               return_type='report',
                                               min_length=min_length,
                                               **kwargs
                                              ))

    return contig_stats

def get_stats_from_contigs(contigs_fasta):
    """
    Use BioPython parser and GC calculator to get contig lengths and GCs from contigs fasta
    """
    
    # initialize lists
    contigs=[]
    lengths=[]
    gcs=[]
    
    # loop over fasta records (this is 2-3 times faster than SeqIO.parse)
    # (and only marginally slower than my custom built parser.)
    with open(contigs_fasta,'r') as CF:
        for title, sequence in SeqIO.FastaIO.SimpleFastaParser(CF):
            # parse title with RegEx
            contig = title.split(None,1)[0]
            length = len(sequence)
            contigs.append(contig)
            lengths.append(length)
            gcs.append(SeqUtils.GC(sequence))
        
    # convert to DataFrame and return
    return pandas.DataFrame({'contig':contigs,'length':lengths,'GC':gcs}).set_index('contig')

def get_samtool_depth_table(depth_file):
    """
    Calculate coverage stats for each contig in an assembly

    Params:
     depth_file: output file from the command:
                    `samtools depth reads.v.contigs.bam`
                 this is a 3 column file with one line per base.
                 columns are:
                     'contig_id base_index base_depth'

    Returns:
     pandas.DataFrame with one row per contig and the three following columns:
            contig  av cov  mx cov
            """
    with open(depth_file,'r') as DEPTHS:
        return get_samtool_depth_table_from_handle(DEPTHS)

def get_samtool_depth_table_from_handle(depth_stream):
    """
    Calculate coverage stats for each contig in an assembly

    Params:
     depth_stream: output file from the command:
                    `samtools depth reads.v.contigs.bam`

                    passed as an open file-like object (aka a file handle)
                 this is a 3 column file with one line per base.
                 columns are:
                     'contig_id base_index base_depth'

    Returns:
     pandas.DataFrame with one row per contig and the three following columns:
            contig  av cov  mx cov
            """

    # reading into lists is a fast way to build a big DataFrame
    contigs, av_covs, mn_covs, mx_covs = [], [], [], []

    # loop over contig bases
    current_contig=None
    for line in depth_stream:
        contig, base, depth = line.split()
        depth=int(depth)
        if contig!=current_contig:
            if current_contig is not None:
                # end of contig, save numbers
                contigs.append(current_contig)
                av_covs.append(depths/bases)
                mn_covs.append(min_depth)
                mx_covs.append(max_depth)
            bases=0
            depths=0
            max_depth=depth
            min_depth=depth
            current_contig = contig

        # update contig numbers with current base
        bases+=1
        depths+=depth
        min_depth=min(depth,min_depth)
        max_depth=max(depth,max_depth)

    # end of final contig, save numbers
    contigs.append(current_contig)
    av_covs.append(depths/bases)
    mn_covs.append(min_depth)
    mx_covs.append(max_depth)

    return pandas.DataFrame({'contig':contigs,'av cov':av_covs,'mx cov':mx_covs,'mn cov':mn_covs},columns=['contig','av cov','mn cov','mx cov']).set_index('contig')

## 
# this evolved from (but now bears little resemblance to) the
# assemlbly_quality_stats.py script by:
# Author: Travis Poulsen
# Date: 09 Feb. 2013
# http://travispoulsen.com/blog/2013/07/basic-assembly-statistics/
# https://gist.github.com/tpoulsen/422b1a19cbd8c0f514fe/raw/assembly_quality_stats.py
def contig_length_stats(contig_stats, return_type=None, 
                                      txt_width=0, 
                                      log=False, 
                                      min_length=0, 
                                      **kwargs):
    """
    Given contig stats table
     * calculate length stats (including N50)
     * optionally plot histogram (use txt_width and backend to select format)
       (if txt_width is greater than 0 (should be at least 40 for a good plot))
     * return_types:
       None: just print text to STDOUT
       'report': return text
       'data': return dictionary of data
    """
    report_data = {"min_length":min_length}
    contig_stats =  contig_stats.loc[contig_stats.length >= min_length]

    if contig_stats.shape[0]==0:
        report_data['Assembly']={'count':0}
    else:

        report_data['Assembly'] = get_N_stats(contig_stats)
        for column,label in {'length':'Contig Lengths', 
                             'read count':'Reads per Contig',
                             'av cov':'Mean Mapped Depth',
                             'mx cov':'Maximum Mapped Depth',
                             'mn cov':'Minimum Mapped Depth',
                             'GC':'GC Content'}.items():
            if column not in contig_stats.columns:
                continue
            report_data[label] = get_column_stats(contig_stats[column])
            if txt_width>0:
                report_data[label]['log'] = "(log)" if log else ""
                report_data[label]['histogram'] = \
                        asciiHistogram(numpy.histogram(contig_stats[column], 
                                                       **kwargs),
                                       log=log,
                                       width=txt_width)

    if return_type=='data':
        return report_data

    report = get_contig_stats_report(report_data)
    if return_type is None:
        print(report)
    else:
        return report

def get_contig_stats_report(report_data):
    """
    return a formatted string summarizing contig length data
    """
    N_stats = report_data.pop("Assembly")
    if N_stats['count']==0:
        return """\
Assembly Summary Stats:
    Contigs: 0
"""

    report = """\
Assembly Summary Stats:
    Contigs: {count}
    N50:     {N50}
    N75:     {N75}
    N90:     {N90}
""".format(**N_stats)

    for column in ['Contig Lengths', 
                   'Reads per Contig',
                   'GC Content',
                   'Mean Mapped Depth',
                   'Maximum Mapped Depth',
                   'Minimum Mapped Depth',
                   ]:
        if column not in report_data:
            continue
        report += """
Summary of {column}:
    Min:    {min}
    Max:    {max}
    Mean:   {mean}
    Median: {median}
    StdDev: {std}
""".format(column=column, **report_data[column])
        if 'histogram' in report_data[column]:
            report += """
Histogram of {column} {log}:
{histogram}
""".format(column=column, **report_data[column])
    return report

def get_N_stats(contig_stats, N_levels=[50,75,90]):
    """
    uses "length" and "cumulength" columns in contig_stats table
    to quickly get N50 and others for given N_levels (def: 50,75,90)
    """
    N_stats = {'count':contig_stats.shape[0]}

    # tota length is last cumulative length value
    total_length = contig_stats['cumul length'].iloc[-1]
    # set up iterator over just these columns
    cumulen_iter = iter(contig_stats[['length','cumul length']].iterrows())

    # Loop over N's. Since they are sorted, we don't need to restart
    #  the length/cumul_length iterator
    cumulength = 0
    for N in sorted(N_levels):
        # looking for N% of the total length
        target = total_length * N / 100
        
        # stop when we get there
        while cumulength < target:
            contig, (length, cumulength) = next(cumulen_iter)

        # Save the contig length that got use here
        N_key = "N{0:02d}".format(N)
        N_stats[N_key]=length

    return N_stats

def get_column_stats(data):
    """
    return a dict of useful stats
    """
    return {'min':data.min(),
            'max':data.max(),
            'mean':data.mean(),
            'median':data.median(),
            'std':data.std()
            }

def asciiHistogram(histogram, log=False, width=60, label='length', maxLabelWidth=10):
    (values,edges)=histogram[:2]
    
    maxValue=max(values)
    
    centers=[int(float(sum(edges[i:i+2]))/2.) for i in range(len(values))]
    largestLabel = max(max([len(str(c)) for c in centers]),len(label))
    if largestLabel<6:
        largestLabel=6
    elif largestLabel>maxLabelWidth:
        largestLabel=maxLabelWidth
    
    plotWidth=width-largestLabel+1
    
    midPoint = numpy.exp((numpy.log(maxValue)-numpy.log(.5))/2) if log else maxValue/2
    output="%s|count%s%s|%s%s|\n" % (rightPad(label,largestLabel),
                                     "".join([" " for i in range(int(plotWidth/2) - len(str(int(midPoint))) - len("count"))]),
                                     str(int(midPoint)),
                                     "".join([" " for i in range(int(numpy.ceil(plotWidth/2.)) - 1 - len(str(int(maxValue))))]),
                                     str(int(maxValue)),
                                     )
    #output+="%s|%s\n" % ("".join(["_" for i in range(largestLabel)]),
    #                     "".join(["_" for i in range(plotWidth)]),
    #                     )
    for i, v in enumerate(values):
        output+="%s|%s\n" % (rightPad(str(centers[i]),largestLabel),getBarString(v, maxValue, plotWidth, log))
    return output

logChars=['-','~','=','#']
def getBarString(value, maxValue, maxWidth, log):
    """
    return string of various signs (-,~,=,#) based on value and scale
    """
    if log:
        value=numpy.log(value)-numpy.log(.5) if value>0 else 0
        maxValue=numpy.log(maxValue)-numpy.log(.5)
    width=maxWidth*(value/float(maxValue))
    if width<1:
        return ''
    char=logChars[0]
    s=char
    while len(s)<width:
        if log:
            #print "s: %s, mw: %s, w: %s" % (s, maxWidth, width)
            char=logChars[int(numpy.ceil(len(logChars)*len(s)/float(maxWidth))-1)]
        s+=char
    return s

def rightPad(name, width):
    if width<6:
        width=6
        logger.warn("Can't force names to be fewer than 6 characters")
    if len(name)>width:
        # remove middle and insert elipsis to get to -width- characters
        return name[:width-4]+'***'+name[-1:]
    while len(name)<width:
        # pad with trailing space to get to 13 characters
        name+=' '
    return name

if __name__ == '__main__':
    main()
