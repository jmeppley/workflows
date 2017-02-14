"""
MEthods particular to pandase

 scan_pandaseq_log(log_file, msg_stream, debug=False):
     look through pandaseq log to find unpaired reads and return ids as a set
"""
import re
from snakemake.logging import logger

overlap_RE = re.compile(r'^\S+\s+INFO\s+BESTOLP\s+(\S+?)(?:\:(?:\d{1,2})|(?:\:[NACTG]+))?\s+(\S+)')
error_RE = re.compile(r'^\S+\s+ERR\s+(\S+)\s+(\S+)(?:\:\d)?\s+')
debug_RE = re.compile(r'^\S+\s+DBG')
read_index_dir_RE = re.compile(r'(#(\d+|[ACTGN]+))?(/[12])?')

def scan_pandaseq_log(pandaseq_log_fie, msg_stream, keep_debug=False):
    """
    scan pandaseq log to get stats and to find unpaired reads

    TODO: check for fatal errors
    """

    unpaired = set()
    line_count = 0
    overlap_count = 0
    error_counts = {}
    with open(pandaseq_log_fie) as stream:
        for line in stream:
            line_count += 1

            # logging
            if msg_stream is not None:
                # only print line if...
                if keep_debug:
                    # debugging is on
                    msg_stream.write(line)
                else:
                    # or this is NOT a DBG line.
                    if debug_RE.match(line) is None:
                        msg_stream.write(line)
                    else:
                        # may as well move to next line since we know this is DBG
                        continue

            # is this an BESTOLP line?
            match = overlap_RE.match(line)
            if match:
                overlap_count += 1
                read = match.group(1)
                score = match.group(2)
                if score == '-1':
                    unpaired.add(read)
                continue

            # is this an ERR line?
            match = error_RE.match(line)
            if match:
                read = match.group(2)
                code = match.group(1)
                error_counts[code] = 1 + error_counts.get(code, 0)
                if keep_debug:
                    msg_stream.write("%s Failed with %s" % (read, code))
                continue

        if keep_debug:
            msg_stream.write("%d overlap lines in %d total lines" % \
                                                    (overlap_count, line_count))
        unpaired_count = len(unpaired)
        msg = """
#===============================
# PandaSeq Complete:
#  Processed: %d
#  Paired:    %d
#  Unpaired:  %d
#  Errors (LOWQ => under quality threshold. These are normal): %s
#===============================\n""" % \
                        (overlap_count,
                         overlap_count - unpaired_count - \
                            sum(error_counts.values()),
                         unpaired_count,
                         sum(error_counts.values()))
        msg_stream.write(msg)

        error_counts['reads'] = overlap_count
        error_counts['paired'] = overlap_count-unpaired_count
        error_counts['unpaired'] = unpaired_count

        logger.debug("%d errors, %d unpaired" % (sum(error_counts.values()),
                                                 unpaired_count))
        return (unpaired, error_counts)


