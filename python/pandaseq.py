overlap_RE=re.compile(r'^\S+\s+INFO\s+BESTOLP\s+(\S+?)(?:\:(?:\d{1,2})|(?:\:[NACTG]+))?\s+(\S+)')
error_RE=re.compile(r'^\S+\s+ERR\s+(\S+)\s+(\S+)(?:\:\d)?\s+')
debug_RE=re.compile(r'^\S+\s+DBG')
read_index_dir_RE=re.compile(r'(#(\d+|[ACTGN]+))?(/[12])?')

def scan_pandaseq_log(pandaseq_log_fie, msg_stream, keep_debug=False):
    """
    scan pandaseq log to get stats and to find unpaired reads

    TODO: check for fatal errors
    """

    unpaired = set()
    lineCount=0
    olpCount=0
    error_counts={}
    with open(pandaseq_log_fie) as stream:
        for line in stream:
            lineCount+=1

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
            m=overlap_RE.match(line)
            if m:
                olpCount+=1
                read=m.group(1)
                score=m.group(2)
                if score == '-1':
                    unpaired.add(read)
                continue

            # is this an ERR line?
            m=error_RE.match(line)
            if m:
                read=m.group(2)
                code=m.group(1)
                error_counts[code]=1+error_counts.get(code,0)
                if keep_debug:
                    msg_stream.write("%s Failed with %s" % (read,code))
                continue

        if keep_debug:
            msg_stream.write("%d overlap lines in %d total lines" % \
                                                    (olpCount, lineCount))
        unpCount=len(unpaired)
        msg="""
#===============================
# PandaSeq Complete:
#  Processed: %d
#  Paired:    %d
#  Unpaired:  %d
#  Errors (LOWQ => under quality threshold. These are normal):
    %s#===============================\n""" % \
                        (olpCount,
                         olpCount-unpCount-sum(error_counts.itervalues()),
                         unpCount,
                         error_counts(error_counts))
        msg_stream.write(msg)

        errCount['reads']=olpCount
        errCount['paired']=olpCount-unpCount
        errCount['unpaired']=unpCount

        logger.debug("%d errors, %d unpaired" % (sum(errCount.values()),
                                                 len(unpaired)))
        return (unpaired,errCount)


def fakePairs(singles,forward,reverse,paired,gap=20,fastq=True,trim=None,inputFormat='fastq',errFile=None,batchSize=10000):
    """
    Given:
        a list of read names,
        forward and revers fastq files,
        file to append to
        gap size
    Pull out forward and reverse reads and join with N's as specified by gap
    WRite combined reads to 'paired' (can be filename or handle)
    """
    trimmingCounts={}

    # make sure singles is a set and remove runs of colons
    prev_len=len(singles)
    singles = set(re.sub(r':+',r':',s) for s in singles)
    if prev_len!=len(singles):
        ## Uh oh, this hack broke something
        raise Exception("The list of unpaired reads from pandaseq has name collisions if multiple colons are ignored. Your reads are using a naming convetion that the author of this software (illuminPrep) didn't anticipate. My apologies")

    # Open input files. Use gzip if needed
    if len(reverse)>3 and reverse[-3:]=='.gz':
        rrecords=SeqIO.parse(gzip.open(reverse), format=inputFormat)
    else:
        rrecords=SeqIO.parse(reverse, format=inputFormat)
    if len(forward)>3 and forward[-3:]=='.gz':
        frecords=SeqIO.parse(gzip.open(forward), format=inputFormat)
    else:
        frecords=SeqIO.parse(forward, format=inputFormat)

    rrecid_set = set()
    rrecid_count=0
    joinCount=0
    fwdCount=0
    revCount=0
    nulCount=0
    totalJoined=0
    totalWritten=0
    if gap>=0:
        gapStr='N'*gap
        gapQ=[0]*gap
    fakedJoins=[]
    if fastq:
        outputFormat=inputFormat
    else:
        outputFormat='fasta'
    for frec in frecords:
        # get matching reverse record and make sure they are still synced
        try:
            rrec=rrecords.next()
        except StopIteration:
            logger.warn("Too few records in reverse fastq: %s" % reverse)
            sys.exit(1)
            if rrec.id[:-1] != frec.id[:-1]:
                logger.warn("Mismatched IDs in forward and reverse fastq (%s :: %s)\n%s\n%s" % (frec.id,rrec.id,forward,reverse))

        # the pandaseq log strips index and direction from the read name
        rrecid_count+=1
        rrecid=read_index_dir_RE.sub('',rrec.id)
        if rrecid in rrecid_set:
            raise Exception("Removing the index and read direction (...#0/2) from the end of read names has caused duplicate names to appear. Your reads are using a naming convetion that the author of this software (illuminPrep) didn't anticipate. My apologies. (%s became %s)" % (rrec.id,rrecid))
        else:
            rrecid_set.add(rrecid)

        # ignore all but singles
        try:
            singles.remove(rrecid)
        except KeyError:
            # not a sing;e, move on to next
            continue

        # optional trim
        if trim is not None:
            (frec,fwhy)=trimEnds(frec,trim)
            trimmingCounts[fwhy]=trimmingCounts.get(fwhy,0)+1
            (rrec,rwhy)=trimEnds(rrec,trim)
            trimmingCounts[rwhy]=trimmingCounts.get(rwhy,0)+1

        # join seqs
        newRecs=[]
        if frec is None:
            if rrec is None:
                logger.debug("Both ends trimmed to oblivion: (%s,%s)" % (fwhy,rwhy))
                nulCount+=1
                continue
            logger.debug("Forward seq trimmed to oblivion (%s)" % (fwhy))
            newRecs.append(revComp(rrec,qual=fastq,suffix=".rev"))
            revCount+=1
        elif rrec is None:
            logger.debug("Reverse seq trimmed to oblivion (%s)" % (rwhy))
            newRecs.append(frec)
            fwdCount+=1
        elif gap>=0:
            newSeq=frec.seq + Seq.Seq(gapStr,frec.seq.alphabet) + rrec.seq.reverse_complement()
            newRec=SeqRecord.SeqRecord(newSeq,id=frec.id,name=frec.name,description="Faked join")
            joinCount+=1
            # join quality
            if fastq:
                newRec.letter_annotations['phred_quality'] = \
                    frec.letter_annotations['phred_quality'] + \
                    gapQ + \
                    list(reversed(rrec.letter_annotations['phred_quality']))
            newRecs.append(newRec)
        else:
            # gap < 0 means don't join...add separately
            newRecs.append(frec)
            newRecs.append(revComp(rrec,qual=fastq,suffix=".rev"))

        fakedJoins.extend(newRecs)
        if len(fakedJoins)>=batchSize:
            numWritten =  SeqIO.write(fakedJoins, paired, format=outputFormat)
            if numWritten != len(fakedJoins):
                logger.warn("Only %d of %d faked joins written!" % (numWritten, len(fakedJoins)))
            totalJoined+=len(fakedJoins)
            totalWritten+=numWritten
            del fakedJoins[:]

    # Confirm that the reverse record iterator is also finished
    try:
        rrec=rrecords.next()
        # should not get here, iterator should be done
        logger.warn("Extra records in reverse fastq (%s):\n%s" % (rrec.id,reverse))
    except StopIteration:
        # this is what we expect
        pass

    numWritten =  SeqIO.write(fakedJoins, paired, format=outputFormat)
    if numWritten != len(fakedJoins):
        logger.warn("Only %d of %d faked joins written!" % (numWritten, len(fakedJoins)))
    totalWritten+=numWritten
    totalJoined+=len(fakedJoins)

    # Report some counts
    msg="#======================\n# Faked joins\n#  Total: %s\n" % (totalJoined)
    if trim is not None:
        msg+="# Joined: %d\n# FwdOnly: %d\n# RevOnly: %d\n# Dropped: %d\n" % (joinCount, fwdCount, revCount, nulCount)
        msg+="#======================\n# End trimming:\n"
        for status,count in trimmingCounts.iteritems():
            msg+="# %s: %d\n" % (status,count)
    msg+="#======================\n"

    if errFile is not None:
        if isinstance(errFile,str):
            with open(errFile,'a') as errstream:
                errstream.write(msg)
        else:
            errFile.write(msg)

    if trim is not None:
        logger.info(msg)
        logger.debug("# Total written: %d" % (totalWritten))
    else:
        logger.info("# Faked joins: %d" % (totalJoined))

    if len(singles)>0:
        raise Exception("Some (%d) unpaired read names were not found! EG: %s" % (len(singles),singles.pop()))

    if len(rrecid_set)!=rrecid_count:
        raise Exception("%d != %d! Removing the index and read direction (...#0/2) from the end of read names has caused duplicate names to appear. Your reads are using a naming convetion that the author of this software (illuminPrep) didn't anticipate. My apologies" % (len(rrecid_set),rrecid_count))

    return numWritten


