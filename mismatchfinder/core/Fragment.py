from numpy import array, union1d, sort, empty, fromiter, fromstring, int32


class Fragment(object):
    """Contains all info available for a dna fragment"""

    def __init__(self, read1, read2=None):

        super(Fragment, self).__init__()

        # get all the relevant info from read1 we would ever want
        r1QuerySeq = array(list(read1.query_sequence))
        r1QueryQual = array(read1.query_qualities, dtype=int32)
        r1RefPos = array(read1.get_reference_positions(full_length=True))
        # build a dictionary mapping the refPos to the readPos
        r1IndDict = dict((k, i) for i, k in enumerate(r1RefPos))
        # remove the None entry
        if None in r1IndDict:
            del r1IndDict[None]

        r1AlStart = read1.query_alignment_start
        r1AlEnd = read1.query_alignment_end

        r1MapQual = read1.mapping_quality
        r1RefName = read1.reference_name

        try:
            r1AlPairs = array(read1.get_aligned_pairs(with_seq=True))
            # i am not smart enough to do this in a one liner, so the more readable version it is
            r1RefSeq = {}
            for (readPos, refPos, seq) in r1AlPairs:
                if not refPos is None:
                    r1RefSeq[int(refPos)] = seq
        except ValueError:
            r1AlPairs = array([])
            r1RefSeq = array([])

        if read2 is not None:
            # get all the info from read2 as well
            r2QuerySeq = array(list(read2.query_sequence))
            r2QueryQual = array(read2.query_qualities, dtype=int32)
            r2RefPos = array(read2.get_reference_positions(full_length=True))
            # build a dictionary mapping the refPos to the readPos
            r2IndDict = dict((k, i) for i, k in enumerate(r2RefPos))
            # remove the None entry
            if None in r2IndDict:
                del r2IndDict[None]

            r2AlStart = read2.query_alignment_start
            r2AlEnd = read2.query_alignment_end

            r2MapQual = read2.mapping_quality
            r2RefName = read2.reference_name

            try:
                r2AlPairs = array(read2.get_aligned_pairs(with_seq=True))
                r2RefSeq = {}
                for (readPos, refPos, seq) in r2AlPairs:
                    if not refPos is None:
                        r2RefSeq[int(refPos)] = seq
            except ValueError:
                r2AlPairs = array([])
                r2RefSeq = array([])

            # now we try to combine the results of the two reads
            # we throw all ref positions into one pot to genrate one hybrid read
            refPosJoined = sort(
                union1d(
                    array(list(r1RefSeq), dtype=int32),
                    array(list(r2RefSeq), dtype=int32),
                )
            )
            querySeqJoined = empty(len(refPosJoined), dtype=str)
            queryQualJoined = empty(len(refPosJoined), dtype=int32)
            refSeqJoined = empty(len(refPosJoined), dtype=str)

            # go through all of the overlaps and decide which of the reads is better
            for i in range(0, len(refPosJoined)):
                refPos = refPosJoined[i]

                if (
                    r1RefName == r2RefName
                    and refPos in r1IndDict
                    and refPos in r2IndDict
                ):
                    # both reads have this position, so we can build a consensus
                    r1IntPos = r1IndDict[refPos]
                    r2IntPos = r2IndDict[refPos]

                    # first we set the things that are equal for the reads
                    refSeqJoined[i] = r1RefSeq[refPos]

                    # if there is a concordance between the two reads, we just adjust the base qual
                    # and are done with it
                    if r1QuerySeq[r1IntPos] == r2QuerySeq[r2IntPos]:
                        querySeqJoined[i] = r1QuerySeq[r1IntPos]
                        queryQualJoined[i] = (
                            r1QueryQual[r1IntPos] + r2QueryQual[r2IntPos]
                        )
                    else:
                        # however, if they arent the same, we take the base from the higher qual
                        # read but also reduce the base quality at that point
                        if r1QueryQuals[r1IntPos] > r2QueryQuals[r2IntPos]:
                            querySeqJoined[i] = r1QuerySeq[r1IntPos]
                            queryQualJoined[i] = r1QueryQuals[r1IntPos] - int(
                                read2Quals[read2IntPos] / 2
                            )

                        elif read1Quals[read1IntPos] < read2Quals[read2IntPos]:
                            querySeqJoined[i] = read2Seq[read2IntPos]
                            queryQualJoined[i] = read2Quals[read2IntPos] - int(
                                read1Quals[read1IntPos] / 2
                            )

                elif refPos in r1IndDict:
                    # only read one has info, so we only take r1
                    r1IntPos = r1IndDict[refPos]
                    querySeqJoined[i] = r1QuerySeq[r1IntPos]
                    queryQualJoined[i] = r1QueryQual[r1IntPos]
                    refSeqJoined[i] = r1RefSeq[refPos]
                elif refPos in r2IndDict:
                    # only read two has info, so we only take r2
                    r2IntPos = r2IndDict[refPos]
                    querySeqJoined[i] = r2QuerySeq[r2IntPos]
                    queryQualJoined[i] = r2QueryQual[r2IntPos]
                    refSeqJoined[i] = r2RefSeq[refPos]

                else:
                    # this shouldnt happen
                    print(read1)
                    print(read2)
                    print(i)
                    print(refPos)
                    raise Exception("Found overlapping readpos with no read attached")

            # Now we just have to assign all fields for later usage
            self.refSeq = refSeqJoined
            self.querySeq = querySeqJoined
            self.queryQual = queryQualJoined
            self.refPos = refPosJoined

            self.mapping_quality = (r1MapQual + r2MapQual) / 2

            self.is_paired = True
        else:
            # here we set all the thing for the single ended case
            self.is_paired = False

            # reference seq aligned to by this fragment
            self.refSeq = empty(len(r1RefSeq), dtype=str)
            self.refPos = empty(len(r1RefSeq), dtype=int32)
            # populate the arrays
            for i, key in enumerate(sorted(r1RefSeq.keys())):
                self.refSeq[i] = r1RefSeq[key]
                self.refPos[i] = key

            # query sequence of this fragment (aligned part)
            self.querySeq = empty(len(r1RefSeq), dtype=str)
            self.queryQual = empty(len(r1RefSeq), dtype=int32)
            for i, key in enumerate(sorted(r1IndDict.keys())):
                ind = r1IndDict[key]
                self.querySeq[i] = r1QuerySeq[ind]
                self.queryQual[i] = r1QueryQual[ind]

            self.mapping_quality = r1MapQual

        # the alignment positions have the last and the first base aligned to already
        self.fragment_start = self.refPos[0]
        self.fragment_end = self.refPos[len(self.refPos) - 1]
        self.refName = r1RefName

    def getMismatches(self, minBQ):
        """get the mismatches contained in this fragment"""

        mutations = []
        # we only go from 1 to n-1 because we cannot build trinucleotides otherwise
        # because python doesnt have the normal from i =0, we need a while
        i = 1
        while i < len(self.refSeq) - 1:

            if self.refSeq[i] == self.querySeq[i]:
                i = i + 1
                continue
            elif self.queryQual[i] < minBQ:
                i = i + 1
                continue
            else:
                # now that we have the standards out of the way, we can continue

                # if there is an indel hidden, the refPosition will not be consecutive
                if (
                    self.refPos[i - 1] != self.refPos[i] - 1
                    or self.refPos[i] != self.refPos[i + 1] - 1
                ):
                    continue

                # if the position before has a mismatch, we can build a DBS instead of a SBS, but
                # if the position ahead also has a mismatch, we dont know what to do, because so
                # far no TBS have been reported. so we can skip 3 positions, in order to get out of
                # this area of mismatches
                type = 1
                if self.refSeq[i - 1] != self.querySeq[i - 1]:

                    # TBS or even more @TODO: add check how far we can jump
                    if self.refSeq[i + 1] != self.querySeq[i + 1]:
                        i = i + 3
                        continue

                    # this is now a DBS instead of a SBS
                    type = 2

                # get the reference context and the query context (slicing is non inclusive for the
                # end )
                queryContext = self.querySeq[i - 1 : i + 2]
                refContext = self.refSeq[i - 1 : i + 2]

                # that would leave us with only single or dinucleotide variants which are shifted to
                # the left e.g.
                # alt: TTG     alt: CTG
                # ref: ccG     ref: CcG
                # which we can easily convert to the classes that we want which are the C>T changes in
                # the dinucleotide context CC,CT,TC because those are specific for melanoma
                mut = (
                    self.refName,
                    self.refPos[i],
                    "".join(queryContext),
                    "".join(refContext),
                    type,
                )

                # finsih this site up
                mutations.append(mut)
                i = i + 1

        return mutations
