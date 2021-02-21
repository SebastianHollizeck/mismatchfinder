from numpy import array, union1d, sort, empty, fromiter, fromstring, int32


class Fragment(object):
    """Contains all info available for a dna fragment"""

    def __init__(self, read1, read2=None):

        super(Fragment, self).__init__()

        # get all the relevant info from read1 we would ever want
        r1QuerySeq = array(list(read1.query_sequence), dtype=str)

        present = True
        try:
            read1.get_tag("MD")
        except ValueError:
            present = False

        r1RefSeq = {}
        if present:
            r1AlPairs = read1.get_aligned_pairs(with_seq=True, matches_only=True)
            for (readPos, refPos, seq) in r1AlPairs:
                r1RefSeq[int(refPos)] = (seq, readPos)

        if read2 is not None:
            # get all the info from read2 as well
            r2QuerySeq = array(list(read2.query_sequence), dtype=str)

            present = True
            try:
                read2.get_tag("MD")
            except ValueError:
                present = False

            r2RefSeq = {}
            if present:
                r2AlPairs = read2.get_aligned_pairs(with_seq=True, matches_only=True)
                for (readPos, refPos, seq) in r2AlPairs:
                    r2RefSeq[int(refPos)] = (seq, readPos)

            # now we try to combine the results of the two reads
            # we throw all ref positions into one pot to genrate one hybrid read
            refPosJoined = sort(
                union1d(
                    array(list(r1RefSeq), dtype=int32),
                    array(list(r2RefSeq), dtype=int32),
                )
            )
            querySeqJoined = empty(len(refPosJoined), dtype=str)
            queryQualsJoined = empty(len(refPosJoined), dtype=int32)
            refSeqJoined = empty(len(refPosJoined), dtype=str)

            # go through all of the overlaps and decide which of the reads is better
            for i in range(0, len(refPosJoined)):
                refPos = refPosJoined[i]

                if (
                    read1.reference_name == read2.reference_name
                    and refPos in r1RefSeq
                    and refPos in r2RefSeq
                ):
                    # both reads have this position, so we can build a consensus
                    refBase, r1IntPos = r1RefSeq[refPos]
                    refBase, r2IntPos = r2RefSeq[refPos]

                    # first we set the things that are equal for the reads
                    refSeqJoined[i] = refBase

                    # if there is a concordance between the two reads, we just adjust the base qual
                    # and are done with it
                    if r1QuerySeq[r1IntPos] == r2QuerySeq[r2IntPos]:
                        querySeqJoined[i] = r1QuerySeq[r1IntPos]
                        queryQualsJoined[i] = (
                            read1.query_qualities[r1IntPos]
                            + read2.query_qualities[r2IntPos]
                        )
                    else:
                        # however, if they arent the same, we take the base from the higher qual
                        # read but also reduce the base quality at that point
                        if (
                            read1.query_qualities[r1IntPos]
                            > read2.query_qualities[r2IntPos]
                        ):
                            querySeqJoined[i] = r1QuerySeq[r1IntPos]
                            queryQualsJoined[i] = read1.query_qualities[r1IntPos] - int(
                                read2.query_qualities[r2IntPos] / 2
                            )

                        elif (
                            read1.query_qualities[r1IntPos]
                            < read2.query_qualities[r2IntPos]
                        ):
                            querySeqJoined[i] = r1QuerySeq[r2IntPos]
                            queryQualsJoined[i] = read2.query_qualities[r2IntPos] - int(
                                read1.query_qualities[r1IntPos] / 2
                            )
                        else:
                            # in this case we just dont look at this site make it reference and
                            # reduce the quality to 0
                            querySeqJoined[i] = r1RefSeq[refPos]
                            queryQualsJoined[i] = 0

                elif refPos in r1RefSeq:
                    # only read one has info, so we only take r1
                    refBase, r1IntPos = r1RefSeq[refPos]
                    querySeqJoined[i] = r1QuerySeq[r1IntPos]
                    queryQualsJoined[i] = read1.query_qualities[r1IntPos]
                    refSeqJoined[i] = refBase
                elif refPos in r2RefSeq:
                    # only read two has info, so we only take r2
                    refBase, r2IntPos = r2RefSeq[refPos]
                    querySeqJoined[i] = r2QuerySeq[r2IntPos]
                    queryQualsJoined[i] = read2.query_qualities[r2IntPos]
                    refSeqJoined[i] = refBase

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
            self.queryQuals = queryQualsJoined
            self.refPos = refPosJoined

            self.mapping_quality = (read1.mapping_quality + read2.mapping_quality) / 2

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
            self.queryQuals = empty(len(r1RefSeq), dtype=int32)
            for i, key in enumerate(sorted(r1RefSeq.keys())):
                refBase, ind = r1IndDict[key]
                self.querySeq[i] = r1QuerySeq[ind]
                self.queryQuals[i] = read1.query_qualities[ind]

            self.mapping_quality = read1.mapping_quality

        # the alignment positions have the last and the first base aligned to already
        self.fragment_start = self.refPos[0]
        self.fragment_end = self.refPos[len(self.refPos) - 1]
        self.refName = read1.reference_name

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
            elif self.queryQuals[i] < minBQ:
                i = i + 1
                continue
            else:
                # now that we have the standards out of the way, we can continue

                # if there is an indel hidden, the refPosition will not be consecutive
                if (
                    self.refPos[i - 1] != self.refPos[i] - 1
                    or self.refPos[i] != self.refPos[i + 1] - 1
                ):
                    i = i + 1
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
