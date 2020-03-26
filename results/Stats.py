class Results(object):
    """This class contains all results accumulated after analysing a bam."""

    sample = None

    nMisMatches = None
    nSites = None
    nReads = None
    nAlignedReads = None
    nAlignedBases = None

    nSomaticMisMatches = None
    nSomaticSites = None

    nConfidentSites = None

    def __init__(self, sample):
        super(Results, self).__init__()
        self.sample = sample
