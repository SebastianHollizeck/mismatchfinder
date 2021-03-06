from datetime import datetime
from logging import error
from sys import stderr

from matplotlib.pyplot import figure, show
from numpy import abs, max, min
from pandas import DataFrame, concat

from mismatchfinder.utils.Misc import DBSorder, SBSorder


class JSONWriter:
    """Writes the results to JSON format."""

    def __init__(self, file):
        super(JsonWriter, self).__init__()
        self.file = file

    def writeResults(results):
        pass


# function to attach the current time to the print statement
def printLog(str, end="\n", addTime=True):
    if addTime:
        str = datetime.now().strftime("%d/%m/%Y %H:%M:%S - ") + str
    print(str, end=end, file=stderr)


def plotStats(tumourDf, normalDf=None):

    # create the drawing plane
    fig = figure(figsize=(8, 8))

    # Add a gridspec with two rows and two columns and a ratio of 2 to 7 between
    # the size of the marginal axes and the main axes in both directions.
    # Also adjust the subplot parameters for a square plot.
    gs = fig.add_gridspec(
        nrows=1,
        ncols=2,
        width_ratios=(7, 2),
        left=0.1,
        right=0.9,
        bottom=0.1,
        top=0.9,
        wspace=0.05,
        hspace=0.05,
    )

    ax = fig.add_subplot(gs[0, 0])
    ax_histy = fig.add_subplot(gs[0, 1], sharey=ax)

    nRowsTumour = len(tumourDf.index)
    tumourX = nRowsTumour * [0]
    tumourY = tumourDf["nMisMatches"] / tumourDf["nAlignedReads"]

    minY = min(tumourY)
    maxY = max(tumourY)

    if not normalDf is None:
        nRowsNormal = len(normalDf.index)
        normalX = nRowsNormal * [0]
        normalY = normalDf["nMisMatches"] / normalDf["nAlignedReads"]
        minYNormal = min(normalY)
        if minYNormal < minY:
            minY = minYNormal

        maxYNormal = max(normalY)
        if maxYNormal > maxY:
            maxY = maxYNormal

    # paint the first picture and store the calculated bin width
    scatter_hist(tumourX, tumourY, ax, ax_histy, label="tumour", range=(minY, maxY))

    # we plot this seperatly so we can get a min and max of all data first
    if not normalDf is None:
        scatter_hist(normalX, normalY, ax, ax_histy, label="normal", range=(minY, maxY))

    ax.legend()
    show()


def scatter_hist(x, y, ax, ax_histy, range=None, label=None, bins=20):
    # no labels
    ax_histy.tick_params(axis="y", labelleft=False)
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)

    # the scatter plot:
    ax.scatter(x, y, label=label)

    # make a pretty histogram
    # first we adjust the range to not start and stop at min and max but a bit lower and higher
    (minY, maxY) = range
    minY = minY - minY / 20
    maxY = maxY + maxY / 20
    range = (minY, maxY)
    # and then we just add the plot in a horizontal orientation next to the scatter
    ax_histy.hist(y, range=range, bins=bins, orientation="horizontal", alpha=0.9)


def createOutputFiles(outFileRoot):

    # create the full path of the files
    SBSFile = outFileRoot.parent / (outFileRoot.name + "_SBScontexts.tsv")
    DBSFile = outFileRoot.parent / (outFileRoot.name + "_DBScontexts.tsv")
    statsFile = outFileRoot.parent / (outFileRoot.name + "_stats.tsv")

    # write the header to the files with one extra column where the name will go
    try:
        with SBSFile.open("w") as fh:
            fh.write("bam\t" + "\t".join(SBSorder) + "\n")

        with DBSFile.open("w") as fh:
            fh.write("bam\t" + "\t".join(DBSorder) + "\n")

        with statsFile.open("w") as fh:
            # we just want to create the file here
            pass
    except:
        error(f"Could not write into folder{outFileRoot.parent}")
        exit(1)


def writeStatsFile(pandas, outFileRoot):

    statsFile = outFileRoot.parent / (outFileRoot.name + "_stats.tsv")
    with statsFile.open("w") as fh:
        pandas.to_csv(fh, sep="\t", index=False)
