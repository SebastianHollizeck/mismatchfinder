from datetime import datetime
from sys import stderr


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
