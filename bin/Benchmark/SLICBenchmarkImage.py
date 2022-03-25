class SLICBenchmarkImage:
    '''An object containing an image (its name / path) and its SLIC
    outputs as attributes Defines order relationships between images
    (based on their sizes)

    '''
    def __init__(self, path, imgType, size, nSteps, slicTime, ccTime):
        self._path = path
        self._type = imgType
        self._size = size
        self._nSteps = nSteps
        self._slicTime = slicTime
        self._ccTime = ccTime

    @property
    def size(self):
        return self._size

    @property
    def get_slicSteps(self):
        return self._nSteps

    @property
    def get_slicTime(self):
        return self._SLICTime

    @property
    def get_ccTime(self):
        return self.CCTime

    @property
    def path(self):
        return self._path

    def __gt__(self, other):
        if(self._size > other._size):
            return True
        else:
            return False

    def __lt__(self, other):
        if(self._size < other._size):
            return True
        else:
            return False

    def __ge__(self, other):
        if(self._size >= other._size):
            return True
        else:
            return False

    def __le__(self, other):
        if(self._size <= other._size):
            return True
        else:
            return False
