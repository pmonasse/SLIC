import numpy as np
import os
from PIL import Image
import subprocess
import SLICBenchmark as sb


class BenchmarkUtility:
    '''A little utility that is able to perform the benchmark, i. e. run
    SLIC on lots of images and store the results in a text file.'''

    def __init__(self,
                 id=np.random.randint(0, 999999999999999),
                 rootPath=(os.path.abspath('.').split(os.path.sep)[0] +
                           os.path.sep),
                 stop=None,
                 destination='slic-benchmark.txt',
                 verbose=False):
        self._id = id
        self._rootPath = rootPath
        self._stop = stop
        self._destination = destination
        self._verbose = verbose
        self._currImg = None  # (PIL.Image object, necessary??)
        self._currImgSize = None
        self._currImgName = None
        self._currImgPath = None
        self._currImgDestName = None
        self._currHeader = None

    def _constructImageHeader(self, imgPath, name, i):
        '''Constructs a little header with image name, adress, number, and
        returns it plus the full destination img path, ready to use for
        SLICing

        '''
        header = f'{i}: '
        header += self._currImgPath
        header += '\n'
        header += f'Size: {self._currImgSize}\n'
        nuName = ".\\SLICed_Images\\SLICed_" + self._currImgName + ".png"
        self._currImgDestName = nuName
        self._currHeader = header

    def getImage(self, pathToImg):
        '''Get the image at pathToImg'''
        # Open image
        currImg = Image.open(pathToImg)
        # Set image size
        self._currImgSize = currImg.size[0] * currImg.size[1]
        # Set name and path
        self._currImgPath = pathToImg
        self._currImgName = pathToImg.split(os.sep)[-1]
        self._constructImageHeader()
        currImg.close()

    def runSlicOnCurrImg():
        '''Runs SLIC on current image and stores the output as as string'''
        output = subprocess.check_output(
            ['../SLIC', '-k 300', '-m 30', '-g 0',
             self._currImgPath, self._currImgDestName]
        )

    def benchmarkSlic()
