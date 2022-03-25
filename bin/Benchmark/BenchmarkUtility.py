import numpy as np
import os
from PIL import Image, UnidentifiedImageError
import subprocess


class BenchmarkUtility:
    '''A little utility that is able to perform the benchmark, i. e. run
    SLIC on lots of images and store the results in a text file.'''

    def __init__(self,
                 rootPath=(os.path.abspath('.').split(os.path.sep)[0] +
                           os.path.sep),
                 stop=None,
                 destination='slic-benchmark.txt',
                 verbose=False,
                 ident=np.random.randint(0, 9999999999)):
        self._id = ident
        self._rootPath = rootPath
        self._stop = stop
        self._count = 0
        self._destination = destination
        self._verbose = verbose
        self._currImg = None  # (PIL.Image object, necessary??)
        self._currImgSize = None
        self._currImgName = None
        self._currImgPath = None
        self._currImgDestName = None
        self._currHeader = None
        self._possibleImgTypes = ['png', 'jpg', 'jpeg']

    def _constructImageHeader(self):
        '''Constructs a little header with image name, adress, number, and
        returns it plus the full destination img path, ready to use for
        SLICing

        '''
        header = f'{self._count}: '
        header += self._currImgPath
        header += '\n'
        header += f'Size: {self._currImgSize}\n'
        nuName = ('.' + os.sep +
                  "SLICed_Images" +
                  os.sep + "SLICed_" +
                  self._currImgName.split('.')[0] +
                  ".png")
        self._currImgDestName = nuName
        self._currHeader = header

    def getImage(self, pathToImg):
        '''Get the image at pathToImg'''
        # Open image
        try:
            currImg = Image.open(pathToImg)
        except UnidentifiedImageError:
            return False
        except FileNotFoundError:
            return False
        # Set image size
        self._currImgSize = currImg.size[0] * currImg.size[1]
        # Set name and path
        self._currImgPath = pathToImg
        self._currImgName = pathToImg.split(os.sep)[-1]
        self._constructImageHeader()
        currImg.close()
        return True

    def _runSlicOnCurrImg(self):
        '''Runs SLIC on current image and returns the output as as string'''
        output = ''
        try:
            output = subprocess.check_output(
                ['../SLIC', '-k 300', '-m 30', '-g 0',
                 self._currImgPath, self._currImgDestName],
                universal_newlines=True
            )
        except subprocess.CalledProcessError:
            os.mkdir('./SLICed_Images')
            output = subprocess.check_output(
                ['../SLIC', '-k 300', '-m 30', '-g 0',
                 self._currImgPath, self._currImgDestName],
                universal_newlines=True
            )
        finally:
            return output

    def _checkRunNow(self):
        '''A little checker for the user, returns True if user types 'y',
        False otherwise'''
        userVerif = ''
        while userVerif != 'y':
            in_str = "This is going to take some time. Run now? (y/n) "
            userVerif = input(in_str)
            if userVerif == 'n':
                break
            elif userVerif != 'y':
                print("Please type 'y' or 'n'.")
        return userVerif == 'y'

    def _checkImgFile(self, fileName):
        try:
            return fileName.split('.')[1] in self._possibleImgTypes
        except IndexError:
            return False

    def runBenchmark(self):
        '''Run SLIC on every image under self._rootPath and stores the
        results.'''
        if not self._checkRunNow():
            raise KeyboardInterrupt("You chose not to run now.")

        # Open file (should I do a with ... construction?)
        txtFile = open(self._destination, 'w')
        # Print ID in it
        txtFile.write(f"Benchmark ID: {self._id}.\n")

        for root, directories, files in os.walk(self._rootPath):
            # Walk through all files
            for name in files:
                # If the number of images has been reached, break
                if self._stop is not None and self._count >= self._stop:
                    print("The specified number of images has been reached")
                    return

                # If this file is not an image, go to next file
                if not self._checkImgFile(name):
                    continue

                # Get the image with PIL.Image
                abspath = os.path.abspath(os.path.join(root, name))
                if not self.getImage(abspath):
                    continue
                # Count it
                self._count += 1
                # Construct header
                self._constructImageHeader()

                # If verbose has been asked, print everything
                if self._verbose:
                    print(self._currHeader)
                else:
                    print('SLICing image', self._count)

                # Run SLIC on current image
                out = self._runSlicOnCurrImg()

                # Again, print everything if verbose
                if self._verbose:
                    print(out)

                # Write results in text file
                txtFile.write(self._currHeader)
                txtFile.write(out)

        # Close file
        txtFile.close()
