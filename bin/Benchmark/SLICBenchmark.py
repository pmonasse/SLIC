import SLICBenchmarkImage as sbi
import matplotlib.pyplot as plt


class SLICBenchmark:
    '''A set of SLICBenchmarkImages with some useful methods for making graphs'''
    def __init__(self, sourceFileName):
        self._sourceFile = sourceFileName
        self._data = self.readSourceFile(sourceFileName)
        self._data.sort()
        self._imgSizes = [piece.size() for piece in self._data]
        self._slicTimes = [piece.get_slicTime() for piece in self._data]
        self._ccTimes = [piece.get_ccTime() for piece in self._data]
        self._steps = [piece.get_slicSteps() for piece in self._data]

    def my_plot(X, Y, col, labl, lStyle, mark, markSize):
        plt.plot(X, Y, color=col, label=labl, ls=lStyle, marker=mark, markersize=markSize)

    def plot_slicTime(self, col='b', labl="SLIC Time", lStyle='', mark='+', markSize=4, silent=False):
        fig = plt.figure()
        ax = plt.gca()
        ax.plot(self._imgSizes, self._slicTimes,
                 color=col, label=labl, ls=lStyle, marker=mark, markersize=markSize)
        ax.legend()
        ax.set_xlabel("Image size (pixels)")
        ax.set_ylabel("Time (seconds)")
        ax.set_title("Time elapsed for SLIC in function of image size")
        plt.savefig("./slicTimes_fct_size.png")
        if not silent:
            plt.show()
        return ax, fig
    
    def plot_ccTime(self, col='r', labl="Connectivity enforcement time", lStyle='', mark='+', markSize=4, silent=False):
        fig = plt.figure()
        ax = plt.gca()
        ax.plot(self._imgSizes,
                 self._ccTimes,
                 color=col, label=labl, ls=lStyle, marker=mark, markersize=markSize)
        ax.legend()
        ax.set_xlabel("Image size (pixels)")
        ax.set_ylabel("Time (seconds)")
        ax.set_title("Time elapsed for connected components algorithm in function of image size")
        plt.savefig("ccTimes_fct_size.png")
        if not silent:
            plt.show()
        return ax, fig

    def plot_iterations(self, col='grey', labl="Number of iterations", lStyle='', mark='+', markSize=4, silent=False):
        fig = plt.figure()
        ax = plt.gca()
        ax.plot(self._imgSizes,
                 self._steps,
                 color=col, label=labl, ls=lStyle, marker=mark, markersize=markSize)
        ax.legend()
        ax.set_xlabel("Image size (pixels)")
        ax.set_ylabel("Numner of iterations")
        ax.set_title("Number of SLIC iterations in function of image size")
        plt.savefig("steps_fct_size.png")
        if not silent:
            plt.show()
        return ax, fig

    def plot_slic_and_cc(self, 
                         slicCol='b', slicLabl="SLIC Time", slicLStyle='', slicMark='+', slicMarkSize=4, 
                         ccCol='r', ccLabl="Connectivity enforcement time", ccLStyle='', ccMark='x', ccMarkSize=4,
                         silent=False):
        fig = plt.figure()
        plt.plot(self._imgSizes, self._slicTimes, col=slicCol, labl=slicLabl, lStyle=slicLStyle, mark=slicMark, markSize=slicMarkSize, silent=True)
        plt.plot(self._imgSizes, self._ccTimes, col=ccCol, labl=ccLabl, lStyle=ccLStyle, mark=ccMark, markSize=ccMarkSize, silent=True)
        plt.savefig("slic_and_cc_fct_time.png")
        if not silent:
            plt.show()
