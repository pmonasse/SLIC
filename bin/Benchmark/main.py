import argparse
import BenchmarkUtility as bu


parser = argparse.ArgumentParser(
    description=("Pass the images under the selected root directory to SLIC" +
                 "and store time results.")
)

parser.add_argument('--root',
                    '-r',
                    dest='rootPath',
                    help=("The root directory from which the benchmark " +
                          "should start (every .png or .jpg file under " +
                          "this directory will be passed to SLIC.exe)."),
                    type=str)

parser.add_argument('--stop',
                    '-s',
                    dest='stop',
                    help="Optional: specify the number of images to test.",
                    type=int,
                    default=None)

parser.add_argument('--destfile',
                    '-d',
                    default="slicBenchmark.txt",
                    dest='destFile',
                    help=("The file where to write the results. " +
                          "/!\\ Existing files will be overwritten. " +
                          "Default: 'slicBenchmark.txt'"))

parser.add_argument('--verbose',
                    action='store_true',
                    dest='verbose',
                    help="Print everything.")

args = parser.parse_args()


benchmarker = bu.BenchmarkUtility(
    args.rootPath,
    args.stop,
    args.destFile,
    args.verbose
)


if __name__ == "__main__":
    benchmarker.runBenchmark()
