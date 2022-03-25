## Python benchmark tool

Usage: put the folder 'Benchmark' in the directory where your SLIC executable is located
and run
`$ cd /path/to/SLIC/dir/Benchmark
$ python3 main.py [args]`

The args are detailed here:

`
main.py [-h] [--root ROOTPATH] [--stop STOP] [--destfile DESTFILE]
               [--verbose]

Pass the images under the selected root directory to SLICand store time
results.

optional arguments:
  -h, --help            show this help message and exit
  --root ROOTPATH, -r ROOTPATH
                        The root directory from which the benchmark should
                        start (every .png or .jpg file under this directory
                        will be passed to SLIC.exe).
  --stop STOP, -s STOP  Optional: specify the number of images to test.
  --destfile DESTFILE, -d DESTFILE
                        The file where to write the results. /!\ Existing
                        files will be overwritten. Default:
                        'slicBenchmark.txt'
  --verbose             Print everything.`

You can get this message with
`$ python3 main.py -h`
or
`$ python3 main.py --help`
