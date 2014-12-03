# Build and run supermodel tests. Optional command line arguments are the test module
# names required in the test, otherwise all test modules are used.

from sys import argv
from glob import glob
import FRUIT

build_command = 'make tests'
driver_name = 'test_all'
test_dir = 'test/src/'
output_dir = 'test/dist/'
testsuffix = '_test'
f90 = '.F90'
setup_source = test_dir + 'setup' + testsuffix + f90
num_procs = 4
mpi_comm = 'PETSC_COMM_WORLD'

if len(argv) > 1: test_sources = [test_dir + arg + testsuffix + f90 for arg in argv[1:]]
else: test_sources = glob(test_dir + '*' + f90)

driver_source = test_dir + driver_name + f90
if driver_source in test_sources: test_sources.remove(driver_source)
if setup_source not in test_sources: test_sources.append(setup_source)

suite = FRUIT.test_suite(test_sources)
suite.build_run(driver_source, build_command, output_dir = output_dir,
                num_procs = num_procs, mpi_comm = mpi_comm)

if suite.built: suite.summary()
else: print 'Failed to build/run tests'
