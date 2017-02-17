"""
Model Intercomparison Study problem 4
"""

import os

from credo.systest import SciBenchmarkTest
from credo.systest import FieldWithinTolTC
from credo.systest import HistoryWithinTolTC

from credo.jobrunner import SimpleJobRunner
from credo.modelresult import ModelResult, HistoryDataResult
from credo.t2model import T2ModelRun, T2ModelResult
from credo.waiwera import WaiweraModelRun

import credo.reporting.standardReports as sReps
from credo.reporting import getGenerators

from mulgrids import mulgrid
import matplotlib.pyplot as plt
import numpy as np

from docutils.core import publish_file

model_name = 'problem4'

WAIWERA_FIELDMAP = {
    'Pressure': 'fluid_pressure',
    'Temperature': 'fluid_temperature',
    'Liquid saturation': 'fluid_liquid_saturation',
}

model_dir = './run'
data_dir = './data'
t2geo_filename = os.path.join(model_dir, 'g' + model_name + '.dat')

num_procs = 1

run_name = 'run'
run_index = 0
test_fields = ["Pressure", "Temperature", "Liquid saturation"]
plot_fields = test_fields
digitised_test_fields = {550: ["Pressure"],
                         1050: ["Pressure", "Liquid saturation"],
                         1550: ["Pressure", "Liquid saturation"],
                         1950: ["Pressure", "Liquid saturation"]}

geo = mulgrid(t2geo_filename)
map_out_atm = range(geo.num_atmosphere_blocks, geo.num_blocks)

problem4_test = SciBenchmarkTest(model_name + "_test", nproc = num_procs)
problem4_test.description = """Model Intercomparison Study problem 4"""

run_base_name = model_name
run_filename = run_base_name + '.json'
model_run = WaiweraModelRun(run_name, run_filename,
                          fieldname_map = WAIWERA_FIELDMAP,
                          simulator = 'waiwera',
                          basePath = os.path.realpath(model_dir))
model_run.jobParams['nproc'] = num_procs
problem4_test.mSuite.addRun(model_run, run_name)

problem4_test.setupEmptyTestCompsList()

run_base_name = model_name
run_filename = os.path.join(model_dir, run_base_name + ".listing")
AUTOUGH2_result = T2ModelResult("AUTOUGH2", run_filename,
                                 geo_filename = t2geo_filename,
                                 ordering_map = map_out_atm)

depths = [550, 1050, 1550, 1950]
obs_cell_index = {}
for depth in depths:
    obs_position = np.array([50., 50., -depth])
    obs_blk = geo.block_name_containing_point(obs_position)
    obs_cell_index[depth] = geo.block_name_index[obs_blk] - geo.num_atmosphere_blocks

data = {}
ref_sim = 'S-Cubed'
for depth in depths:
    for field in digitised_test_fields[depth]:
        obspt = str(depth)
        data_filename = '_'.join((model_name, obspt, field, ref_sim))
        data_filename = data_filename.lower().replace(' ', '_')
        data_filename = os.path.join(data_dir, data_filename + '.dat')
        data[field, obs_cell_index[depth]] = np.loadtxt(data_filename)
ref_result = HistoryDataResult(ref_sim, data)

for depth in depths:
    obspt = str(depth)
    problem4_test.addTestComp(run_index, "AUTOUGH2 z = " + obspt,
                          HistoryWithinTolTC(fieldsToTest = test_fields,
                                             defFieldTol = 2.e-3,
                                             expected = AUTOUGH2_result,
                                             testCellIndex = obs_cell_index[depth]))

    problem4_test.addTestComp(run_index, ref_sim + " z = " + obspt,
                              HistoryWithinTolTC(fieldsToTest = digitised_test_fields[depth],
                                                 defFieldTol = 2.e-2,
                                                 expected = ref_result,
                                                 testCellIndex = obs_cell_index[depth],
                                                 orthogonalError = True))

jrunner = SimpleJobRunner(mpi = True)
testResult, mResults = problem4_test.runTest(jrunner, createReports = True)

# plots:
scale = {"Pressure": 1.e6, "Temperature": 1., "Liquid saturation": 1.}
unit = {"Pressure": "MPa", "Temperature": "$^{\circ}$C", "Liquid saturation": ""}
symbol = {ref_sim: 'o'}
yr = 365. * 24. * 60. * 60.

for depth in depths:
    obspt = str(depth)
    obs_position = np.array([50., 50., -depth])
    obs_blk = geo.block_name_containing_point(obs_position)
    obs_cell_index = geo.block_name_index[obs_blk] - geo.num_atmosphere_blocks
    tc_name = "AUTOUGH2 z = " + obspt

    for field_name in digitised_test_fields[depth]:

        t, var = problem4_test.mSuite.resultsList[run_index].\
                 getFieldHistoryAtCell(field_name, obs_cell_index)
        plt.plot(t / yr, var / scale[field_name], '-', label = 'Waiwera')

        t, var = AUTOUGH2_result.getFieldHistoryAtCell(field_name, obs_cell_index)
        plt.plot(t / yr, var / scale[field_name], '+', label = 'AUTOUGH2')

        t, var = ref_result.getFieldHistoryAtCell(field_name, obs_cell_index)
        plt.plot(t / yr, var / scale[field_name], symbol[ref_sim], label = ref_sim)

        plt.xlabel('time (years)')
        plt.ylabel(field_name + ' (' + unit[field_name] + ')')
        plt.legend(loc = 'best')
        plt.title(' '.join((model_name, field_name.lower(),
                            'results at depth', obspt, 'm')))
        img_filename_base = '_'.join((model_name, obspt, field_name))
        img_filename_base = img_filename_base.replace(' ', '_')
        img_filename = os.path.join(problem4_test.mSuite.runs[run_index].basePath,
                                    problem4_test.mSuite.outputPathBase,
                                    img_filename_base + '.png')
        plt.tight_layout(pad = 3.)
        plt.savefig(img_filename)
        plt.clf()
        problem4_test.mSuite.analysisImages.append(img_filename)

# generate report:

for rGen in getGenerators(["RST"], problem4_test.outputPathBase):
    report_filename = os.path.join(problem4_test.outputPathBase,
                     "%s-report.%s" % (problem4_test.testName, rGen.stdExt))
    sReps.makeSciBenchReport(problem4_test, mResults, rGen, report_filename)
    html_filename = os.path.join(problem4_test.outputPathBase,
                     "%s-report.%s" % (problem4_test.testName, 'html'))
    html = publish_file(source_path = report_filename,
                        destination_path = html_filename,
                        writer_name = "html")

