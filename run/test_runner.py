#!/usr/bin/env python
# This is the object that locates all *.enzo_test directories

import imp
import optparse
import os.path
import os
import shutil
import signal
import subprocess
import sys
import time
import tarfile
import logging

known_categories = [
    "Cooling",
    "Cosmology",
    "DrivenTurbulence3D",
    "FLD",
    "GravitySolver",
    "Hydro",
    "MHD",
    "RadiationTransport",
    "RadiationTransportFLD",
]

import numpy
numpy.seterr(all = "ignore")

import nose
from nose.loader import TestLoader
from nose.plugins import Plugin
from nose.plugins import debug
from nose.plugins.manager import PluginManager

from yt.config import ytcfg
ytcfg["yt","suppressStreamLogging"] = "True"
ytcfg["yt","__command_line"] = "True" 
from yt.mods import *

from yt.utilities.command_line import get_yt_version
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.logger import \
    disable_stream_logging, ufstring
disable_stream_logging()
from yt.utilities.answer_testing.framework import \
    AnswerTesting

try:
    yt_version = get_yt_version()
except:
    print "ERROR: cannot get yt version, install yt in develop mode if you want the test runner to record the yt version."
    yt_version = None

# Test keyword types and default values.
varspec = dict(
    name = (str, ''),
    answer_testing_script = (str, None),
    nprocs = (int, 1),
    runtime = (str, 'short'),
    hydro = (bool, False),
    mhd = (bool, False),
    gravity = (bool, False),
    cosmology = (bool, False),
    chemistry = (bool, False),
    cooling = (bool, False),
    AMR = (bool, False),
    dimensionality = (int, 1),
    author = (str, ''),
    max_time_minutes = (float, 1),
    radiation = (str, None),
    quicksuite = (bool, False),
    pushsuite = (bool, False),
    fullsuite = (bool, False),
    problematic = (bool, False)
)

known_variables = dict( [(k, v[0]) for k, v in varspec.items()] )
variable_defaults = dict( [(k, v[1]) for k, v in varspec.items()] )

# Machine specific job scripts.
run_template_dir = 'run_templates'
machines = {'local':       dict(script = 'local.run',
                                command = 'bash'),

	    'local_nompi': dict(script = 'local_nompi.run',
                                command = 'bash'),

            'nics_kraken': dict(script = 'nics_kraken.run',
                                command = 'qsub')}

# Map between job script variables and test keywords.
template_vars = {'N_PROCS'   : 'nprocs',
                 'PAR_FILE'  : 'run_par_file',
                 'TEST_NAME' : 'name',
                 'WALL_TIME' : 'run_walltime'}

results_filename = 'test_results.txt'
version_filename = 'version.txt'

# Files to be included when gathering results.
results_gather = ['results', version_filename]

# If we are able to, let's grab the ~/.enzo/machine_config.py file.
try:
    f, filename, desc = imp.find_module("machine_config",
                          [os.path.expanduser("~/.enzo/")])
    machine_config = imp.load_module("machine_config", f, filename, desc)
except ImportError:
    machine_config = None

def _get_hg_version(path):
    print "Getting current revision."
    from mercurial import hg, ui, commands 
    u = ui.ui() 
    u.pushbuffer() 
    repo = hg.repository(u, path) 
    commands.identify(u, repo) 
    return u.popbuffer()

def _to_walltime(ts):
    return "%02d:%02d:%02d" % \
        ((ts / 3600), ((ts % 3600) / 60), 
         (ts % 60))

def add_files(my_list, dirname, fns):
    my_list += [os.path.join(dirname, fn) for
                fn in fns if fn.endswith(".enzotest")]


def version_swap(repository, changeset, jcompile):
    """Updates *repository* to *changeset*,
    then does make; make -j *jcompile* enzo"""
    print repository, changeset
    from mercurial import hg, ui, commands, util
    options.repository = os.path.expanduser(options.repository)
    u = ui.ui() 
    u.pushbuffer()
    repo = hg.repository(u, options.repository)
    u.popbuffer()
    commands.update(u,repo,changeset)
    command = "cd %s/src/enzo; pwd; "%options.repository
    command += "make clean && make -j %d enzo.exe"%jcompile
    status = os.system(command)
    return status


def bisector(options,args):
    print options.good
    from mercurial import hg, ui, commands, util
    # get current revision
    options.repository = os.path.expanduser(options.repository)
    u = ui.ui() 
    u.pushbuffer()
    repo = hg.repository(u, options.repository)
    u.popbuffer()
    test_directory = os.getcwd()
    command = "cd %s/src/enzo;"%options.repository
    command += "make clean;"
    command += "make -j %d enzo.exe &&"%int(options.jcompile)
    command += "cd %s;"%test_directory
    command += "./test_runner.py --problematic=True " #we only run the problematic tests with bisector
    command += "--output-dir=%s "%options.output_dir
    command += "--repo=%s "%options.repository
    command += "--compare-dir=%s "%options.compare_dir
    def bisection_default_corrector(key,value):
        """mercurial.commands.bisection has bad default values.  This corrects these values."""
        correct_defaults={'good':False,'bad':False,'skip':False,'extend':False,                                                         'command':False,'extra':None,'reset':False}
        correct_defaults[key]=value
        return correct_defaults

    #should read commands.bisect(u,repo,reset=True), etc.,
    # but the defaults for this command do not work. This should be updated
    # when mercurial is.  
    commands.bisect(u,repo,**bisection_default_corrector("reset",True))
    commands.bisect(u,repo,rev=options.good,**bisection_default_corrector("good",True))
    commands.bisect(u,repo,rev=options.bad,**bisection_default_corrector("bad",True))
    commands.bisect(u,repo,**bisection_default_corrector("command",command))

class ResultsSummary(Plugin):
    name = "results_summary"
    score = 10000
    enabled = True

    def options(self, parser, env):
        super(ResultsSummary, self).options(parser, env)

    def configure(self, options, conf):
        super(ResultsSummary, self).configure(options, conf)
        if not self.enabled:
            return
        self.errors = []
        self.failures = []
        self.successes = []

    def addError(self, test, err):
        self.errors.append("%s: ERROR %s" % (test, err))

    def addFailure(self, test, err):
        self.failures.append("%s: FAILURE %s" % (test, err))

    def addSuccess(self, test):
        self.successes.append("%s: PASS" % (test))

    def finalize(self, result, outfile=None):
        print 'Testing complete.'
        print 'Number of errors: %i' % len(self.errors)
        print 'Number of failures: %i' % len(self.failures)
        print 'Number of successes: %i' % len(self.successes)
        if outfile is not None:
            outfile.write('Test Summary\n')
            outfile.write('Tests Passed: %i\n' % len(self.successes))
            outfile.write('Tests Failed: %i\n' % len(self.failures))
            outfile.write('Tests Errored: %i\n' % len(self.errors))
            outfile.write('\n\n')

            outfile.write('Tests that passed: \n')
            for suc in self.successes: 
                outfile.write(suc)
                outfile.write('\n')
            outfile.write('\n')

            outfile.write('Tests that failed:\n')
            for fail in self.failures: 
                outfile.write(fail)
                outfile.write('\n')
            outfile.write('\n')

            outfile.write('Tests that errored:\n')
            for err in self.errors: 
                outfile.write(err)
                outfile.write('\n')
            outfile.write('\n')

class EnzoTestCollection(object):
    def __init__(self, tests = None, verbose=True, args = None,
                 plugins = None):
        self.verbose = verbose
        if args is None: args = sys.argv[:1]
        self.args = args
        if plugins is None: plugins = []
        self.plugins = plugins
        if tests is None:
            # Now we look for all our *.enzotest files
            fns = []
            for cat in known_categories:
                os.path.walk(cat, add_files, fns)
            self.tests = []
            for fn in sorted(fns):
                if self.verbose: print "HANDLING", fn
                self.add_test(fn)
        else:
            self.tests = tests
        self.test_container = []

    def go(self, output_dir, interleaved, machine, exe_path, sim_only=False,
           test_only=False):
        go_start_time = time.time()
        self.output_dir = output_dir
        total_tests = len(self.tests)
        if interleaved:
            for i, my_test in enumerate(self.tests):
                print "Preparing test: %s." % my_test['name']
                self.test_container.append(EnzoTestRun(output_dir, my_test, 
                                                       machine, exe_path,
                                                       args=self.args,
                                                       plugins=self.plugins))
                if not test_only:
                    print "Running simulation: %d of %d." % (i+1, total_tests)
                    self.test_container[i].run_sim()
                if not sim_only:
                    print "Running test: %d of %d." % (i+1, total_tests)
                    self.test_container[i].run_test()
        else:
            self.prepare_all_tests(output_dir, machine, exe_path)
            if not test_only: self.run_all_sims()
            if not sim_only: self.run_all_tests()
        if not sim_only: self.save_test_summary()
        go_stop_time = time.time()
        print "\n\nComplete!"
        print "Total time: %f seconds." % (go_stop_time - go_start_time)
        print "See %s/%s for a summary of all tests." % \
            (self.output_dir, results_filename)

    def prepare_all_tests(self, output_dir, machine, exe_path):
        print "Preparing all tests."
        for my_test in self.tests:
            print "Preparing test: %s." % my_test['name']
            self.test_container.append(EnzoTestRun(output_dir, my_test, 
                                                   machine, exe_path,
                                                   args = self.args,
                                                   plugins = self.plugins))

    def run_all_sims(self):
        total_tests = len(self.test_container)
        print "Running all simulations."
        for i, my_test in enumerate(self.test_container):
            print "Running simulation: %d of %d." % (i+1, total_tests)
            my_test.run_sim()

    def run_all_tests(self):
        total_tests = len(self.test_container)
        print "Running all tests."
        for i, my_test in enumerate(self.test_container):
            print "Running test: %d of %d." % (i+1, total_tests)
            my_test.run_test()

    def add_test(self, fn):
        # We now do something dangerous: we exec the file directly and grab
        # its environment variables from it.
        local_vars = {}
        execfile(fn, {}, local_vars)
        test_spec = variable_defaults.copy()
        test_spec['fullpath'] = fn
        test_spec['fulldir'] = os.path.dirname(fn)
        test_spec['run_par_file'] = os.path.basename(test_spec['fulldir']) + ".enzo"
        test_spec['run_walltime'] = _to_walltime(60 * test_spec['max_time_minutes'] * 
                                                 options.time_multiplier)
        for var, val in local_vars.items():
            if var in known_variables:
                caster = known_variables[var]
                if val == "False": val = False
                test_spec[var] = caster(val)
                if val == "None": test_spec[var] = None
            else:
                print "%s UNRECOGNIZED VARIABLE %s" % ( fn, var)
        self.tests.append(test_spec)

    def unique(self, param):
        pp = set()
        for t in self.tests:
            pp.add(t.get(param, "Key Missing"))
        return pp

    def params(self):
        pp = set()
        for t in self.tests:
            pp.update(set(t.keys()))
        return pp

    def select(self, **kwargs):
        pp = []
        for t in self.tests:
            include = True
            for param, value in kwargs.items():
                if value == "None": value = None
                if value == "False": value = False
                if t.get(param, "Key Missing") != value:
                    include = False
                    break
            if include == True: pp.append(t)
        #raise RuntimeError
        return EnzoTestCollection(tests = pp, args = self.args,
                                  verbose = self.verbose,
                                  plugins = self.plugins)

    def summary(self):
        for param in sorted(self.params()):
            if param.startswith("full"): continue
            print param
            for v in self.unique(param):
                print "     %s" % (v)
            print
        print
        print "NUMBER OF TESTS", len(self.tests)

    def save_test_summary(self):
        all_passes = all_failures = 0
        run_passes = run_failures = 0
        dnfs = default_test = 0
        f = open(os.path.join(self.output_dir, results_filename), 'w')
        self.plugins[1].finalize(None, outfile=f)
        # for my_test in self.test_container:
        #     default_only = False
        #     if my_test.run_finished:
        #         if my_test.test_data['answer_testing_script'] == 'None' or \
        #                 my_test.test_data['answer_testing_script'] is None:
        #             default_only = True
        #             default_test += 1
        #         t_passes = 0
        #         t_failures = 0
        #         for t_result in my_test.results.values():
        #             t_passes += int(t_result)
        #             t_failures += int(not t_result)
        #         f.write("%-70sPassed: %4d, Failed: %4d" % (my_test.test_data['fulldir'], 
        #                                                    t_passes, t_failures))
        #         if default_only:
        #             f.write(" (default tests).\n")
        #         else:
        #             f.write(".\n")
        #         all_passes += t_passes
        #         all_failures += t_failures
        #         run_passes += int(not (t_failures > 0))
        #         run_failures += int(t_failures > 0)
        #     else:
        #         dnfs += 1
        #         f.write("%-70sDID NOT FINISH\n" % my_test.test_data['fulldir'])

        # f.write("\n")
        # f.write("%-70sPassed: %4d, Failed: %4d.\n" % ("Total", 
        #                                               all_passes, all_failures))
        # f.write("Runs finished with all tests passed: %d.\n" % run_passes)
        # f.write("Runs finished with at least one failure: %d.\n" % run_failures)
        # f.write("Runs failed to complete: %d.\n" % dnfs)
        # f.write("Runs finished with only default tests available: %d.\n" % default_test)
        f.close()
        if all_failures > 0 or dnfs > 0:
            self.any_failures = True
        else:
            self.any_failures = False

class EnzoTestRun(object):
    def __init__(self, test_dir, test_data, machine, exe_path,
                 args = None, plugins = None):
        if args is None: args = sys.args[:1]
        self.args = args
        if plugins is None: plugins = []
        self.plugins = plugins
        self.machine = machine
        self.test_dir = test_dir
        self.test_data = test_data
        self.exe_path = exe_path
        self.results = {}
        if self.exe_path is None:
            self.local_exe = None
        else:
            self.local_exe = os.path.basename(exe_path)

        self.run_dir = os.path.join(self.test_dir, self.test_data['fulldir'])

        self._copy_test_files()
        self._create_run_script()

    def _copy_test_files(self):
        # Check for existence
        if os.path.exists(self.run_dir):
            if options.clobber:
                print "%s exists, but clobber == True, so overwriting it." % self.test_data['name']
                shutil.rmtree(self.run_dir)
                shutil.copytree(self.test_data['fulldir'], self.run_dir)
                # Copy version file into run directory
                shutil.copy(os.path.join(self.test_dir, version_filename),
                            os.path.join(self.run_dir, version_filename))
                if self.exe_path is not None:
                    os.symlink(os.path.realpath(self.exe_path), 
                               os.path.join(self.run_dir, self.local_exe))
            else:
                print "%s already exists. Skipping directory." % self.test_data['name']
        else:
            shutil.copytree(self.test_data['fulldir'], self.run_dir)
            # Copy version file into run directory
            shutil.copy(os.path.join(self.test_dir, version_filename),
                        os.path.join(self.run_dir, version_filename))
            if self.exe_path is not None:
                os.symlink(os.path.realpath(self.exe_path), 
                           os.path.join(self.run_dir, self.local_exe))

    def _create_run_script(self):
        template_path = os.path.join(os.path.dirname(__file__), 
                                     run_template_dir, machines[self.machine]['script'])
        template_dest = os.path.join(self.run_dir, machines[self.machine]['script'])
        f = open(template_path, 'r')
        template = f.read()
        f.close()
        for var in template_vars.keys():
            template = template.replace(('${%s}' % var), 
                                        str(self.test_data[template_vars[var]]))
        template = template.replace('${EXECUTABLE}', "./%s" % self.local_exe)
        for var, value in sorted(getattr(machine_config, self.machine, {}).items()):
            template = template.replace('${%s}' % var, value)
        f = open(template_dest, 'w')
        f.write(template)
        f.close()

    def run_sim(self):
        print "Running test simulation: %s." % self.test_data['fulldir']
        cur_dir = os.getcwd()
        # Check for existence
        if os.path.exists(os.path.join(self.run_dir, 'RunFinished')):
            print "%s run already completed, continuing..." % self.test_data['name']
            return
        
        os.chdir(self.run_dir)
        command = "%s %s" % (machines[self.machine]['command'], 
                             machines[self.machine]['script'])
        sim_start_time = time.time()
        # Run the command.
        proc = subprocess.Popen(command, shell=True, close_fds=True, 
                                preexec_fn=os.setsid)

        print "Simulation started on %s with maximum run time of %d seconds." % \
            (time.ctime(), (self.test_data['max_time_minutes'] * 60 *
                            options.time_multiplier))
        running = 0
        # Kill the script if the max run time exceeded.
        while proc.poll() is None:
            if running > (self.test_data['max_time_minutes'] * 60 *
                          options.time_multiplier):
                print "Simulation exceeded maximum run time."
                os.killpg(proc.pid, signal.SIGUSR1)
            running += 1
            time.sleep(1)
        
        sim_stop_time = time.time()
        if os.path.exists(os.path.join(self.run_dir, 'RunFinished')):
            f = open(os.path.join(self.run_dir, 'run_time'), 'w')
            f.write("%f seconds.\n" % (sim_stop_time - sim_start_time))
            f.close()
            print "Simulation completed in %f seconds." % \
                (sim_stop_time - sim_start_time)
        os.chdir(cur_dir)

    def run_test(self):
        rf = os.path.join(self.run_dir, 'RunFinished')
        self.run_finished = os.path.exists(rf)
        tl = TestLoader()
        tl.config.plugins = PluginManager(plugins = self.plugins)
        suite = tl.loadTestsFromDir(self.run_dir)
        nose.run(argv=self.args, suite=suite)

class UnspecifiedParameter(object):
    pass
unknown = UnspecifiedParameter()

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("--clobber", dest='clobber', default=False,
                      action="store_true", 
                      help="Recopies tests and tests from scratch.")
    parser.add_option("--interleave", action='store_true', dest='interleave', 
                      default=False,
                      help="Option to interleave preparation, running, and testing.")
    parser.add_option("-m", "--machine", dest='machine', default='local', metavar='str',
                      help="Machine to run tests on.")
    parser.add_option("-o", "--output-dir", dest='output_dir', metavar='str',
                      help="Where to place the run directory")
    parser.add_option("--repo", dest='repository', default="../",
                      help="Path to repository being tested.", metavar='str')
    parser.add_option("--sim-only", dest='sim_only', action="store_true", 
                      default=False, help="Only run simulations.")
    parser.add_option("--test-only", dest='test_only', action="store_true", 
                      default=False, help="Only perform tests.")
    parser.add_option("--time-multiplier", dest='time_multiplier',
                      default=1.0, type=float, metavar='int',
                      help="Multiply simulation time limit by this factor.")
    parser.add_option("-v", "--verbose", dest='verbose', action="store_true",
                      default=False, help="Slightly more verbose output.")
    parser.add_option("--pdb", action="store_true", dest="pdb",
                      default=False, help="Drop into debugger on errors")
    parser.add_option("-b", "--bisect", dest="bisect", action="store_true",
                      default=False, help="Run bisection on test. Requires revisions" +
                      "--good and --bad.  Best if --repo is different from location of test_runner.py."+
                      "Runs  --problematic suite.  See README for more info")
    parser.add_option("--good", dest="good", default=None, metavar='str', 
                      help="For bisection, most recent good revision")
    parser.add_option("--bad", dest="bad", default=None, metavar='str',
                      help="For bisection, most recent bad revision")
    parser.add_option("-j", "--jcompile", dest="jcompile", type="int", default=1, 
                      metavar='int',
                      help="number of processors with which to compile when running bisect")
    parser.add_option("--changeset", dest="changeset", default=None, metavar='str',
                      help="Changeset to use in simulation repo.  If supplied, make clean && make is also run")
    parser.add_option("--run-suffix", dest="run_suffix", default=None, metavar='str',
                      help="An optional suffix to append to the test run directory. Useful to distinguish multiple runs of a given changeset.")

    answer_plugin = AnswerTesting()
    answer_plugin.enabled = True
    answer_plugin.options(parser)
    reporting_plugin = ResultsSummary()
    reporting_plugin.enabled = True
    pdb_plugin = debug.Pdb()

    all_suites = ['quick', 'push', 'full']
    suite_vars = [suite+"suite" for suite in all_suites]
    testproblem_group = optparse.OptionGroup(parser, "Test problem selection options")
    testproblem_group.add_option("", "--suite",
                                 dest="test_suite", default=unknown,
                                 help="quick: 37 tests in ~25 minutes, push: 48 tests in ~90 minutes, full: 96 tests in ~18 hours.",
                                 choices=all_suites, metavar=all_suites)

    for var, caster in sorted(known_variables.items()):
        if var not in suite_vars:
            testproblem_group.add_option("", "--%s" % (var),
                                         type=str, default = unknown,
                                         metavar=caster.__name__)
    parser.add_option_group(testproblem_group)
    options, args = parser.parse_args()

    if options.pdb:
        pdb_plugin.enabled = True
        pdb_plugin.enabled_for_failures = True

    # Get information about the current repository, set it as the version in
    # the answer testing plugin.
    options.repository = os.path.expanduser(options.repository)
    hg_current = _get_hg_version(options.repository)
    rev_hash = hg_current.split()[0]

    if options.run_suffix:
        rev_hash += options.run_suffix

    answer_plugin._my_version = rev_hash

    answer_plugin.configure(options, None)
    reporting_plugin.configure(options, None)

    # Break out if output directory not specified.
    if options.output_dir is None:
        print 'Please enter an output directory with -o option'
        sys.exit(1)

    if options.changeset is not None:
        status = version_swap(options.repository, options.changeset, options.jcompile)
        if status:
            sys.exit(status)

    if options.bisect:
        bisector(options,args)
        sys.exit(0)
    etc = EnzoTestCollection(verbose=options.verbose, args=sys.argv[:1],
                             plugins = [answer_plugin, reporting_plugin, pdb_plugin])

    construct_selection = {}
    if options.test_suite != unknown:
        suite_var = str(options.test_suite) + "suite"
        print suite_var
        construct_selection[suite_var] = \
          known_variables[suite_var](options.test_suite)
    for var, caster in known_variables.items():
        if var in suite_vars: continue
        if getattr(options, var) != unknown:
            val = getattr(options, var)
            if val == 'None': val = None
            if val == "False": val = False
            construct_selection[var] = caster(val)
    print
    print "Selecting with:"
    for k, v in sorted(construct_selection.items()):
        print "     %s = %s" % (k, v)
    etc2 = etc.select(**construct_selection)
    print
    print "\n".join(list(etc2.unique('name')))
    print "Total: %s" % len(etc2.tests)

    # Gather results and version files for all test and tar them.
    # get current revision
    options.output_dir = os.path.join(options.output_dir, rev_hash)

    if not os.path.exists(options.output_dir): os.makedirs(options.output_dir)
    f = open(os.path.join(options.output_dir, version_filename), 'w')
    f.write('Enzo: %s' % hg_current)
    f.write('yt: %s\n' % yt_version)
    f.close()

    # the path to the executable we're testing
    exe_path = os.path.join(options.repository, "src/enzo/enzo.exe")

    # Make it happen
    etc2.go(options.output_dir, options.interleave, options.machine, exe_path,
            sim_only=options.sim_only, 
            test_only=options.test_only)

    # Store the results locally or in the cloud.
    answer_plugin.finalize()
    #reporting_plugin.finalize(None, res_file = )

    try:
        import json
    except ImportError:
        json = None
    if json is not None:
        f = open("results.js", "w")
        results = []
        for test in etc2.test_container:
            # This is to avoid any sorting code in JS
            vals = test.results.items()
            vals.sort()
            results.append( dict(name = test.test_data['name'],
                             results = vals) )
        f.write("test_data = %s;\n" % (json.dumps(results, indent=2)))
        f.write("current_set = '%s';\n" % (hg_current.strip()))

    if etc2.any_failures:
        sys.exit(1)
    else:
        sys.exit(0)
