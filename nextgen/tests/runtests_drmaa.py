#!/usr/bin/env python
"""
Runs nosetests in a DRMAA-compatible cluster environment.
"""
import drmaa
import os
import yaml

def main():
    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.remoteCommand = 'nosetests'
    #Run the standard tests
    jt.args = ['-v', '-s', '--with-xunit', '-a', 'standard']

    jt.jobName = "testing_pipeline"
    jt.workingDirectory = drmaa.JobTemplate.HOME_DIRECTORY+'/opt/bcbb/nextgen/tests'
    jt.outputPath = ":"+drmaa.JobTemplate.HOME_DIRECTORY+'/opt/bcbb/nextgen/tests/tests_results.out'
    jt.nativeSpecification = "-A a2010002 -p devel -t 00:30:00"

    jobid = s.runJob(jt)
    print 'Your job has been submitted with id ' + jobid
    print 'The results are being written in ~/opt/bcbb/nextgen/tests/tests_results.out'

    print 'Cleaning up'
    s.deleteJobTemplate(jt)
    s.exit()
    exit(0)

if __name__ == "__main__":
    main()

