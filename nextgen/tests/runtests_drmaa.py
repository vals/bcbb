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
    jt.args = ['-v', '-s', '--with-xunit', '-a', 'standard']

    jt.jobName = "testing_pipeline"
    jt.workingDirectory = drmaa.JobTemplate.HOME_DIRECTORY+'/opt/bcbb/nextgen/tests'
    jt.outputPath = ":"+drmaa.JobTemplate.HOME_DIRECTORY+'/opt/bcbb/nextgen/tests/tests_results.out'
    jt.nativeSpecification = "-A a2010002 -p node -t 04:00:00"
    jt.email=["guillermo.carrasco@scilifelab.se"]

    jobid = s.runJob(jt)
    print 'Your job has been submitted with id ' + jobid
    print 'The results are being written in ~/opt/bcbb/nextgen/tests/tests_results.out'
    retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
    print 'Job: ' + str(retval.jobId) + ' finished with status ' + str(retval.hasExited)

    print 'Cleaning up'
    s.deleteJobTemplate(jt)
    s.exit()
    exit(0)

if __name__ == "__main__":
    main()

