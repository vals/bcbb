#!/usr/bin/env python
"""
This script runs nosetests in a DRMAA-compatible cluster environment.
"""
import drmaa
import os
import yaml

def main():
    s = drmaa.Session()
    s.initialize()

    jt = s.createJobTemplate()
    jt.remoteCommand = 'nosetests'
    jt.args = ['-v', '-s', '--with-xunit']

    jt.job_name = "nosetests"
    #jt.nativeSpecification = config[distributed][platform_args]
    jt.nativeSpecification = "-A a2010002 -p devel -t 00:30:00"

    jobid = s.runJob(jt)
    print 'Your job has been submitted with id ' + jobid

    s.deleteJobTemplate(jt)
    s.exit()
    exit(0)

if __name__ == "__main__":
    main()
