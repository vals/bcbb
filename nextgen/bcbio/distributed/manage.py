"""Manage processes on a cluster

Automate:
 - starting working nodes to process the data
 - kicking off an analysis
 - cleaning up nodes on finishing

Currently works on LSF and SGE managed clusters; it's readily generalizable to
other architectures as well.
"""
import time
import math
import subprocess

from bcbio.log import logger2 as log


def run_and_monitor(config, config_file, args, workers_needed=None,
                    task_module=None, queues=None):
    """Run a distributed analysis in s cluster environment, monitoring outputs.
    """
    cp = config["distributed"]["cluster_platform"]
    cluster = __import__("bcbio.distributed.{0}".format(cp), fromlist=[cp])
    jobids = []
    try:
        # If the manager is going to be run on a cluster, submit the job
        manager_id = None
        local_manager = config["distributed"].get("run_process_program_locally", False)
        if not local_manager:
            log.info("Starting manager")
            manager_id = start_analysis_manager(cluster, args, config)
            jobids.append(manager_id)

        log.info("Starting cluster workers")
        jobids.extend(start_workers(cluster, config, config_file, workers_needed,
                                    task_module, queues))
        while not(cluster.are_running(jobids)):
            time.sleep(5)

        # If manager should run locally, run it as a regular subprocess rather than monitor its status on the cluster
        log.info("Running analysis")
        if not local_manager:
            monitor_analysis(cluster, manager_id)

        else:
            start_analysis_manager(cluster, args, config)

    finally:
        log.info("Cleaning up cluster workers")
        stop_workers(cluster, jobids)


def start_workers(cluster, config, config_file, workers_needed=None,
                  task_module=None, queues=None):
    """Initiate worker nodes on cluster, returning jobs IDs for management.
    """
    # we can manually specify workers or dynamically get as many as needed
    num_workers = config["distributed"].get("num_workers", None)
    if num_workers in [None, "all"]:
        cores_per_host = config["distributed"].get("cores_per_host", 1)
        if cores_per_host == 0:
            raise ValueError("Set num_workers or cores_per_host in YAML config")

        assert workers_needed is not None, \
               "Supply workers needed if not configured in YAML"
        num_workers = int(math.ceil(float(workers_needed) / cores_per_host))

    program_cl = [config["analysis"]["worker_program"], config_file]
    if task_module:
        program_cl.append("--tasks={0}".format(task_module))

    if queues:
        program_cl.append("--queues={0}".format(queues))

    args = config["distributed"]["platform_args"].split()

    return [cluster.submit_job(args, program_cl) for _ in range(num_workers)]


def start_analysis_manager(cluster, args, config):
    """Start analysis manager node on cluster.
    """
    cluster_args = config["distributed"]["platform_args"].split()
    program_cl = [config["analysis"]["process_program"]] + args
    if not config["distributed"].get("run_process_program_locally", False):
        job_id = cluster.submit_job(cluster_args, program_cl)

    else:
        cl = [str(i) for i in program_cl]
        subprocess.check_call(cl)
        job_id = None

    # wait for job to start
    # Avoid this for systems where everything queues as batches
    #while not(cluster.are_running([job_id])):
    #    time.sleep(5)
    return job_id


def monitor_analysis(cluster, job_id):
    """Wait for manager cluster job to finish
    """
    while cluster.are_running([job_id]):
        time.sleep(5)


def stop_workers(cluster, jobids):
    for jobid in jobids:
        try:
            cluster.stop_job(jobid)
        except:
            pass
