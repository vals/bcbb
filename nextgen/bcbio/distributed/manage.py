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


def run_on_slurm(config, config_file, args, workers_needed=None,
                    task_module=None, queues=None):
    "Run distributed analysis in one SLURM job with multiple nodes."
    settings = config["distributed"]["platform_args"]
    assert type(settings) == dict, \
    "SLURM jobs needs a dict of platform arguments, unlike other platforms."

    num_workers = config["distributed"].get("num_workers", None)
    if num_workers in [None, "all"]:
        # TODO: Assertions and such
        cores_per_host = config["distributed"].get("cores_per_host", 1)
        num_workers = int(math.ceil(float(workers_needed) / cores_per_host))

    # Provide filenames for files to send to the batch job system.
    srun_script_name = "srun.sh"
    sbatch_script_name = "run_srun.sh"

    # Create a shell script which executes different programs depending on
    # which node the program is running.
    analysis = config["analysis"]
    manager_case = "0) {process_program} ".format(**analysis) + args
    worker_case = "%i) {worker_program}".format(**analysis) + config_file

    with open(srun_script_name, "w") as script:
        script.write("#!/bin/sh\n\n")
        script.write(r"case $SLURM_NODEID in")
        script.write("\n")
        script.write("  " + manager_case + " ;;\n")
        for i in range(num_workers):
            script.write("  " + worker_case % (i + 1,) + " ;;\n")
        script.write("esac\n")

    subprocess.call(["chmod", "+x", srun_script_name])

    # Create an sbatch script to be sent to SLURM
    with open(sbatch_script_name, "w") as script:
        script.write("#!/bin/sh\n")
        s = "#SBATCH"
        script.write(" ".join([s, "-N", str(num_workers + 1)]) + "\n")
        for prop, arg in settings.items():
            script.write(" ".join([s, prop, arg]) + "\n")
        script.write("\n" + " ".join(["srun", srun_script_name]) + "\n")

    # Send things to the job system
    subprocess.call(["sbatch", sbatch_script_name])


def run_and_monitor(config, config_file, args, workers_needed=None,
                    task_module=None, queues=None):
    """Run a distributed analysis in s cluster environment, monitoring outputs.
    """
    cp = config["distributed"]["cluster_platform"]
    cluster = __import__("bcbio.distributed.{0}".format(cp), fromlist=[cp])
    jobids = []
    try:
        print "Starting manager"
        manager_id = start_analysis_manager(cluster, args, config)
        print "Starting cluster workers"
        jobids.extend(start_workers(cluster, config, config_file, workers_needed,
                                    task_module, queues))
        jobids.append(manager_id)
        while not(cluster.are_running(jobids)):
            time.sleep(5)
        print "Running analysis"
        monitor_analysis(cluster, manager_id)
    finally:
        print "Cleaning up cluster workers"
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
    job_id = cluster.submit_job(cluster_args, program_cl)
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
