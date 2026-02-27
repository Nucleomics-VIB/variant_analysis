#!/bin/bash

cleanup() {
    echo "Cleaning up..."
    echo $(date) >> slurm_resources.txt
    
    sacct -j $SLURM_JOB_ID --format=JobID,JobName,State,ExitCode,Elapsed,CPUTime,TotalCPU,ReqCPUS,AllocCPUS,AveRSS,MaxRSS,ReqMem,AveVMSize,MaxVMSize \
      >> resources.txt
    
    echo "Cleanup complete. Exiting."
    exit 0
}

cleanup

exit 0
