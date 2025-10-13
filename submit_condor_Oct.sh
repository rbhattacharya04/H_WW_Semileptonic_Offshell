#!/bin/bash
# filepath: /afs/cern.ch/work/s/sverma/private/SEMILEP_FRAMEWORK/H_WW_Semileptonic_Offshell/submit_condor.sh
# HTCondor submission script with argument parsing

# Default values
DATASET=""
OUTPUT_DIR="output"
JOBS=10
HOURS=10

# Help function
show_help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -d DATASET    Dataset to process (required)"
    echo "  -o DIR        Output directory (default: output)"
    echo "  -j NUM        Number of jobs (default: 10)"
    echo "  -t HOURS      Wall time in hours (default: 10)"
    echo "  -h            Show this help"
    echo "Example: $0 -d ggH_sonly_off -o results -j 20 -t 12"
    exit 1
}

# Parse arguments
while getopts "d:o:j:t:h" opt; do
    case $opt in
        d) DATASET=$OPTARG ;;
        o) OUTPUT_DIR=$OPTARG ;;
        j) JOBS=$OPTARG ;;
        t) HOURS=$OPTARG ;;
        h) show_help ;;
        *) show_help ;;
    esac
done

# Check required arguments
if [ -z "$DATASET" ]; then
    echo "Error: Dataset (-d) is required"
    show_help
fi

# Setup
MAX_RUNTIME=$((HOURS * 3600))
mkdir -p $OUTPUT_DIR condor_logs

# Create HTCondor submission file
cat > job_${DATASET}.sub << EOF
universe = vanilla
executable = run_job.sh
arguments = \$(Process) $DATASET $OUTPUT_DIR $JOBS
output = condor_logs/${DATASET}_\$(Process).out
error = condor_logs/${DATASET}_\$(Process).err
log = condor_logs/${DATASET}_\$(Process).log
request_cpus = 1
request_memory = 2GB
request_disk = 5GB
+MaxRuntime = $MAX_RUNTIME
transfer_input_files = Dataset.py,Analysis.py,SemiLeptonic.h
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
getenv = True
queue $JOBS
EOF

# Create job script
cat > run_job.sh << 'EOF'
#!/bin/bash
JOB_ID=$1
DATASET=$2
OUTPUT_DIR=$3
JOBS=$4

# Ensure Python finds modules in the current directory
export PYTHONPATH=.:$PYTHONPATH

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

python3 -c "
import sys
import os
from Dataset import dataset
import Analysis
import pickle

# Ensure output directory exists
os.makedirs('$OUTPUT_DIR', exist_ok=True)

files = dataset['$DATASET']
chunk_size = len(files) // $JOBS + 1
start = $JOB_ID * chunk_size
end = min(start + chunk_size, len(files))
chunk = files[start:end]

print(f'Processing {len(chunk)} files for job $JOB_ID')
results = Analysis.makeRDF(chunk, True)

output_file = '$OUTPUT_DIR/histograms_${DATASET}_$JOB_ID.pkl'
print(f'Saving results to {output_file}')
with open(output_file, 'wb') as f:
    pickle.dump(results, f)

print(f'Job $JOB_ID completed successfully')
"
EOF

chmod +x run_job.sh
condor_submit job_${DATASET}.sub
echo "Submitted $JOBS jobs for dataset $DATASET with $HOURS hours runtime"
echo "Check job status with: condor_q"
