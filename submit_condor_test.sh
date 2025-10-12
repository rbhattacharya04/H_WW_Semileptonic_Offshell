#!/bin/bash
# filepath: /afs/cern.ch/work/s/sverma/private/SEMILEP_FRAMEWORK/H_WW_Semileptonic_Offshell/submit_condor.sh

# Default values
DATASET=""
OUTPUT_DIR="output"
HOURS=10

# Help function
show_help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -d DATASET_OR_FILE    Dataset name or single file path (required)"
    echo "  -o DIR                Output directory (default: output)"
    echo "  -t HOURS              Wall time in hours (default: 10)"
    echo "  -h                    Show this help"
    echo "Example: $0 -d ggH_sonly_off -o results -t 12"
    echo "         $0 -d /path/to/file.root -o results -t 12"
    exit 1
}

# Parse arguments
while getopts "d:o:t:h" opt; do
    case $opt in
        d) DATASET=$OPTARG ;;
        o) OUTPUT_DIR=$OPTARG ;;
        t) HOURS=$OPTARG ;;
        h) show_help ;;
        *) show_help ;;
    esac
done

# Check required arguments
if [ -z "$DATASET" ]; then
    echo "Error: Dataset (-d) or file path is required"
    show_help
fi

# Setup
MAX_RUNTIME=$((HOURS * 3600))
mkdir -p $OUTPUT_DIR condor_logs

# Sanitize dataset/file name for submission/log file names
SAFE_NAME=$(basename "$DATASET")
SAFE_NAME=${SAFE_NAME//\//_}

# Create HTCondor submission file (one job per dataset or file)
cat > job_${SAFE_NAME}.sub << EOF
universe = vanilla
executable = run_job.sh
arguments = $DATASET $OUTPUT_DIR
output = condor_logs/${SAFE_NAME}.out
error  = condor_logs/${SAFE_NAME}.err
log    = condor_logs/${SAFE_NAME}.log
request_cpus = 1
request_memory = 2GB
request_disk = 5GB
+MaxRuntime = $MAX_RUNTIME
transfer_input_files = Dataset.py,Analysis.py,SemiLeptonic.h
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
getenv = True
queue 1
EOF

# Create job script
cat > run_job.sh << 'EOF'
#!/bin/bash
DATASET=$1
OUTPUT_DIR=$2

export PYTHONPATH=.:$PYTHONPATH
mkdir -p $OUTPUT_DIR

python3 Analysis.py $DATASET
EOF

chmod +x run_job.sh
condor_submit job_${SAFE_NAME}.sub
echo "Submitted 1 job for $DATASET with $HOURS hours runtime"
echo "Check job status with: condor_q"