if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <GPU_NUMBER>"
    exit 1
fi

# Assign the provided arguments to variables
GPU_NUMBER="$1"

# Build the container image
podman build . -t sm_plants_precursors

# Run the container with the specified Python script as the command
podman run -v /home/jcapela/SMPrecursorPredictor/models_and_datasets/version_3/:/workspace/:z -d --device nvidia.com/gpu="$GPU_NUMBER" --security-opt=label=disable --name="$FOLDER_NAME" sm_plants_precursors /bin/bash -c "cd /workspace/ && python train_models.py > output.txt"
