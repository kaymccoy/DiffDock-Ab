NUM_FOLDS=1  # number of seeds to try, default 5
SEED=0  # initial seed
CUDA=0  # will use GPUs from CUDA to CUDA + NUM_GPU - 1
NUM_GPU=1
BATCH_SIZE=1  # split across all GPUs
NUM_SAMPLES=3

NAME="db5_esm_noST"  # change to name of config file
RUN_NAME="db5_esm_noST"
CONFIG="config/${NAME}.yaml"

SAVE_PATH="ckpts/${RUN_NAME}"
VISUALIZATION_PATH="${RUN_NAME}"

echo SAVE_PATH: $SAVE_PATH

python src/main.py \
    --mode "train" \
    --config_file $CONFIG \
    --run_name $RUN_NAME \
    --save_path $SAVE_PATH \
    --batch_size $BATCH_SIZE \
    --num_folds $NUM_FOLDS \
    --num_gpu $NUM_GPU \
    --gpu $CUDA --seed $SEED \
    --project "DiffDock Tuning" \
    --visualize_n_val_graphs 0 \
    --visualization_path $VISUALIZATION_PATH \
