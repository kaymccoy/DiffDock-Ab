NUM_FOLDS=1  # number of seeds to try, default 5
SEED=0  # initial seed
CUDA=0  # will use GPUs from CUDA to CUDA + NUM_GPU - 1
NUM_GPU=4
BATCH_SIZE=1  # split across all GPUs
NUM_SAMPLES=1

NAME="aadam_notstrict_allLD_val"  #***
RUN_NAME="aadam_notstrict_allLD_val" #***
CONFIG="config/${NAME}.yaml"

SAVE_PATH="ckpts/${RUN_NAME}"
VISUALIZATION_PATH="${RUN_NAME}"

echo SAVE_PATH: $SAVE_PATH

python src/main.py \
    --mode "train" \
    --config_file $CONFIG \
    --run_name $RUN_NAME \
    --save_path $SAVE_PATH \
    --checkpoint_path $SAVE_PATH \
    --batch_size $BATCH_SIZE \
    --num_folds $NUM_FOLDS \
    --num_gpu $NUM_GPU \
    --save_model_every 2 \
    --gpu $CUDA --seed $SEED \
    --project "DiffDock Tuning" \
    --visualize_n_val_graphs 0 \
    --visualization_path $VISUALIZATION_PATH \
