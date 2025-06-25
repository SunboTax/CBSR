#!/bin/bash

cd ../build
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so
# When block_num=1(using pure Topchain), some datasets 
# may need a larger stack size for query processing.
# To make the comparison fair, we set the stack size to 800MB.
# in all experiments.
ulimit -s 819200

mkdir -p ../memory_stats
mkdir -p ../result
mkdir -p ../log

echo "Dataset,BlockNum,Latency(s),MaxMemory(GB)" > ../memory_stats/memory_summary.csv

# define the block number for each dataset
declare -A dataset_block_map
dataset_block_map=(
    ["wikiquote_pl"]=8
    ["stackexchange-stackoverflow"]=8
    ["linux-kernel"]=8
    ["citeulike-ti"]=4
    ["bibsonomy"]=16
    ["twitter"]=8
    ["amazon"]=4
    ["epinions-rating"]=4
    ["lastfm_song"]=2
    ["edit-jawiki"]=4
    ["edit-dewiki"]=4
)

# run the computation for each dataset
for dataset in "${!dataset_block_map[@]}"; do
    block_num=${dataset_block_map[$dataset]}
    for delta in 86400 172800 345600 691200 1382400 2764800;do
        echo "=== Running: ${dataset}, block_num=${block_num}, delta=${delta} ==="

        /usr/bin/time -v ./Computation \
            --inputPath=../dataset/bigraph/${dataset}.txt \
            --queryFilePath=../dataset/query/${dataset}.query \
            --minBlock=${block_num} \
            --algorithm=TopChain \
            --delta=${delta} \
            --outputPath=../log/json-${dataset}_t_${block_num}_delta_${delta}.json \
            > ../log/${dataset}_${block_num}_delta_${delta}.log \
            2> ../memory_stats/${dataset}_${block_num}_delta_${delta}_time.log

        if [ -f "../memory_stats/${dataset}_${block_num}_delta_${delta}_time.log" ]; then
            max_mem_kb=$(grep "Maximum resident set size" "$logfile" | grep -o '[0-9]\+')

            if [ -n "$max_mem_kb" ]; then
                max_mem_gb=$(echo "scale=3; $max_mem_kb / 1024 / 1024" | bc)
            else
                max_mem_gb="N/A"
            fi

            echo "${dataset},${block_num},${delta},${max_mem_gb}" >> ../memory_stats/memory_summary.csv
            echo "Memory Usage: ${max_mem_gb} GB"
        fi

        echo "---"
    done
done
                   