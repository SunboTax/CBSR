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

# a latency of 20 years in seconds which is far longer than the datasets' time span
delta=630720000 

for dataset in wikiquote_pl stackexchange-stackoverflow linux-kernel citeulike-ti bibsonomy twitter amazon epinions-rating lastfm_song edit-jawiki edit-dewiki;
do
    for block_num in  1 2 8 16 32 64;
    do
        echo "=== Running: ${dataset}, block_num=${block_num}, delta=${delta} ==="
        
        # 使用 /usr/bin/time -v 统计详细资源使用情况
        /usr/bin/time -v ./Computation \
            --inputPath=../dataset/bigraph/${dataset}.txt \
            --queryFilePath=../dataset/query/${dataset}.query \
            --minBlock=$block_num \
            --algorithm=TopChain \
            --delta=${delta} \
            --outputPath=../log/json-${dataset}_t_${block_num}_delta_${delta}.json \
            > ../log/${dataset}_${block_num}_delta_${delta}.log \
            2> ../memory_stats/${dataset}_${block_num}_delta_${delta}_time.log
        
        if [ -f "../memory_stats/${dataset}_${block_num}_delta_${delta}_time.log" ]; then
            max_mem_kb=$(grep "Maximum resident set size" ../memory_stats/${dataset}_${block_num}_delta_${delta}_time.log | grep -o '[0-9]*')
            
            if [ -n "$max_mem_kb" ]; then
                max_mem_gb=$(echo "scale=3; $max_mem_kb / 1024 / 1024" | bc)
            else
                max_mem_mb="N/A"
            fi
            
            echo "${dataset},${block_num},${delta},${max_mem_gb}" >> ../memory_stats/memory_summary.csv
            echo "Memory Usage: ${max_mem_gb} GB"
        fi
        
        echo "---"
    done
done                   