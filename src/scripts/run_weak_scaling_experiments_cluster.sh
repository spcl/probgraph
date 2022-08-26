export problem=${1} #one of tc, 4c, jp-jc, jp-cn, jp-ov
export BENCH_PATH=${2} #the directory where the executables (e.g. tc_base, tc_1h) are
export RESULTS_DIR=${3} #output directory. Each experiment will be stored in a different folder inside RESULTS_DIR

THREADS=( 1 2 4 8 16 32 )

#TC has two more baselines: doulion and colorful.
if [ "${problem}" = "tc" ]; then
  export BINARIES=(${problem}_base ${problem}_bf ${problem}_1h ${problem}_doulion ${problem}_colorful)
else
  export BINARIES=(${problem}_base ${problem}_bf ${problem}_1h)
fi

#weak scaling done as follow:
#for T threads, we generate a kroneker graph with 2^20 nodes and (T^2)*EDGE_MULTIPLIER edges
EDGE_MULTIPLIER=4
KRON_SIZE=20

#good parameters for BF and 1H
OH_TH=0.001 
BF_TH=0.1

CLUSTER=0.1 #only used for clustering experiments (jp-jc, jp-cn, jp-ov)


function create_launch {
  PREFIX_OUTCOME=results

  script_body="#!/bin/bash
#SBATCH --job-name=PB_${BINARY}___${T}_${GRAPH_NAME}_${PARAM_TH}_${PARAM_K}
#SBATCH --time=00:10:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${T}

export OMP_NUM_THREADS=${T}

export out_fname=${PREFIX_OUTCOME}_${BINARY}___G-${GRAPH_NAME}___TH-${PARAM_TH}___K-${PARAM_K}___T-${T}

${BENCH_PATH}/${BINARY} ${ARGS} > \${out_fname}.out 2>&1

sleep 1s

grep 'SSS' \${out_fname}.out >> __to_gather__\${out_fname}.out
grep 'RRR' \${out_fname}.out >> __to_gather__\${out_fname}.out
"
      
  script_name=___tmp_script_${BINARY}_${T}_${GRAPH_NAME}_${PARAM_TH}_${PARAM_K}___

  script_folder=${script_name}___exec_dir

  if [[ -d ${script_folder} ]]; then
    echo "experiment exists already"
  else

    mkdir ${script_folder}

    echo "${script_body}" > ${script_folder}/${script_name}.sbatch

    cd ${script_folder}

    sbatch ${script_name}.sbatch
    
    cd ..
  fi
}


#main experiment loop. 
for T in ${THREADS[@]}; do
  export OMP_NUM_THREADS=${T}
  KRON_EDGE=$((${T} * ${T} * ${EDGE_MULTIPLIER}))
  echo "======= Processing graph Kron ${KRON_SIZE} avg. deg ${KRON_EDGE}"

  for BINARY in ${BINARIES[@]}; do
      approx_scheme=$(echo ${BINARY} | cut -d'_' -f 2)
      
      if [ "${approx_scheme}" = "base" ]; then
        export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -a -n 1 -y ${CLUSTER}"
        create_launch
      elif [ "${approx_scheme}" = "colorful" ]; then
        export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -a -n 1 -p 0.5"
        create_launch
      elif [ "${approx_scheme}" = "doulion" ]; then
        export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -a -n 1 -p 0.8"
        create_launch
      elif [ "${approx_scheme}" = "1h" ]; then
        export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -t ${OH_TH} -n 1 -y ${CLUSTER}"
        create_launch
      elif [ "${approx_scheme}" = "bf" ]; then
        export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -t ${BF_TH} -b -1 -n 1 -y ${CLUSTER}"
        create_launch          
      fi
  done
done


