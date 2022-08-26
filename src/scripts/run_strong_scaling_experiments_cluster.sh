export problem=${1}
export BENCH_PATH=${2}
export RESULTS_FILE=${4}

THREADS=( 1 2 4 8 16 32 )

if [ "${problem}" = "tc" ]; then
  export BINARIES=(${problem}_base ${problem}_bf ${problem}_1h ${problem}_doulion ${problem}_colorful)
else
  export BINARIES=(${problem}_base ${problem}_bf ${problem}_1h)
fi

#PARAM_THs=(0.01 0.03 0.1 0.3 0.5 1.0 2.0)
OH_TH=0.001
BF_TH=0.1


KRON_SIZEs=(20)
KRON_EDGEs=(256)
CLUSTER=0.1 #for clustering experiments


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

for T in ${THREADS[@]}; do
  export OMP_NUM_THREADS=${T}
  
  for KRON_SIZE in ${KRON_SIZEs[@]}; do
    for KRON_EDGE in ${KRON_EDGEs[@]}; do
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
  done
done


