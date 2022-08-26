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


#this will be called one time for each combination of parameter/#thread/graph.
function create_launch {
  PREFIX_OUTCOME=results

  export OMP_NUM_THREADS=${T}

  script_name=___tmp_script_${BINARY}_${T}_G-K${KRON_SIZE}_ED-${KRON_EDGE}_${PARAM_TH}_${PARAM_K}___

  script_folder=${RESULTS_DIR}/${script_name}___exec_dir


  out_fname=${PREFIX_OUTCOME}_${BINARY}___G-K${KRON_SIZE}___ED-${KRON_EDGE}___TH-${PARAM_TH}___K-${PARAM_K}___T-${T}


  is_new="TRUE"
  if [[ -d "${script_folder}" ]]; then
    if [[ -f "${script_folder}/__to_gather__${out_fname}.out" ]]; then
      if [[ -s "${script_folder}/__to_gather__${out_fname}.out" ]]; then 
        echo "experiment ${out_fname} exists already. Skipping"
        is_new="FALSE"
      fi
    fi
  fi

  if [[ "${is_new}" = "TRUE" ]]; then
    mkdir ${script_folder}
    echo ./${BENCH_PATH}/${BINARY} ${ARGS}
    ./${BENCH_PATH}/${BINARY} ${ARGS} > ${script_folder}/${out_fname}.out

    grep 'SSS' ${script_folder}/${out_fname}.out >> ${script_folder}/__to_gather__${out_fname}.out
    grep 'RRR' ${script_folder}/${out_fname}.out >> ${script_folder}/__to_gather__${out_fname}.out
  fi  
}



#slightly different function, needed to sort out some naming inconsistencies and avoid rerunning experiments
function create_launch_baseline {
  PREFIX_OUTCOME=results

  export OMP_NUM_THREADS=${T}



  #check that the experiments does not exists already
  is_new="TRUE"
  script_name=___tmp_script_${BINARY}_${T}_G-K${KRON_SIZE}_ED-${KRON_EDGE}_0____
  script_folder=${RESULTS_DIR}/${script_name}___exec_dir
  out_fname=${PREFIX_OUTCOME}_${BINARY}___G-K${KRON_SIZE}___ED-${KRON_EDGE}___TH-0___K-___T-${T}
  if [[ -d "${script_folder}" ]]; then
    if [[ -f "${script_folder}/__to_gather__${out_fname}.out" ]]; then
      if [[ -s "${script_folder}/__to_gather__${out_fname}.out" ]]; then 
        echo "experiment ${out_fname} exists already. Skipping"
        is_new="FALSE"
      fi
    fi
  fi

  script_name=___tmp_script_${BINARY}_${T}_G-K${KRON_SIZE}_ED-${KRON_EDGE}_0____
  script_folder=${RESULTS_DIR}/${script_name}___exec_dir
  out_fname=${PREFIX_OUTCOME}_${BINARY}___G-K${KRON_SIZE}___ED-${KRON_EDGE}___TH-0___K-___T-${T}


  #create folder, run experiment, grep important results.  
  if [[ "${is_new}" = "TRUE" ]]; then
    mkdir ${script_folder}
    echo ./${BENCH_PATH}/${BINARY} ${ARGS}
    ./${BENCH_PATH}/${BINARY} ${ARGS} > ${script_folder}/${out_fname}.out

    grep 'SSS' ${script_folder}/${out_fname}.out >> ${script_folder}/__to_gather__${out_fname}.out
    grep 'RRR' ${script_folder}/${out_fname}.out >> ${script_folder}/__to_gather__${out_fname}.out
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
        create_launch_baseline
      elif [ "${approx_scheme}" = "colorful" ]; then
        export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -a -n 1 -p 0.5"
        create_launch_baseline
      elif [ "${approx_scheme}" = "doulion" ]; then
        export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -a -n 1 -p 0.8"
        create_launch_baseline
      elif [ "${approx_scheme}" = "1h" ]; then
        export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -t ${OH_TH} -n 1 -y ${CLUSTER}"
        create_launch
      elif [ "${approx_scheme}" = "bf" ]; then
        export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -t ${BF_TH} -b -1 -n 1 -y ${CLUSTER}"
        create_launch          
      fi
  done
done
