export problem=${1}
export BENCH_PATH=${2}
export RESULTS_DIR=${3}

THREADS=( 32 )

if [ "${problem}" = "tc" ]; then
	export BINARIES=(${problem}_base ${problem}_bf ${problem}_1h ${problem}_doulion ${problem}_colorful ${problem}_pgp ${problem}_redex)
else
  export BINARIES=(${problem}_base ${problem}_bf ${problem}_1h ${problem}_pgp ${problem}_redex)
fi

REDEX_THs=(0.5)
BF_THs=(0.5)
OH_THs=(0.01)
KRON_SIZEs=(16 18 20 22)
KRON_EDGEs=(256 1024 4096)
CLUSTER=0.1

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

function create_launch_baseline {
  PREFIX_OUTCOME=results

  export OMP_NUM_THREADS=${T}


  #needed to sort out some naming inconsistencies and avoid rerunning experiments
  is_new="TRUE"
  for PARAM_TH in ${PARAM_THs[@]}; do
    script_name=___tmp_script_${BINARY}_${T}_G-K${KRON_SIZE}_ED-${KRON_EDGE}_${PARAM_TH}_${PARAM_K}___
    script_folder=${RESULTS_DIR}/${script_name}___exec_dir
    out_fname=${PREFIX_OUTCOME}_${BINARY}___G-K${KRON_SIZE}___ED-${KRON_EDGE}___TH-${PARAM_TH}___K-${PARAM_K}___T-${T}
    if [[ -d "${script_folder}" ]]; then
      if [[ -f "${script_folder}/__to_gather__${out_fname}.out" ]]; then
        if [[ -s "${script_folder}/__to_gather__${out_fname}.out" ]]; then 
          echo "experiment ${out_fname} exists already. Skipping"
          is_new="FALSE"
        fi
      fi
    fi
  done

    script_name=___tmp_script_${BINARY}_${T}_G-K${KRON_SIZE}_ED-${KRON_EDGE}__${PARAM_K}___
    script_folder=${RESULTS_DIR}/${script_name}___exec_dir
    out_fname=${PREFIX_OUTCOME}_${BINARY}___G-K${KRON_SIZE}___ED-${KRON_EDGE}___TH-___K-${PARAM_K}___T-${T}
    if [[ -d "${script_folder}" ]]; then
      if [[ -f "${script_folder}/__to_gather__${out_fname}.out" ]]; then
        if [[ -s "${script_folder}/__to_gather__${out_fname}.out" ]]; then 
          echo "experiment ${out_fname} exists already. Skipping"
          is_new="FALSE"
        fi
      fi
    fi


  script_name=___tmp_script_${BINARY}_${T}_G-K${KRON_SIZE}_ED-${KRON_EDGE}_0_${PARAM_K}___
  script_folder=${RESULTS_DIR}/${script_name}___exec_dir
  out_fname=${PREFIX_OUTCOME}_${BINARY}___G-K${KRON_SIZE}___ED-${KRON_EDGE}___TH-0___K-${PARAM_K}___T-${T}


  if [[ "${is_new}" = "TRUE" ]]; then
    mkdir ${script_folder}
    echo ./${BENCH_PATH}/${BINARY} ${ARGS}
    ./${BENCH_PATH}/${BINARY} ${ARGS} > ${script_folder}/${out_fname}.out

    grep 'SSS' ${script_folder}/${out_fname}.out >> ${script_folder}/__to_gather__${out_fname}.out
    grep 'RRR' ${script_folder}/${out_fname}.out >> ${script_folder}/__to_gather__${out_fname}.out
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
          	create_launch_baseline
        elif [ "${approx_scheme}" = "colorful" ]; then
          	export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -a -n 1 -p 0.5"
          	create_launch_baseline
        elif [ "${approx_scheme}" = "doulion" ]; then
	        export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -a -n 1 -p 0.8"
          	create_launch_baseline
	elif [ "${approx_scheme}" = "pgp" ] || [ "${approx_scheme}" = "redex" ]; then
		for PARAM_TH in ${REDEX_THs[@]}; do
			export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -a -n 1 -t ${PARAM_TH} -y ${CLUSTER}"
			create_launch
		done
	elif [ "${approx_scheme}" = "approx1" ] || [ "${approx_scheme}" = "approx2" ]; then
		for PARAM_TH in ${AUTO_THs[@]}; do
			export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -a -n 1 -t ${PARAM_TH} -y ${CLUSTER}"
			create_launch
		done
        elif [ "${approx_scheme}" = "1h" ]; then
	  	for PARAM_TH in ${OH_THs[@]}; do
            		export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -t ${PARAM_TH} -n 1 -y ${CLUSTER}"
            		create_launch
    	  	done
	elif [ "${approx_scheme}" = "bf" ]; then
		for PARAM_TH in ${BF_THs[@]}; do 
            		export ARGS="-s -g ${KRON_SIZE} -k ${KRON_EDGE} -t ${PARAM_TH} -b -1 -n 1 -y ${CLUSTER}"
			create_launch
		done          
        fi
      done
    done
  done
done

