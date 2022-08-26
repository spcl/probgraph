PROBLEMS=(tc jp-jc jp-cn jp-ov)

#real graphs experiments
for PROB in ${PROBLEMS[@]}; do
	mkdir real_graphs_results/${PROB}
	./src/scripts/run_generic_problem_real_graphs.sh ${PROB} ./src/ ./src/graphs/ ./real_graphs_results/${PROB}/
done

#strong scaling experiments
mkdir scaling_results/strong_results/
for PROB in ${PROBLEMS[@]}; do
	mkdir scaling_results/strong_results/${PROB}
	./src/scripts/run_strong_scaling_experiments.sh ${PROB} ./src/ ./scaling_results/strong_results/${PROB}/
done

#weak scaling experiments
mkdir scaling_results/weak_results/
for PROB in ${PROBLEMS[@]}; do
	mkdir scaling_results/weak_results/${PROB}/
	./src/scripts/run_weak_scaling_experiments.sh ${PROB} ./src/ ./scaling_results/weak_results/${PROB}/
done




