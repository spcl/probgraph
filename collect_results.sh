PROBLEMS=(tc jp-jc jp-cn jp-ov)

for PROB in ${PROBLEMS[@]}; do
	echo collecting kronecker results for ${PROB}
	python src/scripts/collect_into_csv.py --input-dir kronecker_graphs_results/${PROB}/ --output-name kronecker_graphs_results/${PROB}-kron.csv
done

for PROB in ${PROBLEMS[@]}; do
	echo collecting real graphs results for ${PROB}
	python src/scripts/collect_into_csv.py --input-dir real_graphs_results/${PROB}/ --output-name real_graphs_results/${PROB}-real.csv
done

for PROB in ${PROBLEMS[@]}; do
	echo collecting strong results for ${PROB}
	python src/scripts/collect_into_csv.py --input-dir scaling_results/strong_results/${PROB}/ --output-name scaling_results/strong_scaling_${PROB}.csv
done

for PROB in ${PROBLEMS[@]}; do
	echo collecting weak results for ${PROB}
	python src/scripts/collect_into_csv.py --input-dir scaling_results/weak_results/${PROB}/ --output-name scaling_results/weak_scaling_${PROB}.csv
done
