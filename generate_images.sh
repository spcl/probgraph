mkdir real_graphs_images
mkdir scaling_images

EXPS=(tc jp-jc jp-cn jp-ov 4c)
for e in ${EXPS[@]}; do
	echo images for $e on real graphs
	python src/scripts/make_comparison_images.py
done

python src/scripts/make_scaling_plots.py
