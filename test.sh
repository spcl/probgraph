PROBLEMS=(tc)
KRON_SIZE=14
KRON_EDGES=64
OPTIONS="-g ${KRON_SIZE} -k ${KRON_EDGES} -y 0.1"

for PROB in ${PROBLEMS[@]}; do	
	printf "\n***************************************************************\n"
	printf "*************RUNNING ${PROB} with BASELINE **********************\n"
	printf "*****************************************************************\n"
	./src/${PROB}_base ${OPTIONS}
	
	printf "\n**************************************************************\n"
	printf "*************RUNNING ${PROB} with ONE HASH**********************\n"
	printf "****************************************************************\n"
	./src/${PROB}_1h ${OPTIONS} -t 0.01
	
	printf "\n**************************************************************\n"
	printf "*************RUNNING ${PROB} with BLOOM FILTERS ****************\n"
	printf "****************************************************************\n"
	./src/${PROB}_bf ${OPTIONS} -t 0.5
done
