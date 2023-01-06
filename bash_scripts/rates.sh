#!/bin/bash

mkdir 0
cd 0
cp ../coord.xml .
cp -r ../ConfGen.py .
cp -r ../SOMA .       
sed -i "/<conversion_rate>/,/<\/conversion_rate>/ s/rate_A/0.0005/g" coord.xml 
python ConfGen.py -i coord.xml
nohup mpiexec -n 1 ./SOMA -c coord.h5 -a coord_ana.h5  -t 200000 -o 0 &
cd ..

for rates in {1..9}
do
	rate=$(echo "scale=5; $rates/1000" | bc -l)
	num_gpu=$((rates%2))
	mkdir $rates
	cd $rates
	cp ../coord.xml .
	cp -r ../ConfGen.py .
	cp -r ../SOMA .       
	sed -i "/<conversion_rate>/,/<\/conversion_rate>/ s/rate_A/$rate/g" coord.xml 
	python ConfGen.py -i coord.xml
	nohup mpiexec -n 1 ./SOMA -c coord.h5 -a coord_ana.h5  -t 200000 -o ${num_gpu} &
	cd ..
done
