#!/bin/bash
for t in 11 20 30 40 50 60 70 80 90 100
do
	echo "c $t"
	mkdir "LengthMonitoring$t"
	cp Volator "LengthMonitoring$t/Volator"
	cp InputVolator.txt "LengthMonitoring$t/InputVolator.txt"

	perl -pi -e "s/MonitoringDuration=11/MonitoringDuration=$t/g" "LengthMonitoring$t/InputVolator.txt"

	cd "LengthMonitoring$t"
	./Volator
	cd ..
done
