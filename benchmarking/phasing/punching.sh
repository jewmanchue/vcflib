#!/bin/bash
for i in $(seq 1 20);
do
	python knockout_GT.py s612.1.vcf 0.3 > knocked_out/s612_$i\_3.vcf &
done