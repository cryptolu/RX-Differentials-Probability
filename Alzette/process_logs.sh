#!/bin/bash -e
cd logs
for f in *time*.txt; do
	python ../process_log_bounds_evolution.py "$f" "../figs/$f.png";
done
for f in *time*.txt; do
	python ../process_log_bounds_evolution.py "$f" "../figs/$f.pdf";
done
