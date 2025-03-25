# Exact Formula for RX-Differential Probability through Modular Addition for All Rotations

This repository contains supporting code for the [paper](https://tosc.iacr.org/index.php/ToSC/article/view/12087) 

> "Exact Formula for RX-Differential Probability through Modular Addition for All Rotations"

from ToSC 2025 (1), by Alex Biryukov, Baptiste Lambin and Aleksei Udovenko.

See also [slides](./slides.pdf) from the presentation.

## Setup

Requirements:

- python3 (tested on v3.9), pypy3 is recommended for good performance
- gcc/g++ compiler
- Gurobi C++ library installed (and license), for the Alzette part only

## Contents

[verifyTheory](./verifyTheory) contains code related to verifying our theory (Sections 2-3 from the paper):
- verifying main lemmas/theorems (Section 2, Section 2.7);
- verifying constant addition theorem (Section 2.6);
- verifying claims regarding previous formula from TOSC 2016 (Section 3).

[Alzette](./Alzette) contains main code and data related to the MILP-based search of Alzette and its verification (Section 5, Section 5.1):
- MILP search for trails (heuristic + full precise) (Section 5);
- search logs and trails (Appendix D);
- experimental verification of trails/distinguishers with logs.

[Malzette](./Malzette) contains code and data related to the Malzette permutation (Section 6):
- generation of Malzette instances (Section 6.1);
- experimental verification of trails/distinguishers with logs (Section 6.2).

## More Contents

[impossibleDiff](./impossibleDiff) contains code related to experiments with impossible RX-differentials (Appendix G in the paper).

[otherTrails](./otherTrails) contains code and some logs related to RX-trail searches for Speck and Salsa and other primitives.