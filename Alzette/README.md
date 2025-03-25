# Alzette - Trail Search (MILP) and Verification (Section 5, Section 5.1)

# Part 1: Trail search (MILP) (Section 5)

First, we need to make the command `runAlzette`. Note that Gurobi version is hardcoded in the makefile, and maybe needs to be adapted:

```
... -lgurobi_g++8.5 -lgurobi110
```

```sh
$ make runAlzette
$ ./runAlzette
Usage: ./runAlzette <constIndex> <rotationAmount> <fullOpt-flag> <timeLimitNEQ> <timeLimitFull>
Timelimits in seconds.
```

Now we can run e.g. constant `c0` rotation `k=3` heuristic NEQ optimization (fullOpt flag 0), with 50 seconds limit for heuristic and 100 seconds limits for future run of full (precise) search.

Note that one has to remove the corresponding trail file if it already exists!

```
$ ./runAlzette 0 1 0 60 120
logging to logs/alzette_constb7e15162_rotation1_neqopt_time60.txt
$ cat logs/alzette_constb7e15162_rotation1_neqopt_time60.txt
...
*17386 12412             102      30.0000000    1.50000  95.0%   113   32s
 19171 13512   20.29167   72  177   30.00000    1.50000  95.0%   113   35s
 22691 15379    8.90741   51  209   30.00000    1.50000  95.0%   115   40s
 26951 17246 infeasible   81        30.00000    1.50000  95.0%   114   45s
 30843 18859   27.75000   80  156   30.00000    1.50000  95.0%   112   50s
 34535 20143 infeasible   80        30.00000    1.55013  94.8%   112   60s
in callback 
nbModAdd = 4
nbRound = 5
=========================================================
======== Solution found with Objective : 30 ========
Number of notAllEqual vars : 30
Trail :
Round 0 : 0x2023c000 0x1021f000 
Round 1 : 0xd80413a7 0x37c1f000 
Round 2 : 0xc827fc23 0x3003f802 
Round 3 : 0x40000f81 0x0044004f 
Round 4 : 0xdc23ac27 0x5fc4044f 
0x2023c000 + 0x2043e000 -> 0x0027e000 (2^-6.41504)
0xd80413a7 + 0xf8001be0 -> 0x10040f84 (2^-14)
0xc827fc23 + 0x3003f802 -> 0x9823fc26 (2^-12)
0x40000f81 + 0x44004f00 -> 0x04005f80 (2^-6.41504)
sumLog = -38.8301
modAddTrail = [[0x2023c000,0x2043e000,0x0027e000],
[0xd80413a7,0xf8001be0,0x10040f84],
[0xc827fc23,0x3003f802,0x9823fc26],
[0x40000f81,0x44004f00,0x04005f80]]
=========================================================

Explored 34545 nodes (3864768 simplex iterations) in 60.01 seconds (59.80 work units)
Thread count was 8 (of 8 available processors)

Solution count 10: 30 31 33 ... 42

Time limit reached
Best objective 3.000000000000e+01, best bound 2.000000000000e+00, gap 93.3333%

User-callback calls 82251, time in user-callback 0.01 sec
***********************************************************************
*** Alzette 4 rounds with k = 1 cst = 0xb7e15162 ***
*** Best solution found with objective 30 ***
Solution obtained after time-out
Number of notAllEqual vars : 30
Trail :
x0 : 0x2023c000 y0 : 0x1021f000
x1 : 0xd80413a7 y1 : 0x37c1f000
x2 : 0xc827fc23 y2 : 0x3003f802
x3 : 0x40000f81 y3 : 0x0044004f
x4 : 0xdc23ac27 y4 : 0x5fc4044f
modAdd :
0x2023c000 + 0x2043e000 -> 0x0027e000 (2^-6.41504)
0xd80413a7 + 0xf8001be0 -> 0x10040f84 (2^-14)
0xc827fc23 + 0x3003f802 -> 0x9823fc26 (2^-12)
0x40000f81 + 0x44004f00 -> 0x04005f80 (2^-6.41504)
sumLog = -38.8301
modAddTrail = [[0x2023c000,0x2043e000,0x0027e000],
[0xd80413a7,0xf8001be0,0x10040f84],
[0xc827fc23,0x3003f802,0x9823fc26],
[0x40000f81,0x44004f00,0x04005f80]]
***********************************************************************

```

The found trail will be saved in the trails folder. Now we can run full search by switching the full opt flag:
```
$ ./runAlzette 0 1 1 60 120
logging to logs/alzette_constb7e15162_rotation1_fullopt_time120.txt
$ cat logs/alzette_constb7e15162_rotation1_fullopt_time120.txt
...
H 1458  1026                      36.6601493   24.36461  33.5%   585   91s
  1471  1035   24.75084   14  457   36.66015   24.75084  32.5%   580   95s
  1486  1045   24.85542   25  531   36.66015   24.85542  32.2%   574  100s
  1498  1053   24.93671   28  514   36.66015   24.93671  32.0%   569  105s
  1517  1067   24.94404    8  495   36.66015   24.94404  32.0%   682  110s
  1536  1079   25.48253   20  593   36.66015   25.48253  30.5%   673  115s

Cutting planes:
  Gomory: 30
  Cover: 210
  MIR: 60
  StrongCG: 2
  Flow cover: 632
  Zero half: 138
  RLT: 3

Explored 1555 nodes (1143665 simplex iterations) in 120.01 seconds (127.52 work units)
Thread count was 8 (of 8 available processors)

Solution count 7: 36.6601 36.6601 36.6601 ... 38.8301

Time limit reached
Best objective 3.666014930918e+01, best bound 2.555151060834e+01, gap 30.3017%

User-callback calls 32326, time in user-callback 0.01 sec
***********************************************************************
*** Alzette 4 rounds with k = 1 cst = 0xb7e15162 ***
*** Best solution found with objective 36.6601 ***
Solution obtained after time-out
Number of notAllEqual vars : 31
Trail :
x0 : 0x400df001 y0 : 0x2032f000
x1 : 0xd80403a7 y1 : 0x07c2f000
x2 : 0xc827fc25 y2 : 0x0003f802
x3 : 0x10000f81 y3 : 0x9044004f
x4 : 0x8c23acb6 y4 : 0xcf55544f
modAdd :
0x400df001 + 0x4065e000 -> 0x0027f000 (2^-7.41504)
0xd80403a7 + 0x780003e1 -> 0x10040f82 (2^-11.415)
0xc827fc25 + 0x0003f802 -> 0xc823fc26 (2^-9.41504)
0x10000f81 + 0x44004f90 -> 0x54005f11 (2^-8.41504)
sumLog = -36.6601
modAddTrail = [[0x400df001,0x4065e000,0x0027f000],
[0xd80403a7,0x780003e1,0x10040f82],
[0xc827fc25,0x0003f802,0xc823fc26],
[0x10000f81,0x44004f90,0x54005f11]]
***********************************************************************
Set parameter TimeLimit to value 1e+100
Set parameter SolutionLimit to value 1
Set parameter MIPFocus to value 1
Gurobi Optimizer version 11.0.2 build v11.0.2rc0 (linux64 - "Linux Mint 22.1")

CPU model: 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz, instruction set [SSE2|AVX|AVX2|AVX512]
Thread count: 4 physical cores, 8 logical processors, using up to 8 threads

Optimize a model with 8928 rows, 1592 columns and 32544 nonzeros
Model fingerprint: 0x8d0b7e51
Model has 8 general constraints
Variable types: 0 continuous, 1592 integer (1592 binary)
Coefficient statistics:
  Matrix range     [1e+00, 1e+00]
  Objective range  [1e+00, 1e+00]
  Bounds range     [1e+00, 1e+00]
  RHS range        [1e+00, 4e+00]
Presolve removed 4971 rows and 1024 columns
Presolve time: 0.05s
Presolved: 3957 rows, 568 columns, 13775 nonzeros
Variable types: 0 continuous, 568 integer (568 binary)

Root relaxation: objective 2.952835e+00, 1062 iterations, 0.05 seconds (0.05 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0    2.95284    0  358          -    2.95284      -     -    0s
     0     0    5.31548    0  308          -    5.31548      -     -    0s
     0     0    5.33899    0  313          -    5.33899      -     -    0s
     0     0    5.43656    0  308          -    5.43656      -     -    0s
     0     0    5.43693    0  304          -    5.43693      -     -    0s
     0     0    5.43693    0  305          -    5.43693      -     -    0s
     0     0    5.93361    0  375          -    5.93361      -     -    0s
     0     0    5.96690    0  352          -    5.96690      -     -    0s
     0     0    5.96835    0  332          -    5.96835      -     -    0s
     0     0    5.96875    0  349          -    5.96875      -     -    0s
     0     0    5.96875    0  349          -    5.96875      -     -    0s
     0     0    6.63940    0  333          -    6.63940      -     -    0s
     0     0    6.80596    0  350          -    6.80596      -     -    0s
     0     0    6.81597    0  349          -    6.81597      -     -    0s
     0     0    6.81604    0  347          -    6.81604      -     -    0s
     0     0    7.53149    0  333          -    7.53149      -     -    0s
     0     0    7.58175    0  324          -    7.58175      -     -    0s
     0     0    7.58175    0  314          -    7.58175      -     -    0s
     0     0    7.83046    0  307          -    7.83046      -     -    1s
     0     2    7.83046    0  307          -    7.83046      -     -    4s
    23    27    8.99972    5  337          -    8.41595      -   408    5s
   834   650   11.88889   18  255          -    9.05844      -   311   10s
H 1243   903                      27.0000000    9.12667  66.2%   313   20s

Cutting planes:
  Gomory: 2
  Implied bound: 4
  MIR: 5
  Flow cover: 28
  Zero half: 19
  RLT: 24
  BQP: 1

Explored 1243 nodes (405289 simplex iterations) in 20.39 seconds (26.55 work units)
Thread count was 8 (of 8 available processors)

Solution count 1: 27 

Solution limit reached
Best objective 2.700000000000e+01, best bound 1.000000000000e+01, gap 62.9630%
*******************************
*** Found pair of plaintext ***
x0 = 0xc19254c2 y0 = 0xf614e542
x1 = 0xc3295984 y1 = 0xcc1b3a85
Resulting input dx dy : 0x400df001 0x2032f000
 Expected input dx dy : 0x400df001 0x2032f000
Resulting output dx dy : 0x8c23acb6 0xcf55544f
 Expected output dx dy : 0x8c23acb6 0xcf55544f
Input differentials matches
Output differentials matches
*******************************
```

Resulting trail can be found in the file with the same filename as the log but in the trails folder.

```sh
cat trails/alzette_constb7e15162_rotation1_fullopt_time120.txt
Alzette_b7e15162_k1_nr4
32 1 2 0  # word size, k (rotation), num state words, num key words
31 -36.660150  # NEQ, probability (log2)
5 # num states
# state differences
400df001 2032f000 
d80403a7 07c2f000 
c827fc25 0003f802 
10000f81 9044004f 
8c23acb6 cf55544f 
0
4
# modadd differences (alpha, beta, Delta)
400df001 4065e000 0027f000 
d80403a7 780003e1 10040f82 
c827fc25 0003f802 c823fc26 
10000f81 44004f90 54005f11 
```

The makefile contains commands used to generate full results from the paper (perhaps with slightly adapted timelimits during the revision):

```sh
runA1:
	# c0..c7 k=1..16
	bash -c 'parallel -j 2 ./runAlzette ::: {0..7} ::: {1..16} ::: 0 ::: 3600 ::: 14400' # 1 hour NEQopt
	bash -c 'parallel -j 2 ./runAlzette ::: {0..7} ::: {1..16} ::: 1 ::: 3600 ::: 14400' # 4 hours FullOpt

runA2:
	# c0 k=1..4
	bash -c 'parallel -j 2 ./runAlzette ::: 0 ::: {1..4} ::: 0 ::: 28800 ::: 144000' # 8 hours NEQopt
	bash -c 'parallel -j 2 ./runAlzette ::: 0 ::: {1..4} ::: 1 ::: 28800 ::: 144000' # 40 hours FullOpt
```

The logs are available in the [logs](./logs/) folder. The corresponding trails are available in the [trails](./trails) folder. Finally, [figs](./figs) contain the corresponding graphs of evolution of lower/upper bounds, generated using the [process_log_bounds_evolution.py](./process_log_bounds_evolution.py) script. All this data should correspond to the tables from Appendix D in the paper.


## Part 2: Verification 

First we compile the program.

```sh
make verifyAlzette
```

Now we can verify any trace file with the given amount of data and seeds (hexadecimal values for reproducibility).

```sh
$ ./verifyAlzette trails/alzette_constb7e15162_rotation1_fullopt_time14400.txt 28 123 456

Verification logs from the paper are available in the logs folder with prefix `verify_`, e.g.:

```sh
$ cat logs/verify_alzette_constb7e15162_rotation15_fullopt_time14400.txt
===========================
Normal rounds testing
===========================
Testing Alzette_b7e15162_k15_nr4, using 2^42 data, seeds 0000000000000123 0000000000000456
Trail (state words):
  9fd00c0e dff00e1f 
  0080029c 0ff80000 
  1fd00842 0e900040 
  10100290 10101044 
  1f504c12 56901044 
Found match #1 929b8389:8b5dea97 (... follows trail? 1)
Found match #2 66bdee34:3a59d0a6 (... follows trail? 1)
Found match #3 402d3ddf:d8a3ca8d (... follows trail? 1)
Found match #4 2d75c1a1:4c1faeb0 (... follows trail? 1)
Found match #5 42833a7d:3adfcb24 (... follows trail? 1)
Found match #6 881ebf82:141dd29c (... follows trail? 1)
Found match #7 6756dc81:b01bd88f (... follows trail? 1)
Found match #8 af4b0c3f:9423e6af (... follows trail? 1)
Found match #9 28c3a6b9:3567a6ab (... follows trail? 1)
Found match #10 4f7e9231:2d95e22a (... follows trail? 1)
Found match #11 9387642f:11e352f7 (... follows trail? 1)
Found match #12 31713598:cb1ddd0b (... follows trail? 1)
Found match #13 98e34ac6:06e96dda (... follows trail? 1)
Found match #14 6e9fbb65:aedb8bfc (... follows trail? 1)
Found match #15 325fc44b:b3172210 (... follows trail? 1)
Found match #16 b8ebcedb:0e5db047 (... follows trail? 1)
Found match #17 08219ba7:f623ee78 (... follows trail? 1)
Found match #18 4f89f457:b7e38acd (... follows trail? 1)
Found match #19 9b386d65:8d99d77a (... follows trail? 1)
Found match #20 7b09936c:0d29eab3 (... follows trail? 1)
Found match #21 73739c23:291feebe (... follows trail? 1)
Found match #22 6dc527a5:8d67abb0 (... follows trail? 1)
Found match #23 bdcda4d9:06ddd787 (... follows trail? 1)
Found match #24 8b518cec:0f9d65ef (... follows trail? 1)
Found match #25 24fbf41d:58d3944f (... follows trail? 1)
Found match #26 7de9b202:045da316 (... follows trail? 1)
Found match #27 4d99c8bd:4fd7f264 (... follows trail? 1)
Found match #28 17c71207:cae7e591 (... follows trail? 1)
Found match #29 b8118615:841de987 (... follows trail? 1)
Found match #30 9b8b0ec0:ae69ae90 (... follows trail? 1)
Found match #31 6967aa16:1421d252 (... follows trail? 1)
  failed follow at round 2/4: 00800350 != 008002d0
  failed follow at round 3/4: 0f200802 != 0f400802
Found match #32 49018537:ae65672c (... follows trail? 0)
Found match #33 4f0bd0e5:af1ba4fc (... follows trail? 1)
Found match #34 b7c32e7d:096bcfa6 (... follows trail? 1)
Found match #35 93dbf66e:b05592ac (... follows trail? 1)
  failed follow at round 2/4: 00800150 != 008002d0
  failed follow at round 3/4: 0fa00802 != 0f400802
Found match #36 3bf9a679:ce9de8a4 (... follows trail? 0)
Found match #37 4f944fd2:5461dc94 (... follows trail? 1)
Found match #38 0743df71:39a39937 (... follows trail? 1)
Found match #39 37893ff2:d3e589af (... follows trail? 1)
Found match #40 0787e31e:586b16d0 (... follows trail? 1)
Found match #41 94ecc070:8859aa71 (... follows trail? 1)
Found match #42 70b39b24:08e16237 (... follows trail? 1)
Found match #43 80714bcb:119ff51a (... follows trail? 1)
Found match #44 430b4d63:d31bf8b3 (... follows trail? 1)
Found match #45 1a41bcb1:669f4f77 (... follows trail? 1)
Found match #46 0ebf1894:75d99e4e (... follows trail? 1)
Found match #47 cee9d608:9667e797 (... follows trail? 1)
Found match #48 ba5ea905:0491e153 (... follows trail? 1)
Alzette_b7e15162_k15_nr4: 48/4398046511104 = 2^-36.42 experimental differential probability
Alzette_b7e15162_k15_nr4: 46/48/4398046511104 = 2^-36.48 experimental trail probability
```