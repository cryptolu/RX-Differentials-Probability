# Correcting [AL16] formula

## Verify 1
The following command runs original "verification" C script from the paper, with hardcoded transitions showing inconsistency of the formula.
```
$ make verify1
```
Execution result:

```
gcc Verify-thm-1.c -lm -O3 -o verify-thm-1
./verify-thm-1
Transition with incorrect probability (factor x2, example in the paper)
counter0: 16384 14.000000

echo

gcc Verify-thm-1.c -lm -O3 -DMISMATCH -o verify-thm-1-mismatch
./verify-thm-1-mismatch
Transition with incorrect fractional part (not .0 or .415)
counter0: 1523712 20.539159

echo
````

## Verify 2

The following command runs quick verification of statements from the paper (up to n=8 exnhaustively). Pypy3 is recommended to run it fast enough.
```
make verify2-fast
```

The following command runs longer verification of statements from the paper (up to n=12 exnhaustively).Pypy3 is recommended to run it fast enough.
```
make verify2-long
```

Example verification logs are given in [verify2-fast.log](./verify2-fast.log) and [verify2-long.log](./verify2-long.log).