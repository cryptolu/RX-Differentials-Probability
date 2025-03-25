# Theory verification

## Part 1: Verify theory (lemmas/theorems) / Section 2

[verifyTheory.cpp](./verifyTheory.cpp) contains code for verifying theory from the paper (described in Section 2.7 of the paper). It can be run with

```
$ make verifyTheory && ./verifyTheory
```

Log of such execution is recorded in [verifyTheory.log](./verifyTheory.log).


## Part 2: Verify add-constant formula  / Section 2.6

The code for verifying the add-constant formula is given in [verify_add_const.py](./verify_add_const.py), mentioned in the last paragraph of Section 2.6:

```
$ pypy3 verify_add_const.py
```

Log of such execution is recorded in [verify_add_const.log](./verify_add_const.log).


## Part 3: Verify claims about the previous theorem (TOSC 2016) / Section 3

### Verify 1
The following command runs original "verification" C script from the paper (described in Section 3, last paragraph), with hardcoded transitions, showing inconsistency of the formula.
```
$ make --silent verify1
```
Execution result:

```
Transition with incorrect probability (factor x2, example in the paper)
counter0: 16384

16384 = 2^14 vs 2^15 (2^32 * 2^-17 prob.) claimed by previous theorem

Transition with incorrect fractional part (not .0 or .415)
counter0: 1523712

1523712 = 2^20.53915881110803 impossible by previous theorem
````

### Verify 2

The following command runs quick verification of statements from the paper (up to $n=8$ exnhaustively). Pypy3 is recommended (and hardcoded) to run it fast enough.
```
make verify2-fast
# or 
python3 verify_theory.py 8
```

The following command runs longer verification of statements from the paper (up to $n=12$ exnhaustively).Pypy3 is recommended to run it fast enough.
```
make verify2-long
# or 
python3 verify_theory.py 12
```

Example verification logs are given in [AL16_verify2-fast.log](./AL16_verify2-fast.log) and [AL16_verify2-long.log](./AL16_verify2-long.log).