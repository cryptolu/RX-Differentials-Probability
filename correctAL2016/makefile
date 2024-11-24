all: verify1 verify2


verify1: Verify-thm-1.c
	gcc Verify-thm-1.c -lm -O3 -o verify-thm-1
	./verify-thm-1
	echo
	gcc Verify-thm-1.c -lm -O3 -DMISMATCH -o verify-thm-1-mismatch
	./verify-thm-1-mismatch
	echo

verify2-fast: verify_theory.py
	echo Verifying up to n=8
	# note: can replace pypy3 with python3 but it is 30x slower
	pypy3 verify_theory.py 8

verify2-long: verify_theory.py
	echo Verifying up to n=12
	# note: can replace pypy3 with python3 but it is 30x slower
	pypy3 verify_theory.py 12


