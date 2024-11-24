import sys
from process_logs import LogFile

if __name__ == '__main__':
    log_file = sys.argv[1]

    try:
        plot_file = sys.argv[2]
    except IndexError:
        plot_file = None

    try:
        sec_limit = int(sys.argv[3])
    except IndexError:
        sec_limit = 2**128

    print("Processing", log_file)
    F = LogFile(log_file=log_file, sec_limit=sec_limit)
    F.plot_evolution(plot_file=plot_file)


'''
import re
import sys

times = []
lower = []
upper = []

timesX = []
lowerX = []

timesY = []
upperY = []

prev = -1

log_file = sys.argv[1]

try:
    plot_file = sys.argv[2]
except IndexError:
    plot_file = None

try:
    sec_limit = int(sys.argv[3])
except IndexError:
    sec_limit = 2**128

print("Processing", log_file)
model_count = 0
for line in open(log_file):
    if line.startswith("Gurobi Optimizer"):
        model_count += 1
        if model_count >= 2:
            break
        continue
    match = re.search(r"\s+[\d.\-]+\s+[\d.\-]+s$", line)
    # Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time
    if not match:
        if line.strip().endswith("s"):
            if line.strip().endswith("threads"): continue
            if line.strip().endswith("nonzeros"): continue
            if line.strip().endswith("constraints"): continue
            if line.strip().endswith("columns"): continue
            if line.startswith("Presolve"): continue
            # make sure we don't miss useful lines
            print("skip", repr(line))
        continue
    # continue
    *_, incumb, bound, gap, spd, time = line.strip().split()
    if incumb == "-":
        incumb = 0.0
    else:
        incumb = float(incumb)
    bound = float(bound)
    time = int(time.rstrip("s"))
    if incumb > 40:
        continue
    # print(time, bound, "...", incumb, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" if bound_improved else "")

    if "neqopt" in log_file:
        bound_improved = (int(prev) != int(bound))
    else:
        bound_improved = not lower or (bound - max(lower)) > 1e-9

    if bound_improved:
        timesX.append(time)
        lowerX.append(bound)

    if not upper or -(incumb - min(upper)) > 1e-9:
        timesY.append(time)
        upperY.append(incumb)

    times.append(time)
    lower.append(bound)
    upper.append(incumb)
    # print(line.strip())
    # print(time, bound, incumb)
    # print()

    prev = bound
    if time > sec_limit:
        break


import matplotlib.pyplot as plt
import numpy as np



fig, ax = plt.subplots()

ax.plot(timesX, lowerX, "x", color="blue")
ax.plot(times, lower, "-", color="blue")

ax.plot(timesY, upperY, "x", color="red")
ax.plot(times, upper, "-", color="red")

ax.set(xlabel='time (s)', ylabel='-log(prob)', title='')
ax.grid()

if plot_file:
    fig.savefig(plot_file)
else:
    fig.savefig("test.png")
    plt.show()
'''
