import re
import sys

import matplotlib.pyplot as plt
import numpy as np


class LogFile:
    def __init__(self, log_file, sec_limit=2**128):
        self.log_file = log_file
        self.times = []
        self.lower = []
        self.upper = []

        self.timesX = []
        self.lowerX = []

        self.timesY = []
        self.upperY = []

        self.has_no_match = None

        prev = -1

        for line in open(self.log_file):
            if "No matching pair" in line:
                assert self.has_no_match is None
                self.has_no_match = True
            elif "Found pair of" in line:
                assert self.has_no_match is None
                self.has_no_match = False

        # print("Processing", self.log_file)
        model_count = 0
        for line in open(self.log_file):
            if line.startswith("Gurobi Optimizer"):
                model_count += 1
                if model_count >= 2:
                    # print("Stop on 2nd Gurobi opt.")
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
            if incumb > 48:
                continue
            # print(time, bound, "...", incumb, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" if bound_improved else "")

            if "neqopt" in log_file:
                bound_improved = (int(prev) != int(bound))
            else:
                bound_improved = not self.lower or (bound - max(self.lower)) > 1e-9

            if bound_improved:
                self.timesX.append(time)
                self.lowerX.append(bound)

            if not self.upper or -(incumb - min(self.upper)) > 1e-9:
                self.timesY.append(time)
                self.upperY.append(incumb)

            self.times.append(time)
            self.lower.append(bound)
            self.upper.append(incumb)

            prev = bound
            if time > sec_limit:
                break

    def plot_evolution(self, plot_file=None):
        fig, ax = plt.subplots()

        ax.plot(self.timesX, self.lowerX, "x", color="blue")
        ax.plot(self.times, self.lower, "-", color="blue")

        ax.plot(self.timesY, self.upperY, "x", color="red")
        ax.plot(self.times, self.upper, "-", color="red")

        ax.set(xlabel='time (s)', ylabel='-log(prob)', title='')
        ax.grid()

        if plot_file:
            fig.savefig(plot_file)
        else:
            fig.savefig("test.png")
            plt.show()

