import time, os
from cascada.differential.difference import XorDiff, RXDiff
from cascada.smt.chsearch import ChModelAssertType, PrintingMode, round_based_ch_search, INCREMENT_NUM_ROUNDS
from cascada.primitives import speck

from cascada.bitvector.core import Constant
from cascada.bitvector.operation import RotateLeft, RotateRight

from cascada.bitvector.ssa import RoundBasedFunction

from multiprocessing import Process
from time import sleep

def print_time():
    # 1 hour every minute
    for _ in range(60):
        os.system("date")
        sleep(60)

    # 4 hours every 5 minutes
    for _ in range(4 * 60 // 5):
        os.system("date")
        sleep(5 * 60)

    # then every hour
    while True:
        os.system("date")
        sleep(3600)

Printer = Process(target = print_time)
Printer.start()


class Alzette(RoundBasedFunction):
    num_rounds = 4
    input_widths = [32, 32]
    output_widths = [32, 32]

    ROT = [
        [31, 24],
        [17, 17],
        [0, 31],
        [24,16],
    ]
    CONST = 0xb7e15162

    @classmethod
    def set_num_rounds(cls, new_num_rounds):
        cls.num_rounds = new_num_rounds

    @classmethod
    def eval(cls, x, y):
        for i in range(cls.num_rounds):
            x += RotateRight(y, cls.ROT[i % 4][0])
            y ^= RotateRight(x, cls.ROT[i % 4][1])
            x ^= cls.CONST

            cls.add_round_outputs(x, y)

        return x, y

CONSTS = [
    0xb7e15162,
    0xbf715880,
    0x38b4da56,
    0x324e7738,
    0xbb1185eb,
    0x4f7c7b57,
    0xcfbfa1c8,
    0xc2b3293d,
]
from cascada.differential.opmodel import RXModelBvAdd

#RXModelBvAdd.precision = 3  # default
#RXModelBvAdd.precision = 5

for const in CONSTS:
    Alzette.CONST = const
    assert_type = ChModelAssertType.ValidityAndWeight
    iterator = round_based_ch_search(Alzette, 2, 4, RXDiff, assert_type, "btor",
        extra_chfinder_args={"exclude_zero_input_prop": False, "printing_mode": PrintingMode.Debug},
        extra_findnextchweight_args={"initial_weight": 0})

    print()
    print("=====================")
    print("Const 0x%08x" % const)
    print("=====================")
    t0 = time.time()
    while True:
        try:
            num_rounds, ch = next(iterator)
        except StopIteration:
            break
        iterator.send(INCREMENT_NUM_ROUNDS)
        print(num_rounds, ":", ch.srepr(), "time", time.time() - t0)

    print("Finished const 0x%08x" % const)

print("All work done!")
Printer.terminate()
