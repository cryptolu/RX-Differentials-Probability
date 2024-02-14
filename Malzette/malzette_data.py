from binteger import Bin

consts_Malzette1 = [
    (0x00000000, 0x4e381c1c),
    (0x2aaaaaaa, 0x36dbe492),
    (0x7fffffff, 0x1236db6c),
    (0x55555555, 0x0763638e),
    (0x2aaaaaaa, 0x1b6d4949),
    (0x55555555, 0x638ef1c7),
    (0x00000000, 0x47638e39),
    (0x2aaaaaaa, 0x5236b6db),
    (0x55555555, 0x4e381c1c),
    (0x7fffffff, 0x638eb1c7),
    (0x7fffffff, 0x47638e39),
    (0x3f2bb31e, 0xb6c004cc),
]

consts_Malzette2 = [
    (0x1c71c924, 0x249cad47),
    (0x49249c71, 0x1249871c),
    (0x6db6c71c, 0x5b127ffe),
    (0x38e39249, 0x152ad249),
    (0x638e36db, 0x649cad55),
    (0x1c71c7ff, 0x471c9492),
    (0x36db6d55, 0x63f1c71d),
    (0x471c7249, 0x36a4ff1c),
    (0x4924938e, 0x5b6c8e47),
    (0x2aab6db6, 0x71c736db),
    (0x6db638e3, 0x55b9c71d),
    (0xfb3d2330, 0xb6da4b61),
]

differential_Malzette1 = [
    (0x7ffffff8, 0x3ffffffe),
    (0xb989d417, 0x0347a2a9),
]
differential_Malzette2 = [
    (0x7fff0003, 0xffffc001),
    (0xdd2bc54b, 0x078b106c),
]

# default Alzette rotations
# repeated 3 times (12 rounds)
# rotations are to the right
rotations = [
    (31, 24),
    (17, 17),
    (0, 31),
    (24, 16),
] * 3

# more inputs can be found easily
# see verify.cpp and change seed
input_Malzette1 = 0x618e890e, 0x37326dd5
input_Malzette2 = 0x3759c889, 0x0bd01faa


def Malzette(l, r, consts):
    for i in range(12):
        l += r.ror(rotations[i][0])
        r ^= l.ror(rotations[i][1])
        l ^= consts[i][0]
        r ^= consts[i][1]
    return l, r


if __name__ == '__main__':
    l, r = [Bin(v, 32) for v in input_Malzette1]
    ll = l.rol(3) ^ differential_Malzette1[0][0]
    rr = r.rol(3) ^ differential_Malzette1[0][1]

    l, r = Malzette(l, r, consts_Malzette1)
    ll, rr = Malzette(ll, rr, consts_Malzette1)

    assert ll == l.rol(3) ^ differential_Malzette1[1][0]
    assert rr == r.rol(3) ^ differential_Malzette1[1][1]
    print("Malzette1 test ok")

    l, r = [Bin(v, 32) for v in input_Malzette2]
    ll = l.rol(3) ^ differential_Malzette2[0][0]
    rr = r.rol(3) ^ differential_Malzette2[0][1]

    l, r = Malzette(l, r, consts_Malzette2)
    ll, rr = Malzette(ll, rr, consts_Malzette2)

    assert ll == l.rol(3) ^ differential_Malzette2[1][0]
    assert rr == r.rol(3) ^ differential_Malzette2[1][1]
    print("Malzette2 test ok")
