"""
For extracting metadata (dipole moments, rotational constants,
and partition function) from .int, .var, and .qpart files, respectively.
"""
from decimal import Decimal


def read_intfile(filein):
    '''Reads in .int file and returns diple moments'''
    file = filein.read().splitlines()
    mu_a, mu_b, mu_c = None, None, None
    for line in file[2:]:
        split_line = line.split()
        if split_line[0] == b'001':
            mu_a = split_line[1].decode()
        elif split_line[0] == b'002':
            mu_b = split_line[1].decode()
        elif split_line[0] == b'003':
            mu_c = split_line[1].decode()
    return mu_a, mu_b, mu_c


def read_varfile(filein):
    '''Reads in .var file and returns rotational constants'''
    file = filein.read().splitlines()
    # a_const is third line second column
    a_const = Decimal(file[3].split()[1].decode())
    # b_const is fourth line second column
    b_const = Decimal(file[4].split()[1].decode())
    # c_const is fifth line second column
    c_const = Decimal(file[5].split()[1].decode())
    return a_const, b_const, c_const


def read_qpartfile(filein):
    '''Reads in .qpart file and returns partition function'''
    file = filein.read().splitlines()
    partition_dict = {}
    for line in file[1:]:
        split_line = line.split()
        partition_dict[split_line[0].decode()] = split_line[1].decode()
    if '300.000' not in partition_dict:
        raise ValueError('Partition function does not contain 300.000 K')
    return partition_dict
