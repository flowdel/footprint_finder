#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from sys import argv
import csv
import random


def create_families_dict(partis_output):
    families = {}
    a = 0
    file1 = csv.reader(open(partis_output, 'r'), delimiter=',')
    for j in file1:
        a += 1
        families[a] = j[0].split(':')
    return families


def create_sequences_dict(sequences):
    seq_dict = {}
    current_id = ''
    with open(sequences, 'r') as sequences:
        # file with sequences исходная fasta
        for seq in sequences:
            if seq[0] == '>':
                current_id = seq[1:-1]
                seq_dict[current_id] = next(sequences)
            else:
                seq_dict[current_id] += seq

    return seq_dict


def choose_from_family(seq_dict):
    list_seq = []
    with open('families.txt', 'r') as f, open('chosen_sequences.fasta', 'w+') as w:
        begin = False

        for line in f:
            if line[0:6] == 'family' and begin:
                seq_id = random.choice(list_seq)
                w.write('>'+seq_id + '\n')
                w.write(seq_dict[seq_id]+'\n')
                list_seq = []
            elif line[0:6] == 'family' and next(f) != 'family':
                begin = True
            elif line[0] == '>':
                list_seq.append(line[1:-1])


def divide_by_family(sequences, partis_output):
    families = create_families_dict(partis_output)
    seq_dict = create_sequences_dict(sequences)
    with open('families.txt', 'w+') as w:
        # file with sequences separated into families
        for family_number in families.keys():
            w.write('family' + '\n')
            for k in families[family_number]:
                if k in seq_dict:
                    w.write('>' + k + '\n' + seq_dict[k] + '\n')

    choose_from_family(seq_dict)


if __name__ == "__main__":

    if len(argv) > 2:
        divide_by_family(argv[1], argv[2])
    else:
        raise Exception("Enter file with sequences and corresponding Partis output")
