#!/usr/bin/env python3
from sys import argv
import csv


def file_parser(filename, directory):

    string = directory + '/'
    in_txt = csv.reader(open(filename, 'r'), delimiter='\t')
    out_csv = csv.writer(open(string + 'N1zones.csv', 'w+'), delimiter='\t')

    for row in in_txt:
        if row and row[0] != 'Sequence number' and row[2] == 'productive':
            out_csv.writerow([row[1], row[3], row[21]])


def weight(sc, tc, matrix):

    alphabet = "ATGC-"

    row = alphabet.find(sc)
    line = alphabet.find(tc)

    return matrix[line][row]


def dist(i, j, mas, s, t, loc, matrix):
    d = 0
    a = mas[i][j-1] + weight('-', t[j-1], matrix)
    b = mas[i - 1][j] + weight(s[i - 1], '-', matrix)
    c = mas[i - 1][j - 1] + weight(s[i - 1], t[j - 1], matrix)

    max_value = 0
    if a > b:
        max_value = a
        loc[i][j] = 1
    else:
        max_value = b
        loc[i][j] = 2

    if c > max_value:
        max_value = c
        loc[i][j] = 3

    if d > max_value:
        max_value = d
        loc[i][j] = 0

    return max_value


# local alignment
def alignment(s, t):

    matrix = [[50, 0, 0, 0, -20000],
              [0, 50, 0, 0, -20000],
              [0, 0, 50, 0, -20000],
              [0, 0, 0, 50, -20000],
              [-20000, -20000, -20000, -20000, -20000]]

    n = len(s) + 1
    m = len(t) + 1

    value = [[0 for x in range(m)] for y in range(n)]
    loc = [[0 for x in range(m)] for y in range(n)]

    for i in range(1, n):
        loc[i][0] = 2

    for j in range(1, m):
        loc[0][j] = 1

    for i in range(1, n):
        for j in range(1, m):
            value[i][j] = dist(i, j, value, s, t, loc, matrix)

    max_value = 0

    for i in range(0, n):
        for j in range(0, m):
            if max_value < value[i][j]:
                max_value = value[i][j]

    if max_value < 50*(len(t) - 1):
        return -1
    else:
        return len(t)


# find all footprints in N1 for every sequence
# write in footprints_in_N1.txt sequence id, N1, found footprint,
# position of footprint, length of N1 and footprint
# if there are several footprints in one sequence, every
# unique pair (sequence, footprint) will be written
# max - 1 mismatch, 0 gap
def inexact_footprint_search(directory):

    string = directory + '/'

    with open(string + 'N1zones.csv', 'r') as sequences, \
            open(string + 'footprints_in_N1.txt', 'w+') as file_out:

        ids = []
        for seq in sequences:
            t = seq.split('\t')
            seq_id = t[0]
            if seq_id not in ids:
                ids.append(seq_id)
            n1 = t[2][: -1].upper()
            if n1 == '':
                continue
            if n1[-1] == '\n':
                n1 = n1[: -1]

            with open('auxiliary_files/Footprints.txt', 'r') as footprints:
                for footprint in footprints:
                    size = alignment(n1, footprint[:-1])
                    if size > -1:
                        file_out.write(seq_id + '\n')
                        file_out.write(n1 + '\n')
                        file_out.write(footprint)
                        file_out.write('length N1 = ' + str(len(n1)) + '\n')
                        file_out.write('length FT = ' + str(len(footprint) - 1) + '\n')
                        file_out.write('pos = ' + "not defined" + '\n\n')


# find all footprints in N1 for every sequence
# write in footprints_in_N1.txt sequence id, N1, found footprint,
# position of footprint, length of N1 and footprint
# if there are several footprints in one sequence, every
# unique pair (sequence, footprint) will be written
def exact_footprint_search(directory):

    string = directory + '/'
    with open(string + 'N1zones.csv', 'r') as sequences, \
            open(string + 'footprints_in_N1.txt', 'w+') as file_out:

        ids = []
        for seq in sequences:
            t = seq.split('\t')
            seq_id = t[0]
            if seq_id not in ids:
                ids.append(seq_id)
            n1 = t[2][: -1].upper()
            if n1 == '':
                continue

            with open('auxiliary_files/Footprints.txt', 'r') as footprints:
                if n1[-1] == '\n':
                    n1 = n1[: -1]
                for footprint in footprints:
                    index = n1.find(footprint[: -1])
                    displacement = 0
                    while index > -1:
                        file_out.write(seq_id + '\n')
                        file_out.write(n1 + '\n')
                        file_out.write(footprint)
                        displacement += index
                        file_out.write('length N1 = ' + str(len(n1)) + '\n')
                        file_out.write('length FT = ' + str(len(footprint) - 1) + '\n')
                        file_out.write('pos = ' + str(displacement + 1) + '\n\n')
                        displacement += 1
                        index = n1[displacement: len(n1)].find(footprint[: -1])


# create a dict: v-gen - footprints it can make
def vgen_footprints_dict():

    vgen_dict = {}
    with open('auxiliary_files/V-gen-footprints.txt', 'r') as file_in:
        current_genes = []
        is_gen = True
        for line in file_in:
            t = line.split()
            for i in range(len(t)):
                if t[i][0] == 'I':
                    if not is_gen:
                        current_genes.clear()
                        is_gen = True
                    vgen_dict.update({t[i]: []})
                    current_genes.append(t[i])

                elif t[i][0] != 'I':
                    if is_gen:
                        is_gen = False
                    for gen in current_genes:
                        new_list = vgen_dict[gen]
                        new_list.append(t[i])
                        vgen_dict.update({gen: new_list})

    return vgen_dict


# in (in)exact_found_footprints.txt write sequence id, footprint,
# v-gen of sequence, possible parent of found footprint
def write_footprint_parent(directory, option):

    string = directory + '/'
    with open(string + 'footprints_in_N1.txt', 'r') as footprints,\
            open(string + option + '_found_footprints.txt', 'w+') as file_out:

        d = vgen_footprints_dict()

        for i, line in enumerate(footprints):
            if i % 7 == 0:
                file_out.write(line)
                seq_id = line[: -1]

                with open(string + 'N1zones.csv', 'r') as n1_zones:
                    for zone in n1_zones:
                        number, tmp = zone.split('\t')[0], zone.split('\t')[1]
                        if number == seq_id:
                            index1 = tmp.find('I')
                            index2 = tmp.find('*')
                            # print(tmp[index1:])
                            if index2 > -1:
                                v_gen = tmp[index1: index2]
                            else:
                                v_gen = tmp[index1:]
                            file_out.write(v_gen + '\n')
                            break

            elif i % 7 == 2:
                file_out.write(line)
                for key in d.keys():
                    if line[: -1] in d[key]:
                        file_out.write(key)
                        file_out.write(' ')
                file_out.write('\n')

            elif i % 7 != 6:
                file_out.write(line)

            else:
                file_out.write('\n')


def possible_parent(other_v_gen, v_gen, undefined_genes, footprint, position_v_genes):

    defined_genes = (other_v_gen not in undefined_genes) and (v_gen not in undefined_genes)
    not_ORgenes = (other_v_gen.find('OR') == -1) and (v_gen.find('OR') == -1)
    if defined_genes and not_ORgenes:
        return (len(footprint) > 5) and (position_v_genes[v_gen] < position_v_genes[other_v_gen])
    return False


def find_genuine_footprints(directory, option):

    string = directory + '/'
    position_v_genes = {}
    # dict for v-gen place in germline
    place = 0
    with open('auxiliary_files/V-genes.csv', 'r') as v_genes:
        for line in v_genes:
            place += 1
            position_v_genes[line[: -1]] = place

    with open("auxiliary_files/Undefined_genes.txt", 'r') as undef_genes_file:
        undefined_genes = []
        for undef_gen in undef_genes_file:
            undefined_genes.append(undef_gen[:-1])

    only_one = {}

    with open(string + option + '_found_footprints.txt', 'r') as file_in, \
            open(string + option + '_genuine-footprints.txt', 'w+') as file_out:

        v_for_footprints = []
        read_list = []
        for i, line in enumerate(file_in):
            # seq_id, v_gen, n1, footprint, list_fp, lenn1, lenft, pos
            if i % 9 != 4 and i % 9 < 8:
                read_list.append(line)

            elif i % 9 == 4:
                v_for_footprints = line[: -2].split()

            elif i % 9 == 8:
                v_gen = read_list[1][:-1]
                for v_for_footprint in v_for_footprints:
                    index = v_for_footprint.find('*')
                    if index == -1:
                        tmp = v_for_footprint
                    else:
                        tmp = v_for_footprint[: index]

                    if tmp[-1] == '\n':
                        other_v_gen = tmp[: -1]
                    else:
                        other_v_gen = tmp
                    # it is candidate for footprint parent

                    if possible_parent(other_v_gen, v_gen, undefined_genes, read_list[3], position_v_genes):
                        if only_one.get((read_list[0], read_list[3], other_v_gen, read_list[6])) is None:
                            file_out.write(read_list[0])
                            file_out.write(read_list[3])
                            file_out.write(read_list[1])
                            file_out.write(other_v_gen + '\n')
                            file_out.write(read_list[2])
                            file_out.write(read_list[4])
                            file_out.write(read_list[5])
                            file_out.write(read_list[6] + '\n')
                            only_one[(read_list[0], read_list[3], other_v_gen, read_list[6])] = 1
                read_list = []


# Return count of sequences with footprint
def count_sequences(directory, option):
    string = directory + '/'
    with open(string + option + '_genuine-footprints.txt', 'r') as f:
        sequences_with_footprints = {}
        for i, line in enumerate(f):
            if i % 9 == 0:
                number = line[: -3]
                if sequences_with_footprints.get(number) is None:
                    sequences_with_footprints[number] = 1

    return len(sequences_with_footprints.keys())


if __name__ == "__main__":
    if len(argv) > 3:
        file = argv[1]
        directory = argv[2]
        option = argv[3]
        if file.find('csv') == -1:
            file_parser(file, directory)
        if option == "exact":
            exact_footprint_search(directory)
        elif option == "inexact":
            inexact_footprint_search(directory)
        else:
            raise Exception("Type of search can be exact or inexact")
        write_footprint_parent(directory, option)
        find_genuine_footprints(directory, option)
    else:
        raise Exception("Enter file, output directory and search type - exact\inexact")

