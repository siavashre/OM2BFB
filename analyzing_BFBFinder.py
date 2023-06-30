from collections import defaultdict
import argparse
import pickle
import numpy as np
from math import*
import os
from numpy.linalg import norm

parser = argparse.ArgumentParser()
parser.add_argument("-bfbf", "--bfbfinder", help="BFBfinder output", required=True)
parser.add_argument("-a", "--answer", help="pipeline answer output", required=True)
parser.add_argument("-o", "--output", help="pipeline answer output", required=True)
parser.add_argument("-arm", "--arm", help="Specifying which chromosome arm it is", required=True)
parser.add_argument("-seg", "--segments", help="Segments dictionary", required=True)
parser.add_argument("-m", "--m", help="Number of time divided by two", required=True)
parser.add_argument("-p", "--p", help="Penalty", required=True)
parser.add_argument("-norm", "--norm", help="Power of normalizer", required=True)
parser.add_argument("-segscore", "--segscore", help="Segmentation Score", required=True)
parser.add_argument("-foldback", "--foldback", help="Foldback", required=True)
parser.add_argument("-name", "--name", help="Foldback", required=True)
parser.add_argument("-centro", "--centro", help="Foldback", required=False)
parser.add_argument("-rcmap", "--rcmap", help="Foldback", required=False)
args = parser.parse_args()


def parse_segments(file):
    d = {}
    with open(file, 'r') as f:
        for line in f:
            if not line.startswith('#Segment'):
                line = line.strip().split('\t')
                segment = line[0]
                right_foldback = float(line[6])
                left_foldback = float(line[8])
                cn = float(line[4])
                d[segment] = {'cn': cn, 'left': left_foldback, 'right': right_foldback}
    return d


def euclidean_distance(x, y):
    return sqrt(sum(pow(a - b, 2) for a, b in zip(x, y)))

def parse_foldbacks(file):
    d = {}
    with open(file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                segment = line[2]
                if segment not in d.keys():
                    d[segment] = {'right': [], 'left': []}
                start = float(line[4])
                end = float(line[5])
                q_id = int(line[0])
                dir = line[1]
                if dir == 'right':
                    d[segment]['right'].append([q_id, start, end])
                else:
                    d[segment]['left'].append([q_id, start, end])
    return d


def reverse(a):
    if a == 'right':
        return 'left'
    return 'right'


def calculate_cn(bfb):
    d = defaultdict(int)
    for c in bfb:
        d[c] += 1
    current = 'right'
    if args.arm == 'p':
        current = 'left'
    d_left_right = defaultdict(lambda: defaultdict(int))
    for i in range(len(bfb) - 1):
        if bfb[i] == bfb[i + 1]:
            d_left_right[bfb[i]][current] += 1
            current = reverse(current)
    return d, d_left_right


def parse_BFBFinderOutput(file):
    ans = []
    with open(file, 'r') as f:
        for line in f:
            if not line.startswith('Total'):
                ans.append(line.strip())
    return ans


def miroring_bfb(bfb, m):
    while m > 0:
        bfb = bfb + bfb[::-1]
        m = m - 1
    return bfb


def calculate_score(segment_data, answer):
    cn, foldbacknumber = calculate_cn(answer)
    merged_segments = 0
    merged_id = []
    for i in range(len(list(segment_data.keys())) - 1):
        seg = list(segment_data.keys())[i]
        next_seg = list(segment_data.keys())[i + 1]
        if cn[seg] == cn[next_seg]:
            if foldbacknumber[seg]['right'] == 0:  # and segment_data[seg]['right'] == 0:
                if foldbacknumber[next_seg]['left'] == 0:  # and segment_data[next_seg]['left'] == 0:
                    merged_segments += 1
                    if len(merged_id) == 0:
                        merged_id.append([i, i + 1])
                    else:
                        if i - merged_id[-1][-1] == 0:
                            merged_id[-1].append(i + 1)
                        else:
                            merged_id.append([i, i + 1])
    cn_score = 0
    obs_foldback_dist = []
    expect_foldback_dist = []
    segment_number = 0
    for i, segment in enumerate(segment_data.keys()):
        segment_number += 1
        find = False
        for cluster in merged_id:
            if i in cluster:
                find = True
                cn_score += abs(segment_data[segment]['cn'] - cn[segment]) / (
                        segment_data[segment]['cn'] * len(cluster))
        if find == False:
            cn_score += abs(segment_data[segment]['cn'] - cn[segment]) / segment_data[segment]['cn']
        expect_foldback_dist.append(segment_data[segment]['left'])
        expect_foldback_dist.append(segment_data[segment]['right'])
        obs_foldback_dist.append(foldbacknumber[segment]['left'])
        obs_foldback_dist.append(foldbacknumber[segment]['right'])
    sum_exp = sum(expect_foldback_dist)
    sum_obs = sum(obs_foldback_dist)
    new_obs_foldback_dist = []
    if sum_obs != 0:
        for i in obs_foldback_dist:
            new_obs_foldback_dist.append(i / sum_obs)
    else:
        new_obs_foldback_dist = obs_foldback_dist
    expect_foldback_dist[:] = [x / sum_exp for x in expect_foldback_dist]
    foldback_score = 0
    counter = 0
    penalty = float(args.p)
    for i in range(len(new_obs_foldback_dist)):
        if new_obs_foldback_dist[i] != 0 and expect_foldback_dist[i] == 0:
            foldback_score += penalty*2
            counter +=1
            if 1 < (segment_number - merged_segments) < 5 and (0 < i <len(new_obs_foldback_dist)-1):
                counter +=1
        #### sohbat konim aya bayad bashse ya nabashe??????????
        elif new_obs_foldback_dist[i] == 0 and expect_foldback_dist[i] != 0:
            foldback_score += abs(new_obs_foldback_dist[i] - expect_foldback_dist[i]) / sum_exp
        elif expect_foldback_dist[i] != 0:
            foldback_score += min(penalty*2,
                                  abs(new_obs_foldback_dist[i] - expect_foldback_dist[i]) / expect_foldback_dist[i])
    cosin_sim = 999999
    if norm(np.array(new_obs_foldback_dist)) * norm(np.array(expect_foldback_dist)) !=0:
        cosin_sim = 10 * (1-  np.dot(np.array(expect_foldback_dist), np.array(new_obs_foldback_dist)) / (
                    norm(np.array(new_obs_foldback_dist)) * norm(np.array(expect_foldback_dist))) )
    euclidean = 7 * euclidean_distance(expect_foldback_dist, new_obs_foldback_dist)
    euclidean = euclidean + penalty * counter
    if (segment_number - merged_segments) == 1 :
        return [4,4,4, 4, segment_data, 4, 4 ]
    return [(euclidean + cn_score + float(args.segscore)) / ((segment_number - merged_segments) ** (float(args.norm))),
            cn_score,
            foldback_score, float(args.segscore), segment_data, cosin_sim, euclidean ]
def convert_string(convertor, bfb):
    ans = ''
    for i in bfb:
        ans = ans + convertor[i]
    return ans

with open(args.segments, 'rb') as f:
    segments = pickle.load(f)
bfbfinder = parse_BFBFinderOutput(args.bfbfinder)
expected = parse_segments(args.answer)
convertor = {}
segment_alphabet = list(expected.keys())
for i in range(len(segment_alphabet)):
    convertor[segment_alphabet[i]] = segment_alphabet[len(segment_alphabet) - i - 1]
structure = []
for bfb in bfbfinder:
    if str(args.arm) == 'p':
        bfb = convert_string(convertor, bfb)
    bfb = miroring_bfb(bfb, int(args.m))
    structure.append((bfb, calculate_score(expected, bfb)))
structure = sorted(structure, key=lambda x: x[1][0])
count = 0
foldback = parse_foldbacks(args.foldback)
with open(args.output, 'w') as f:
    while structure and count < 20:
        a = structure.pop(0)
        if count == 0:
            f.write('Structure\tFinalScore\tCNScore\tFNScore\tSegmentScore\n')
        f.write(''.join(i for i in a[0]) + '\t' + str(a[1][0]) + '\t' + str(a[1][1]) + '\t' + str(a[1][2]) + '\t' + str(
            a[1][3]))
        f.write('\n')
        if count == 0:
            with open(args.bfbfinder[:-13] + 'score.csv', 'a') as file:
                file.write('Name\tpenalty' + '\t' + 'normalization' + '\tStructure\tFinalScore\tCNScore\tFNScore\tSNScire\tCosineSim\tEuclidianDistance\n')
                file.write(args.name + '\t'+
                    args.p + '\t' + args.norm + '\t' + ''.join(i for i in a[0]) + '\t' + str(a[1][0]) + '\t' + str(
                        a[1][1]) + '\t' + str(a[1][2]) + '\t' + str(a[1][3]) + '\t' + str(a[1][5]) + '\t' + str(a[1][6]) + '\n')
                file.close()
        count += 1
visualization_cmd = 'python3 BFB_vis.py -sg {sg} -sc {sc} -centro {centro} -rcmap {rcmap} -foldback {foldback} -o {output} '.format(
    sg = args.answer,foldback = args.foldback,sc = args.output, centro = args.centro, rcmap = args.rcmap, output = args.output[:-4])
print(visualization_cmd)
os.system(visualization_cmd)  

