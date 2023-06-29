import matplotlib.pyplot as plt
import argparse
from collections import defaultdict
import pandas as pd
from parsers import SmapEntry, parse_rmcap, parse_smap, parse_xmap, parse_centro, generate_mol_contig, molecule_support_bp, parse_occ, parse_segmentation
from operator import itemgetter
import csv, math
import os, shutil
import string
import pickle
import math
import numpy as np

global xmap
global fig_counter
global contigs
global c_index

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--rcmap", help="rcamp dir", required=True)
parser.add_argument("-c", "--centro", help="hg38 centromeric region", required=True)
parser.add_argument("-n", "--name", help="Name", required=True)
parser.add_argument("-o", "--output", help="Output_dir", required=True)
parser.add_argument("-s", "--smap", help="smap_dir", required=False)
parser.add_argument("-f", "--fsv", help="Fandom_sv dir", required=False)
parser.add_argument("-x", "--xmap", help="xmap dir", required=False)
parser.add_argument("-cmap", "--cmap", help="cmap dir", required=False)
parser.add_argument("-fol", "--folderalignment", help="Folder alignment", required=False)
parser.add_argument("-cov", "--coverage", help="Bionano sample Coverage default is 77", required=False)
parser.add_argument("-bfbfinder", "--bfbfinder", help="Path to BFBFinder Jar file", required=True)
args = parser.parse_args()
bfb_mode = False
if args.folderalignment is not None:
    bfb_mode = True

fig_counter = 1
c_index = 1 
OCC = 77 
amplicon_threshold = 3 #Minimum CN for detecting amplicon
if args.coverage is not None:
    OCC = int(float(args.coverage))


#################################### parsing FaNDOM SV file ##########################################
def parse_fandom_sv(sv_dir):
    bfb_count = {}
    bfb_count = defaultdict(lambda: [], bfb_count)
    with open(sv_dir, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                sv_type = line[6]
                if sv_type == 'duplication_inverted':
                    chrom = 'chr' + line[0]
                    start_pos = float(line[1]) - 4 #what this 4 is?
                    end_pos = float(line[4]) + 4
                    bfb_count[chrom].append([start_pos, end_pos])
    return bfb_count


#################################### generating bfb coverage on reference ##########################################
def generate_bfb(bfb_count, p_cop):
    bfb = {}
    bfb = defaultdict(lambda: {}, bfb)
    for k in p_cop.keys():
        bfb[k] = {i: [] for i in p_cop[k].keys()}
    for k in bfb_count.keys():
        for b in bfb_count[k]:
            mol_list = molecule_support_bp(b.xmap_id1, b.xmap_id2, contigs, xmap)
            for i in bfb[k].keys():
                if b.ref_start <= i <= b.ref_end:
                    bfb[k][i].extend(mol_list)
    for k in bfb.keys():
        remove = []
        for i in bfb[k].keys():
            if len(bfb[k][i]) == 0:
                remove.append(i)
            else:
                bfb[k][i] = len(set(bfb[k][i]))
        for i in remove:
            bfb[k].pop(i)
    return bfb


#################################### for drawing zoom bfb version cluster neighborhood bfbs and return them as a list ##########################################
def cluster_bfb_maps(bfb_maps):
    bfb_maps.sort(key=lambda x: x.ref_start)
    ans = []
    ans_clusters = []
    max_limit = 2800000  # Normal 2000000
    start = bfb_maps[0].ref_start
    for b in bfb_maps:
        if start <= b.ref_end <= start + max_limit:
            ans.append(b)
        else:
            ans.sort(key=lambda x: x.q_id)
            ans_clusters.append(ans)
            ans = [b]
            start = b.ref_start
    ans.sort(key=lambda x: x.q_id)
    ans_clusters.append(ans)
    return ans_clusters

def detect_fold_back_segment(start, direction, segments_cordinate):
    select = ''
    dist = 9999999999
    if direction == 'right':
        for k in segments_cordinate.keys():
            if len(segments_cordinate[k]) > 0:
                if abs(segments_cordinate[k][1] - start) < dist:
                    dist = abs(segments_cordinate[k][1] - start)
                    select = k
    else:
        for k in segments_cordinate.keys():
            if len(segments_cordinate[k]) > 0:
                if abs(segments_cordinate[k][0] - start) < dist:
                    dist = abs(segments_cordinate[k][0] - start)
                    select = k
    return select


#####################################################################################################
def detect_extract_foldback(bfb_cluster, chrom, p_cop, out, segments_cordinate,
                    right_foldback_contig, right_number, left_foldback_contig, left_number, translocation,
                    segmentation_score):
    global fig_counter
    colors = ['red', 'olive', 'teal', 'purple', 'orange', 'dodgerblue', 'sienna', 'coral', 'tomato', 'orchid']
    window_limit = 300000
    min_x = bfb_cluster[0].ref_start - window_limit
    max_x = bfb_cluster[-1].ref_end + window_limit
    cn_dict = {}
    for i in p_cop[chrom].keys():
        if min_x <= int(i) <= max_x:
            cn_dict[i] = p_cop[chrom][i]
    max_cn = max(list(cn_dict.values()))
    min_coverage_for_fold_back = 99999
    if len([i for i in list(left_number.values()) if i > 0]) > 0:
        min_coverage_for_fold_back = min(min([i for i in list(left_number.values()) if i > 0]),
                                         min_coverage_for_fold_back)
    if len([i for i in list(right_number.values()) if i > 0]) > 0:
        min_coverage_for_fold_back = min(min([i for i in list(right_number.values()) if i > 0]),
                                         min_coverage_for_fold_back)
    if min_coverage_for_fold_back == 0:
        min_coverage_for_fold_back = OCC
    with open(args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_foldback_coordinate.txt', 'w') as file:
        file.write('#contig_id\tdirection\tsegment\tchromosome\tstart\tend\n')
        for k in left_foldback_contig.keys():
            if len(left_foldback_contig[k]) > 0:
                p_start, p_end = 9999999999, 0
                for foldback in left_foldback_contig[k]:
                    q_id = foldback.q_id
                    alignment1 = xmap[str(foldback.xmap_id1)]
                    alignment2 = xmap[str(foldback.xmap_id2)]
                    if min(p_end, max(alignment1['RefEndPos'], alignment2['RefEndPos'])) < max(p_start, min(
                            alignment1['RefStartPos'], alignment2['RefStartPos'])):
                        counter = float(segments_cordinate[k][2]) + max(0.5, int(0.05 * max_cn))
                        diff_counter = max(0.5, 0.08 * (max_cn + 3))
                        if min(alignment1['QryStartPos'], alignment1['QryEndPos']) > min(alignment2['QryStartPos'],
                                                                                         alignment2['QryEndPos']):
                            alignment1, alignment2 = alignment2, alignment1
                        if min(alignment1['RefStartPos'], alignment2['RefStartPos']) < min_x:
                            min_x = min(alignment1['RefStartPos'], alignment2['RefStartPos']) - 250000
                        if max(alignment1['RefEndPos'], alignment2['RefEndPos']) > max_x:
                            max_x = max(alignment1['RefEndPos'], alignment2['RefEndPos']) + 250000
                        if alignment1['Orientation'] == '+':
                            seg = detect_fold_back_segment(min(alignment1['RefStartPos'], alignment2['RefStartPos']),
                                                           'left', segments_cordinate)
                            file.write(
                                "{id}\t{dir}\t{segment}\t{chromosome}\t{start}\t{end}\n".format(id=q_id, dir='left',
                                                                                                chromosome=chrom,
                                                                                                segment=seg,
                                                                                                start=alignment1[
                                                                                                    'RefEndPos'],
                                                                                                end=alignment2[
                                                                                                    'RefEndPos']))
                        else:
                            seg = detect_fold_back_segment(min(alignment1['RefStartPos'], alignment2['RefStartPos']),
                                                           'left', segments_cordinate)
                            file.write(
                                "{id}\t{dir}\t{segment}\t{chromosome}\t{start}\t{end}\n".format(id=q_id, dir='left',
                                                                                                chromosome=chrom,
                                                                                                segment=seg,
                                                                                                start=alignment1[
                                                                                                    'RefStartPos'],
                                                                                                end=alignment2[
                                                                                                    'RefStartPos']))
                    p_start = min(p_start, min(alignment1['RefStartPos'], alignment2['RefStartPos']))
                    p_end = max(p_end, max(alignment1['RefEndPos'], alignment2['RefEndPos']))
        
        for k in right_foldback_contig.keys():
            if len(right_foldback_contig[k]) > 0:
                p_start, p_end = 9999999999, 0
                for foldback in right_foldback_contig[k]:
                    q_id = foldback.q_id
                    alignment1 = xmap[str(foldback.xmap_id1)]
                    alignment2 = xmap[str(foldback.xmap_id2)]
                    if min(p_end, max(alignment1['RefEndPos'], alignment2['RefEndPos'])) < max(p_start, min(
                            alignment1['RefStartPos'], alignment2['RefStartPos'])):
                        if len(left_foldback_contig[k]) > 0:
                            counter = float(segments_cordinate[k][2]) + max(0.5,
                                                                            int(0.05 * max_cn)) + 1.5 * max(
                                0.5, 0.08 * (max_cn + 3))
                        else:
                            counter = float(segments_cordinate[k][2]) + max(0.5, int(0.05 * max_cn))
                        diff_counter = 0.08 * (max_cn + 3)
                        if min(alignment1['QryStartPos'], alignment1['QryEndPos']) > min(alignment2['QryStartPos'],
                                                                                         alignment2['QryEndPos']):
                            alignment1, alignment2 = alignment2, alignment1
                        if min(alignment1['RefStartPos'], alignment2['RefStartPos']) < min_x:
                            min_x = min(alignment1['RefStartPos'], alignment2['RefStartPos']) - 250000
                        if max(alignment1['RefEndPos'], alignment2['RefEndPos']) > max_x:
                            max_x = max(alignment1['RefEndPos'], alignment2['RefEndPos']) + 250000
                        if alignment1['Orientation'] == '+':
                            seg = detect_fold_back_segment(max(alignment1['RefEndPos'], alignment2['RefEndPos']),
                                                           'right', segments_cordinate)
                            file.write(
                                "{id}\t{dir}\t{segment}\t{chromosome}\t{start}\t{end}\n".format(id=q_id, dir='right',
                                                                                                segment=seg,
                                                                                                chromosome=chrom,
                                                                                                start=alignment1[
                                                                                                    'RefEndPos'],
                                                                                                end=alignment2[
                                                                                                    'RefEndPos']))
                        else:
                            seg = detect_fold_back_segment(max(alignment1['RefEndPos'], alignment2['RefEndPos']),
                                                           'right', segments_cordinate)
                            file.write(
                                "{id}\t{dir}\t{segment}\t{chromosome}\t{start}\t{end}\n".format(id=q_id, dir='right',
                                                                                                segment=seg,
                                                                                                chromosome=chrom,
                                                                                                start=alignment1[
                                                                                                    'RefStartPos'],
                                                                                                end=alignment2[
                                                                                                    'RefStartPos']))




#################################### generate map between contigs and molecules in each baleb in contigs ##########################################
def check_overlap(a, b, c, d):
    if max(a, c) <= min(b, d):
        return True
    return False


def check_fold_back_dir(xmap1, xmap2):
    read1 = xmap[str(xmap1)]
    read2 = xmap[str(xmap2)]
    if min(float(read1['QryStartPos']), float(read1['QryEndPos'])) > min(float(read2['QryStartPos']),
                                                                         float(read2['QryEndPos'])):
        read2, read1 = read1, read2
    if read1['Orientation'] == '+':
        return 'right'
    else:
        return 'left'


def alignment_CN(xmap1, xmap2):
    read1 = xmap[str(xmap1)]
    contig_id = int(read1['QryContigID'])
    mol_list = molecule_support_bp(int(xmap1), int(xmap2), contigs, xmap)
    return len(mol_list) / OCC


################## return all fold-back in specific region
class FoldBack():
    xmap_id1 = ''
    xmap_id2 = ''
    cn1 = 0
    cn2 = 0
    q_id = ''
    diff = 100000000


def fold_back_overlap(chrom, start, end):
    ans_left = []
    ans_right = []
    for bfb_pair in bfb_count[chrom]:
        if check_overlap(start, end, bfb_pair.ref_start, bfb_pair.ref_end):
            CN = alignment_CN(bfb_pair.xmap_id1, bfb_pair.xmap_id2)
            f = FoldBack()
            f.xmap_id1 = bfb_pair.xmap_id1
            f.cn1 = CN
            f.xmap_id2 = bfb_pair.xmap_id2
            f.cn2 = CN
            f.q_id = bfb_pair.q_id
            if check_fold_back_dir(bfb_pair.xmap_id1, bfb_pair.xmap_id2) == 'left':
                ans_left.append(f)
            else:
                ans_right.append(f)
    return ans_left, ans_right


def detect_closest_segment(segments, point, orientation):
    distance = 1000000000
    ans = '-'
    for k in segments:
        if len(segments[k]) > 0:
            if orientation == 'right':
                if abs(point - segments[k][1]) < distance and abs(point - segments[k][1]) < 0.7 * (
                        segments[k][1] - segments[k][0]):
                    distance = abs(point - segments[k][1])
                    ans = k
            if orientation == 'left':
                if abs(point - segments[k][0]) < distance and abs(point - segments[k][0]) < 0.7 * (
                        segments[k][1] - segments[k][0]):
                    distance = abs(point - segments[k][0])
                    ans = k
    return ans


def calculate_segments_fold_back_support(segments, foldback, orientation):
    radios = 300000
    ans = {}
    for k in segments.keys():
        if len(segments[k]) > 0:
            ans[k] = []
            for f in foldback:
                if orientation == 'right':
                    pos1 = int(float(xmap[str(f.xmap_id1)]['RefEndPos']))
                    pos2 = int(float(xmap[str(f.xmap_id2)]['RefEndPos']))
                    f.diff = abs(pos2 - pos1)
                    if detect_closest_segment(segments, (max(pos2, pos1)), 'right') == k:
                        ans[k].append(f)
                if orientation == 'left':
                    pos1 = int(float(xmap[str(f.xmap_id1)]['RefStartPos']))
                    pos2 = int(float(xmap[str(f.xmap_id2)]['RefStartPos']))
                    f.diff = abs(pos2 - pos1)
                    if detect_closest_segment(segments, (min(pos2, pos1)), 'left') == k:
                        ans[k].append(f)
    segments_to_number_support = {}
    for k in ans.keys():
        mol_list = []
        ans[k] = sorted(ans[k], key=lambda x: x.diff)
        for f in ans[k]:
            mol_list.extend(molecule_support_bp(f.xmap_id1, f.xmap_id2, contigs, xmap))
        segments_to_number_support[k] = len(set(mol_list))
    return ans, segments_to_number_support


def merge_segmentation(segments):
    ans = dict.fromkeys(string.ascii_uppercase, [])
    m_min = 99999999999
    m_max = 0
    for k in segments.keys():
        if len(segments[k]) > 0:
            m_max = max(m_max, max(segments[k][0], segments[k][1]))
            m_min = min(m_min, min(segments[k][0], segments[k][1]))
    region_length = m_max - m_min
    seg_length = abs(segments['A'][0] - segments['A'][1])
    i = 0
    j = 0
    merge_prev = False
    while len(segments[list(segments.keys())[i]]) > 0:
        k = list(segments.keys())[i]
        seg_length = abs(segments[k][0] - segments[k][1])
        merge = False
        if len(segments[list(segments.keys())[i + 1]]) > 0:
            next_seg = segments[list(segments.keys())[i + 1]]
        else:
            next_seg = [0, 1, 999999]
        if merge_prev:
            merge_prev = False
            prev_seg = segments[list(segments.keys())[i - 1]]

            ans[list(segments.keys())[j]] = [min(prev_seg[0], segments[k][0]), max(prev_seg[1], segments[k][1]),
                                             segments[k][2]]
            j += 1
        else:
            if j == 0:
                prev_seg = [9999999999, 1, 99999]
            else:
                prev_seg = ans[list(segments.keys())[j - 1]]
            if seg_length < min(0.05 * region_length, 150000):# or seg_length < 40000:
                merge = True
            elif seg_length < min(0.07 * region_length, 200000):
                if prev_seg[2] > segments[k][2] > next_seg[2] or prev_seg[2] < segments[k][2] < next_seg[2]:
                    merge = True
            if merge:
                if abs(segments[k][2] - prev_seg[2]) < abs(segments[k][2] - next_seg[2]):
                    ans[list(segments.keys())[j - 1]] = [min(prev_seg[0], segments[k][0]),
                                                         max(prev_seg[1], segments[k][1]), prev_seg[2]]
                else:
                    merge_prev = True
            else:
                ans[list(segments.keys())[j]] = segments[list(segments.keys())[i]]
                j += 1
        i += 1
    return ans


def extend_candidate_region(chrom, candidate, left_foldback, right_foldback):
    max_x = max(candidate)
    min_x = min(candidate)
    for l in left_foldback + right_foldback:
        if float(xmap[str(l.xmap_id1)]['RefEndPos']) > max_x:
            max_x = float(xmap[str(l.xmap_id1)]['RefEndPos'])
        if float(xmap[str(l.xmap_id2)]['RefEndPos']) > max_x:
            max_x = float(xmap[str(l.xmap_id2)]['RefEndPos'])
        if float(xmap[str(l.xmap_id1)]['RefStartPos']) < min_x:
            min_x = float(xmap[str(l.xmap_id1)]['RefStartPos'])
        if float(xmap[str(l.xmap_id2)]['RefStartPos']) < min_x:
            min_x = float(xmap[str(l.xmap_id2)]['RefStartPos'])
    candidate = []
    for i in p_cop[chrom].keys():
        if i >= min_x and i <= max_x:
            candidate.append(i)
    return list(set(candidate))


def extract_translocation_orientation(bp):
    alignment1 = xmap[str(bp.xmap_id1)]
    alignment2 = xmap[str(bp.xmap_id2)]
    orient1 = alignment1['Orientation']
    if alignment2['Orientation'] == '+':
        orient2 = '-'
    else:
        orient2 = '+'
    return orient1, orient2


def tune_cn(segments_cordinate_k, chrom):
    score = 1000000
    ans = 0
    for cn in range(math.floor(0.8 * segments_cordinate_k[2]), math.ceil(1.2 * segments_cordinate_k[2]) + 1, math.ceil(
            (math.ceil(1.2 * segments_cordinate_k[2]) + 1 - math.floor(0.8 * segments_cordinate_k[2])) / 10)):
        sum_dif = 0
        for ii in range(len(list(p_cop[chrom].keys())) - 1):
            pos = list(p_cop[chrom].keys())[ii]
            n_pos = list(p_cop[chrom].keys())[ii + 1]
            if segments_cordinate_k[0] <= pos < segments_cordinate_k[1]:
                sum_dif += abs(p_cop[chrom][pos] - cn) * abs(n_pos - pos)
            if pos >= segments_cordinate_k[1]:
                if cn != 0:
                    sum_dif = sum_dif / (
                            math.sqrt(cn) * (
                            segments_cordinate_k[1] - segments_cordinate_k[0]))
                else:
                    sum_dif = 99999
                if sum_dif < score:
                    ans = cn
                    score = sum_dif
                break
    return ans


def reconstruct_bfb(chrom, candidate, left_foldback, right_foldback, deletions, translocation):
    arm = ''
    print('\n# New Candidate')
    candidate = extend_candidate_region(chrom, candidate, left_foldback, right_foldback)
    print(
        'SequenceEdge:\tStartPosition,\tEndPosition,\tPredictedCopyCount,\tAverageCoverage,\tSize,\tNumberReadsMapped')
    segment_list = []
    for l in left_foldback + right_foldback:
        segment_list.append(int(float(xmap[str(l.xmap_id1)]['RefStartPos'])))
        segment_list.append(int(float(xmap[str(l.xmap_id2)]['RefStartPos'])))
        segment_list.append(int(float(xmap[str(l.xmap_id1)]['RefEndPos'])))
        segment_list.append(int(float(xmap[str(l.xmap_id2)]['RefEndPos'])))
    for l in deletions:
        segment_list.append(int(l.ref_start))
        segment_list.append(int(l.ref_end))
    sequence_save = [[]]
    prev_pos = min(candidate)
    with open(args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Chrom', 'Start', 'End', 'CN', 'Length'])
        for i in range(len(list(p_cop[chrom].keys()))):
            pos = list(p_cop[chrom].keys())[i]
            if min(candidate) < float(pos) <= max(candidate):
                if p_cop[chrom][prev_pos] != p_cop[chrom][pos] or int(float(pos)) in segment_list or pos == max(
                        candidate):
                    if prev_pos == min(candidate):
                        prev_pos -= 1
                    weight_affine = 0
                    if (pos - prev_pos - 1) < 10000:
                        weight_affine = float((pos - prev_pos - 1) / 10000)
                    else:
                        weight_affine = 1
                    csvwriter.writerow([chrom, prev_pos + 1, pos, p_cop[chrom][pos], weight_affine])
                    sequence_save.append([prev_pos + 1, pos])
                    prev_pos = pos
    r_cmd = 'Rscript segmentation.R {input} {output}'.format(
        input=args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '.csv',
        output=args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_Rsegmentation.txt')
    os.system(r_cmd)
    segments = parse_segmentation(args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_Rsegmentation.txt')
    segments_cordinate = dict.fromkeys(string.ascii_uppercase, [])
    for i in string.ascii_lowercase:
        segments_cordinate[i] = []
    for i in range(len(segments)):
        s = segments[i]
        print(
            "sequence\t{chr1}:{pos1}-\t{chr1}:{pos2}+\t{cn}\t{cv}\t1\t1".format(chr1=chrom, pos1=sequence_save[s[0]][0],
                                                                                pos2=sequence_save[s[1]][1], cn=s[2],
                                                                                cv=int(s[2]) * OCC))
        segments_cordinate[list(segments_cordinate.keys())[i]] = [sequence_save[s[0]][0], sequence_save[s[1]][1], s[2]]
    while segments_cordinate != merge_segmentation(segments_cordinate):
        segments_cordinate = merge_segmentation(segments_cordinate)
    segmentation_score = 0
    for k in segments_cordinate.keys():
        if len(segments_cordinate[k]) > 0:
            segments_cordinate[k][2] = tune_cn(segments_cordinate[k], chrom)
            sum_dif = 0
            for ii in range(len(list(p_cop[chrom].keys())) - 1):
                pos = list(p_cop[chrom].keys())[ii]
                n_pos = list(p_cop[chrom].keys())[ii + 1]
                if segments_cordinate[k][0] <= pos < segments_cordinate[k][1]:
                    sum_dif += abs(p_cop[chrom][pos] - segments_cordinate[k][2]) * abs(n_pos - pos)
                if pos >= segments_cordinate[k][1]:
                    sum_dif = sum_dif / (
                            math.sqrt(segments_cordinate[k][2]) * (
                            segments_cordinate[k][1] - segments_cordinate[k][0]))
                    segmentation_score += sum_dif
                    break
    print(
        'BreakpointEdge:\tStartPosition->EndPosition,\tPredictedCopyCount,\tNumberOfReadPairs,\tHomologySizeIfAvailable(<0ForInsertions),\tHomology/InsertionSequence')
    for l in left_foldback:
        if xmap[str(l.xmap_id1)]['Orientation'] == '+':
            temp = l.xmap_id1
            l.xmap_id1 = l.xmap_id2
            l.xmap_id2 = temp
            # l = (l[2],l[3],l[0],l[1],l[4])
        print("discordant\t{chr1}:{pos1}-->{chr1}:{pos2}-\t{cn}\t1\tNone\tNone".format(chr1=chrom, pos1=int(
            float(xmap[str(l.xmap_id1)]['RefStartPos'])),
                                                                                       pos2=int(float(
                                                                                           xmap[str(l.xmap_id2)][
                                                                                               'RefStartPos'])),
                                                                                       cn=l.cn1))
    for l in right_foldback:
        if xmap[str(l.xmap_id1)]['Orientation'] == '-':
            # l = (l[2],l[3],l[0],l[1],l[4])
            temp = l.xmap_id1
            l.xmap_id1 = l.xmap_id2
            l.xmap_id2 = temp
        print("discordant\t{chr1}:{pos1}+->{chr1}:{pos2}+\t{cn}\t1\tNone\tNone".format(chr1=chrom, pos1=int(
            float(xmap[str(l.xmap_id1)]['RefEndPos'])),
                                                                                       pos2=int(float(
                                                                                           xmap[str(l.xmap_id2)][
                                                                                               'RefEndPos'])),
                                                                                       cn=l.cn1))
    for l in deletions:
        print(
            "discordant\t{chr1}:{pos1}+->{chr1}:{pos2}-\t{cn}\t1\tNone\tNone".format(chr1=chrom, pos1=int(l.ref_start),
                                                                                     pos2=int(l.ref_end), cn=1))
    for trans in translocation:
        orient1, orient2 = extract_translocation_orientation(trans)
        print("discordant\tchr{chr1}:{pos1}{o1}->chr{chr2}:{pos2}{o2}\t{cn}\t1\tNone\tNone".format(chr1=trans.ref_c_id1,
                                                                                                   chr2=trans.ref_c_id2,
                                                                                                   pos1=int(
                                                                                                       trans.ref_start),
                                                                                                   pos2=int(
                                                                                                       trans.ref_end),
                                                                                                   cn=1, o1=orient1,
                                                                                                   o2=orient2))
    bfb_fold = []
    m_min = 9999999999999
    m_max = 0
    for i in range(1, len(sequence_save)):
        if sequence_save[i][0] < m_min:
            m_min = sequence_save[i][0]
        if sequence_save[i][0] > m_max:
            m_max = sequence_save[i][1]
    for l in left_foldback + right_foldback:
        fold = SmapEntry()
        fold.ref_start, fold.ref_end = sequence_save[1][0], sequence_save[-1][1]
        fold.q_id = l.q_id
        fold.xmap_id1, fold.xmap_id2 = l.xmap_id1, l.xmap_id2
        bfb_fold.append(fold)
    left_foldback_contig, left_number = calculate_segments_fold_back_support(segments_cordinate, left_foldback, 'left')
    right_foldback_contig, right_number = calculate_segments_fold_back_support(segments_cordinate, right_foldback,
                                                                               'right')
    window_limit = 300000
    max_x = bfb_fold[-1].ref_end + window_limit
    if max_x < min(centro[str(chrom)]):
        arm = 'p'
    else:
        arm = 'q'
    cn_list = []
    with open(args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_ans.txt', 'w') as f:
        print('#Segment\tChromosome\tStartPos\tEndPos\tCN\tRightFoldIds\tRightCount\tLeftFoldIds\tLeftCount')
        f.write('#Segment\tChromosome\tStartPos\tEndPos\tCN\tRightFoldIds\tRightCount\tLeftFoldIds\tLeftCount\n')
        for k in segments_cordinate.keys():
            if len(segments_cordinate[k]) > 0:
                print(
                    '{segment}\t{chrom}\t{start}\t{end}\t{CN}\t{right_supp_ids}\t{rightcount}\t{left_supp_ids}\t{leftcount}'.format(
                        segment=k, chrom=chrom, start=segments_cordinate[k][0], end=segments_cordinate[k][1],
                        CN=segments_cordinate[k][2],
                        right_supp_ids=max('-', ','.join(str(i.q_id) for i in right_foldback_contig[k])),
                        rightcount=int(right_number[k]) / OCC,
                        left_supp_ids=max('-', ','.join(str(i.q_id) for i in left_foldback_contig[k])),
                        leftcount=int(left_number[k]) / OCC))
                cn_list.append(round(float(segments_cordinate[k][2])))
                f.write(
                    '{segment}\t{chrom}\t{start}\t{end}\t{CN}\t{right_supp_ids}\t{rightcount}\t{left_supp_ids}\t{leftcount}\n'.format(
                        segment=k, chrom=chrom, start=segments_cordinate[k][0], end=segments_cordinate[k][1],
                        CN=segments_cordinate[k][2],
                        right_supp_ids=max('-', ','.join(str(i.q_id) for i in right_foldback_contig[k])),
                        rightcount=int(right_number[k]) / OCC,
                        left_supp_ids=max('-', ','.join(str(i.q_id) for i in left_foldback_contig[k])),
                        leftcount=int(left_number[k]) / OCC))
    f.close()
    times = 0
    while sum(cn_list) > 50:
        cn_list = [i / 2 for i in cn_list]
        times += 1
    cn_list = [round(i) for i in cn_list]
    if arm == 'p':
        cn_list = cn_list[::-1]
    with open(args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_segments_dict.pkl', 'wb') as seg_dict_file:
        pickle.dump(segments_cordinate, seg_dict_file)
    bfbfinder_cmd = 'java -jar {bfb_dir} -a -s -e=PoissonErrorModel -w=0.63 [{cn}] > {outputdir}'.format(
        bfb_dir=args.bfbfinder,
        cn=','.join(str(i) for i in cn_list),
        outputdir=args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_ans_BFBFinder.txt')
    print(bfbfinder_cmd)
    os.system(bfbfinder_cmd)
    detect_extract_foldback(bfb_fold, chrom, p_cop, args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_ans.png',
                    segments_cordinate, right_foldback_contig, right_number, left_foldback_contig,
                    left_number, translocation, segmentation_score)
    for pen in [i * 0.5 for i in range(4, 5)]:
        for norm in [i * 0.2 for i in range(5, 6)]:
            # os.mkdir(
            #     args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_' + str(float(pen)) + '_' + str(float(norm))[
            #                                                                                         :3])
            os.mkdir(args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom))
            result_analysis_cmd = 'python3 analyzing_BFBFinder.py -foldback {fold} -bfbf {bfbfinder} -a {expected} -o {out} -arm {arm} -seg {seg} -m {mul} -p {pen} -norm {norm} -segscore {segscore} -name {name} -centro {centro} -rcmap {rcmap}'.format(
                bfbfinder=args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_ans_BFBFinder.txt',
                expected=args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_ans.txt',
                out=args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_Final_ans.txt', arm=arm,
                seg=args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_segments_dict.pkl', mul=times, pen=pen,
                norm=norm, segscore=segmentation_score,
                fold=args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_foldback_coordinate.txt',
                name=args.name + '_' + 'amplicon' + str(c_index) + str(chrom), centro=args.centro, rcmap=args.rcmap)
            print(result_analysis_cmd)
            os.system(result_analysis_cmd)
    remove_cm = 'rm '+args.output + '/' + 'amplicon' + str(c_index)+'_' + str(chrom) + '_segments_dict.pkl'
    os.system(remove_cm)


def detect_deletions(chrom, start, end):
    ans = []
    for bp in smap:
        if bp.sv_type.startswith('dele'):
            if chrom == 'chr' + str(bp.ref_c_id1):
                if start <= bp.ref_start <= end and start <= bp.ref_end <= end:
                    ans.append(bp)

    return ans


def merge_translocation(translocation_list):
    ans = []
    for i in range(len(translocation_list)):
        find = False
        t1 = translocation_list[i]
        for j in range(i + 1, len(translocation_list)):
            t2 = translocation_list[j]
            if (t1.ref_c_id1 == t2.ref_c_id1 and t1.ref_c_id2 == t2.ref_c_id2 and abs(
                    t1.ref_start - t2.ref_start) < 30000 and abs(t1.ref_end - t2.ref_end) < 30000):
                find = True
                break
        if not find:
            ans.append(t1)
    return ans


def detect_translocation(chrom, start, end):
    ans = []
    for bp in smap:
        if bp.sv_type.startswith('trans') and bp.confidence >= 0.05:
            if (chrom == 'chr' + str(bp.ref_c_id1) and start <= bp.ref_start <= end) or (
                    chrom == 'chr' + str(bp.ref_c_id2) and start <= bp.ref_end <= end):
                ans.append(bp)
    ans = merge_translocation(ans)
    return ans


def find_border_fold_back(foldbakcs):
    min_m = 99999999999999
    max_m = 0
    for f in foldbakcs:
        x1 = xmap[str(f.xmap_id1)]
        x2 = xmap[str(f.xmap_id2)]
        max_m = max(max_m, max(x1['RefStartPos'], x1['RefEndPos'], x2['RefStartPos'], x2['RefEndPos']))
        min_m = min(min_m, min(x1['RefStartPos'], x1['RefEndPos'], x2['RefStartPos'], x2['RefEndPos']))
    return min_m, max_m


def calculate_mean_cn(p_cop, chrom, start, end):
    sum = 0
    count = 0
    for i in p_cop[chrom].keys():
        if start <= i <= end:
            sum += p_cop[chrom][i]
            count += 1
    return sum / count


################## detect bfb candidate
def detect_bfb_candidate(bfb, p_cop):
    global c_index
    max_bfb_length = 1500000  # Important to change
    detected_region = {}
    for chrom in bfb.keys():
        detected_region[chrom] = []
        for region in bfb[chrom].keys():
            if bfb[chrom][region] >= normal_coverage and p_cop[chrom][region] >= amplicon_threshold:  # Important to change
                if len(detected_region[chrom]) == 0:
                    detected_region[chrom].append([region])
                else:
                    if abs(detected_region[chrom][-1][-1] - region) < max_bfb_length or calculate_mean_cn(p_cop, chrom,
                                                                                                          detected_region[
                                                                                                              chrom][
                                                                                                              -1][-1],
                                                                                                          region) > 7:
                        detected_region[chrom][-1].append(region)
                    else:
                        detected_region[chrom].append([region])
    for chrom in detected_region.keys():
#         if chrom == 'chr7':
        for candidate in detected_region[chrom]:
            # if max(candidate)> 90000000:
            # for k in p_cop[chrom].keys():
            #     if 39398385<=int(k)<=40268743 and int(k) not in candidate:
            #         candidate.append(int(k))
            # candidate = sorted(candidate)
            left_foldback, right_foldback = fold_back_overlap(chrom, min(candidate), max(candidate))
            # for i in left_foldback + right_foldback:
            #     print(i.xmap_id1, i.xmap_id2, i.q_id)
            min_m, max_m = find_border_fold_back(left_foldback + right_foldback)
            deletions = detect_deletions(chrom, min_m, max_m)
            translocation = detect_translocation(chrom, min_m - 10000, max_m + 10000)
            if len(candidate) > 1:
                reconstruct_bfb(chrom, candidate, left_foldback, right_foldback, deletions, translocation)
                c_index += 1


def extract_contig_id_to_xmap_entry(xmap):
    ans = {}
    ans = defaultdict(lambda: [], ans)
    for x in xmap.keys():
        ans[xmap[x]['QryContigID']].append(x)
    return ans


def find_xmap_entry(xmap_id):
    return [float(xmap[xmap_id]['QryStartPos']), float(xmap[xmap_id]['QryEndPos']), float(xmap[xmap_id]['RefStartPos']),
            float(xmap[xmap_id]['RefEndPos']), xmap[xmap_id]['Orientation']]


def check_bfb_quality(bfb):
    q_id = str(bfb.q_id)
    alignment1 = xmap[str(bfb.xmap_id1)]
    alignment2 = xmap[str(bfb.xmap_id2)]
    if min(alignment1['QryStartPos'], alignment1['QryEndPos']) > min(alignment2['QryStartPos'],
                                                                     alignment2['QryEndPos']):
        alignment1, alignment2 = alignment2, alignment1
    diff_q_len = abs(int((min(alignment2['QryStartPos'], alignment2['QryEndPos']) - max(alignment1['QryEndPos'],
                                                                                        alignment1[
                                                                                            'QryStartPos'])) / 1000))
    if alignment1['Orientation'] == '+':  # its right foldback
        d_prime = int(abs(alignment1['RefEndPos'] - alignment2['RefEndPos']) / 1000)
    else:  # it is left fold back
        d_prime = int(abs(alignment1['RefStartPos'] - alignment2['RefStartPos']) / 1000)

    if diff_q_len > 120 and abs(d_prime - diff_q_len) > 60:
        return False
    if d_prime > 120 and abs(d_prime - diff_q_len) > 60:
        return False
    if abs(d_prime - diff_q_len) > 100:
        return False
    for xmap_id in contig_to_xmap_entry[q_id]:
        if xmap_id != str(bfb.xmap_id1) and xmap_id != str(bfb.xmap_id2):
            alignment = xmap[str(xmap_id)]
            s = min(alignment1['QryStartPos'], alignment1['QryEndPos'], alignment2['QryStartPos'],
                    alignment2['QryEndPos'])
            e = max(alignment1['QryStartPos'], alignment1['QryEndPos'], alignment2['QryStartPos'],
                    alignment2['QryEndPos'])
            if (s <= alignment['QryStartPos'] <= e or s <= alignment['QryEndPos'] <= e) and alignment[
                'Confidence'] > min(alignment1['Confidence'], alignment2['Confidence']):
                return False
    return True


def remove_unconfident_bfb(bfb_count):
    ans = {}
    for chrom in bfb_count.keys():
        ans[chrom] = []
        for bfb in bfb_count[chrom]:
            if check_bfb_quality(bfb):
                ans[chrom].append(bfb)
    return ans


for f in os.listdir(args.output):
    if os.path.isfile(os.path.join(args.output, f)):
        os.remove(os.path.join(args.output, f))
    elif os.path.isdir(os.path.join(args.output, f)):
        shutil.rmtree(os.path.join(args.output, f))

# ============================================================================================
normal_coverage = 10  # If raw molecule greater than this TH will candidate a foldback
if args.coverage is not None:
    normal_coverage = int (int(float(args.coverage))/10)
#### extracting centromere regions
centro = parse_centro(args.centro)
##################
xmap = parse_xmap(
    args.xmap)  # Output format: return dict of dict with this info : ["XmapEntryID", "QryContigID", "RefContigID", "Orientation", "Confidence", "QryLen", "RefLen", "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", "Alignment"]
##################
contig_to_xmap_entry = extract_contig_id_to_xmap_entry(
    xmap)  # Output format return dict of contig_id to list of xmap id of that contig
##################
bionano_cn_dir = args.rcmap
p_cov, p_cop = parse_rmcap(bionano_cn_dir)
if args.smap is not None:
    bfb_count, smap, _ = parse_smap(
        args.smap)  # return dict to chrom and in each dict list of duplication inverted bfb_count[chrom].append([start_pos, end_pos, q_id, xmap_id1, xmap_id2])
    # bfb_count = remove_unconfident_bfb(bfb_count)
contig_id_list = []
for chrom in bfb_count.keys():
    for i in bfb_count[chrom]:
        contig_id_list.append(int(i.q_id))
contigs = generate_mol_contig(contig_id_list, args.folderalignment)
if args.fsv is not None:
    bfb_count = parse_fandom_sv(args.fsv)
if bfb_mode:
    bfb = generate_bfb(bfb_count, p_cop)
    detect_bfb_candidate(bfb, p_cop)
# generate_bfb_plot(bfb_count, p_cop)
################## ploting main plot###################
fig = plt.figure(figsize=(20, 25))
fig.suptitle(args.name, fontsize=48)
rows = 6
column = 4
grid = plt.GridSpec(rows, column, wspace=.25, hspace=.25)
i = 0
for k in p_cop.keys():
    a = pd.DataFrame(p_cop[k].items(), columns=['Pos', 'Cov'])
    exec(f"plt.subplot(grid{[i]})")
    plt.scatter(list(a['Pos']), list(a['Cov']), s=0.01, label='Copy Number')
    plt.axvline(min(centro[k]), color='green')
    plt.axvline(max(centro[k]), color='green')
    if max(list(a['Cov'])) > 15:
        ax = plt.gca()
        ax.set_yscale('symlog', base=10)
        plt.ylim([-0.5, max(max(list(a['Cov'])), 0) + 50])
        plt.yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150])
    else:
        plt.ylim([-0.5, min(15, max(list(a['Cov']))) + 1])
    ax = plt.gca()
    if bfb_mode:
        ax2 = ax.twinx()
        ax2.scatter(list(bfb[k].keys()), list(bfb[k].values()), color='red', s=5, label='BFB coverage')
        ax2.set_ylim([0, 1400])
    if (i % 4 == 3) and bfb_mode:
        ax2.set_ylabel("bfb count", color="black", labelpad=1, fontsize=15)
    if (i % 4 == 0):
        ax.set_ylabel("CopyNumber", color="black", labelpad=1, fontsize=15)
    plt.title(str(k))
    i += 1
lines, labels = fig.axes[-1].get_legend_handles_labels()
fig.legend(lines, labels, loc='lower center', ncol=2)
plt.savefig(args.output + '/' + args.name + '.png', dpi=300)
plt.close()
###########################################
font = {'family' : 'Arial',
        'size'   : 10}
plt.rc('font', **font)
for k in p_cop.keys():
    plt.clf()
    fig, ax = plt.subplots()
    fig.set_size_inches(3.5, 3)
    a = pd.DataFrame(p_cop[k].items(), columns=['Pos', 'Cov'])
    plt.scatter(list(a['Pos']), list(a['Cov']), s=0.01, label='Copy Number')
    plt.axvline(min(centro[k]), color='green')
    plt.axvline(max(centro[k]), color='green')
    if max(list(a['Cov'])) > 15:
        ax = plt.gca()
        ax.set_yscale('symlog', base=10)
        plt.ylim([-0.5, max(max(list(a['Cov'])), 0) + 50])
        plt.yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150])
    else:
        plt.ylim([-0.5, min(15, max(list(a['Cov']))) + 1])
    ax = plt.gca()
    if bfb_mode:
        ax2 = ax.twinx()
        ax2.scatter(list(bfb[k].keys()), list(bfb[k].values()), color='red', s=5, label='BFB coverage')
        ax2.set_ylim(bottom=0)
    ax2.set_ylabel("#foldback molecules", color="black", labelpad=1, fontsize=11)
    ax.set_ylabel("CopyNumber", color="black", labelpad=1,fontsize = 11)
    ax.set_xlabel('Genomic Coordinate')
    plt.title(str(k))
    plt.savefig(args.output + '/' + args.name  + str(k) + '.png', dpi=300,bbox_inches='tight')
    plt.close()
