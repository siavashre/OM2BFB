import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
import argparse
font = {'family' : 'Arial',
        'size'   : 10}
plt.rc('font', **font)
parser = argparse.ArgumentParser()
parser.add_argument("-sg", "--sg", help="segments coordinates dir", required=True)
parser.add_argument("-sc", "--sc", help="Scores dir", required=True)
parser.add_argument("-centro", "--centro", help="centromere region", required=True)
parser.add_argument("-o", "--output", help="Output_dir", required=True)
parser.add_argument("-rcmap", "--rcmap", help="rcamp dir", required=True)
parser.add_argument("-foldback", "--foldback", help="Foldback Directory", required=True)
args = parser.parse_args()

def parse_segment_coordinates(file_dir):
    segments_coordinates = {}
    with open(file_dir, 'r') as f:        
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                name = line[0]
                chrom =  line[1]
                start = int(line[2])
                end = int(line[3])
                cn = int(line[4])
                segments_coordinates[name] = {'chrom':chrom, 'start':start, 'end':end,'cn':cn}
    return segments_coordinates
def parse_scores2(file_dir):
    with open(file_dir, 'r') as f: 
        for line in f:
            if not line.startswith('Name'):
                line = line.strip().split('\t')
                final_score = float(line[4])
                CN_score = float(line[5])
                FN_score = float(line[6])
                SN_score = float(line[7])
                Cosin_sime = float(line[8])
                Euclidian_sime = float(line[9])
                structure = line[3].strip()
                scores = {'Structure':structure,'Final_score':final_score,'CN_score':CN_score,'FN_score':FN_score,'SN_score':SN_score , 'CosineSim':Cosin_sime,'Euclidian':Euclidian_sime}
    return scores

def parse_scores(file_dir):
    ans = []
    with open(file_dir, 'r') as f: 
        for line in f:
            if not line.startswith('Structure'):
                line = line.strip().split('\t')
                final_score = float(line[1])
                CN_score = float(line[2])
                FN_score = float(line[3])
                SN_score = float(line[4])
                structure = line[0].strip()
                scores = {'Structure':structure,'Final_score':final_score,'CN_score':CN_score,'FN_score':FN_score,'SN_score':SN_score }
                ans.append(scores)
    return ans

def parse_centro(centro_dir):
    r = {}
    r = defaultdict(lambda: [], r)
    with open(centro_dir, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            key = line[0]
            pos1 = int(line[1])
            pos2 = int(line[2])
            r[key].append(pos1)
            r[key].append(pos2)
    return r


def parse_rmcap(cmap_dir):
    cov = {}
    cov = defaultdict(lambda: {}, cov)
    cop = {}
    cop = defaultdict(lambda: {}, cop)
    fcop = {}
    fcop = defaultdict(lambda: {}, fcop)
    with open(cmap_dir, 'r') as f:
        for line in f:
            if line.startswith("#"):
                head = line[1:].rstrip().rsplit()
            if not line.startswith('#'):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                chrom = fD['CMapId']
                pos = int(float(fD['Position']))
                cover = int(fD['Coverage'])
                copynumber = int(fD['CopyNumber'])
                fcn = float(fD['fractionalCopyNumber'])
                cov['chr' + chrom][pos] = cover
                cop['chr' + chrom][pos] = copynumber
                fcop['chr' + chrom][pos] = fcn
    return cov, cop,fcop

def parse_foldback_coordinate(file_dir):
    foldbacks = []
    with open(file_dir, 'r') as f:        
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                fb_type = line[1]
                segment_name = line[2]
                start = float(line[4])
                end = float(line[5])
                chrom = line[3]
                foldbacks.append({'fb_type':fb_type,'segment_name':segment_name,'start':start, 'end':end,'chrom':chrom})
    return foldbacks        

def filter_foldbacks(foldbacks,segments_coordinates):
    filtered_fb = []
    for segment in segments_coordinates.keys():
        for fb_type in ['left','right']:
            dist = 9999999999
            ans = ''
            for fb in foldbacks:
                if fb['fb_type'] == fb_type and fb['segment_name'] == segment:
                    avg_point = (fb['start'] + fb['end'])/2
                    if fb_type == 'left':
                        if abs(avg_point - segments_coordinates[segment]['start']) < dist:
                            dist = abs(avg_point - segments_coordinates[segment]['start'])
                            ans = fb
                    if fb_type == 'right':
                        if abs(avg_point - segments_coordinates[segment]['end']) < dist:
                            dist = abs(avg_point - segments_coordinates[segment]['end'])
                            ans = fb
            if ans !='':
                filtered_fb.append(ans)
    return filtered_fb

def detect_start_end(segments):
    chrom = ''
    start = 9999999999
    end = 0
    for k in segments.keys():
        chrom = segments[k]['chrom']
        if segments[k]['start'] < start:
            start = segments[k]['start']
        if segments[k]['end'] > end:
            end = segments[k]['end']
    return chrom, start-300000, end + 300000

def extract_fcna(fcop,chrom,start,end):
    x = []
    y = []
    for k in fcop[chrom].keys():
        if start<= k <= end:
            x.append(k)
            y.append(fcop[chrom][k])
    return x, y

def plot_segments(segments_coordinates):
    for s in segments_coordinates.values():
        ax.hlines(y = s['cn'], xmin = s['start'], xmax = s['end'],alpha=0.7,color='black', linewidth=1)
    for i, s in enumerate(list(segments_coordinates.keys())):
        prev_s = list(segments_coordinates.keys())[i-1]
        if i >0:
            x1 = segments_coordinates[prev_s]['end']
            y2 = segments_coordinates[prev_s]['cn']
            x2 = segments_coordinates[s]['start']
            y1 = segments_coordinates[s]['cn']
            ax.vlines(ymin =min(y1,y2),ymax = max(y1,y2), x = x1,alpha=0.7,color='black', linewidth=1)

def plot_foldbacks(foldbacks_coordinate,max_y,max_x,start_x):
    prop = dict(arrowstyle="-|>,head_width=0.1,head_length=0.17",
            shrinkA=0,shrinkB=0,color = 'darkblue',alpha = 0.7, linewidth = 0.5)
    arrow_length = 0.03 * max_x
    for foldback in foldbacks_coordinate:
        ax.plot([foldback['start'], foldback['end']], [max_y*1.08,max_y*1.13],color='darkblue',linestyle='--', alpha = 0.6, linewidth=0.5)
        if foldback['fb_type'] == 'left':
            ax.annotate("", xy=(foldback['start'],max_y*1.08), xytext=(foldback['start']+arrow_length,max_y*1.08), arrowprops=prop)
            ax.annotate("", xy=(foldback['end']+arrow_length,max_y*1.13), xytext=(foldback['end'],max_y*1.13), arrowprops=prop)
        if foldback['fb_type'] == 'right':
            ax.annotate("", xy=(foldback['start'],max_y*1.08), xytext=(foldback['start']-arrow_length,max_y*1.08), arrowprops=prop)
            ax.annotate("", xy=(foldback['end']-arrow_length,max_y*1.13), xytext=(foldback['end'],max_y*1.13), arrowprops=prop)

def plot_rectangle_plot(segments_coordinates,max_cn):
    for s in segments_coordinates:
        ax.add_patch(plt.Rectangle((segments_coordinates[s]['start'], 0), segments_coordinates[s]['end'] - segments_coordinates[s]['start'],
                                      - 0.07 * max(1.2 * max_cn, max_cn + 3), edgecolor='r', facecolor='none',
                                      clip_on=False, alpha=1))
        x = 1.1 * ((segments_coordinates[s]['start'] + segments_coordinates[s]['end']) / 2)
        ax.annotate(s, xy=(x - 0.15 * (segments_coordinates[s]['end'] - segments_coordinates[s]['start']),
                                    - 0.07 * max(1.2 * max_cn, max_cn + 3)), weight="bold")
def plot_segments_border(segments_coordinates,max_cn,max_y):
    ax.vlines(ymin = 0 , ymax = max_cn, x = list(segments_coordinates.values())[0]['start'], linewidth = 1, alpha = 0.8 ,linestyle='--', color = 'gray')
    for s in segments_coordinates:
        ax.vlines(ymin = 0 , ymax = max_cn, x = segments_coordinates[s]['end'], linewidth = 1, alpha = 0.8 ,linestyle='--', color = 'gray')
        x = (segments_coordinates[s]['start'] + segments_coordinates[s]['end']) / 2
        ax.annotate(s, xy=(x ,0.02 * max(1.2 * max_y, max_y + 3)), weight="bold", ha = 'center')

def find_in_foldback(segment, direction, foldbacks):
    for f in foldbacks:
        if f['segment_name'] == segment and f['fb_type'] == direction:
            return True
    return False

def plot_structure(score, segments_coordinates, arm,max_y,max_x, foldbacks):
    structure = score['Structure']
    times = 0
    arrow_length = 0.03 * max_x
    prop = dict(arrowstyle="-|>,head_width=0.05,head_length=0.15",
            shrinkA=0,shrinkB=0,color = 'darkblue',alpha = 0.6, linewidth = 0.5)
    rectangle_width = 0.008
    y_increase = 0.028
    color = '#0072b2'
    edgecolor = 'black'
    while(len(structure)>40):
        times += 1 
        structure = structure[0:int(len(structure)/2)]
    direction = 'right'
    prev = ''
    y = max_y* 1.2
    if arm == 'p':
        direction = 'left'
    for i, s in enumerate(structure):
        if s != prev:
            ax.add_patch(
                plt.Rectangle((segments_coordinates[s]['start'], y), segments_coordinates[s]['end'] - segments_coordinates[s]['start'],
                              rectangle_width* max_y, edgecolor=edgecolor, facecolor=color,
                              clip_on=False, alpha=0.6, linewidth=0.3))
            prev = s
        else:
            y = y + y_increase * max_y
            ax.add_patch(
                plt.Rectangle((segments_coordinates[s]['start'], y), segments_coordinates[s]['end'] - segments_coordinates[s]['start'],
                              rectangle_width* max_y, edgecolor=edgecolor, facecolor=color,
                              clip_on=False, alpha=0.6, linewidth=0.3))   
            if direction == 'right':
                linestyle = '--'
                color2 = 'red'
                prop = dict(arrowstyle="-|>,head_width=0.05,head_length=0.15",
            shrinkA=0,shrinkB=0,color = 'red',alpha = 0.6, linewidth = 0.5)
                if find_in_foldback(s,direction,foldbacks):
                    linestyle = '-'
                    color2 = 'darkblue'
                    prop = dict(arrowstyle="-|>,head_width=0.05,head_length=0.15",
            shrinkA=0,shrinkB=0,color = 'darkblue',alpha = 0.6, linewidth = 0.5)
                direction = 'left'
                x_point = [segments_coordinates[s]['end'],segments_coordinates[s]['end'] + arrow_length,segments_coordinates[s]['end'] + arrow_length]
                y_point = [y - y_increase * max_y + (rectangle_width/2) * max_y, y - y_increase * max_y + (rectangle_width/2) * max_y, y + (rectangle_width/2) * max_y]
                ax.plot(x_point, y_point, color=color2, alpha=0.6, linewidth=0.5, linestyle = linestyle)
                ax.annotate("", xy=(segments_coordinates[s]['end'],y + (rectangle_width/2) * max_y), xytext=(segments_coordinates[s]['end'] + arrow_length,y + (rectangle_width/2) * max_y), arrowprops=prop)
            else:
                linestyle = '--'
                color2 = 'red'
                prop = dict(arrowstyle="-|>,head_width=0.05,head_length=0.15",
            shrinkA=0,shrinkB=0,color = 'red',alpha = 0.6, linewidth = 0.5)
                if find_in_foldback(s,direction,foldbacks):
                    linestyle = '-'
                    color2 = 'darkblue'
                    prop = dict(arrowstyle="-|>,head_width=0.05,head_length=0.15",
            shrinkA=0,shrinkB=0,color = 'darkblue',alpha = 0.6, linewidth = 0.5)
                direction = 'right'
                x_point = [segments_coordinates[s]['start'],segments_coordinates[s]['start'] - arrow_length,segments_coordinates[s]['start'] - arrow_length]
                y_point = [y - y_increase * max_y + (rectangle_width/2) * max_y, y - y_increase * max_y + (rectangle_width/2) * max_y, y + (rectangle_width/2) * max_y]
                ax.plot(x_point, y_point, color=color2, alpha=0.6, linewidth=0.5, linestyle = linestyle)
                ax.annotate("", xy=(segments_coordinates[s]['start'],y + (rectangle_width/2) * max_y), xytext=(segments_coordinates[s]['start'] - arrow_length,y + (rectangle_width/2) * max_y), arrowprops=prop)
    ax.annotate('x'+str(2**times), xy = (ax.get_xlim()[1]-0.08*max_x, 0.90*ax.get_ylim()[1]),weight = 'bold',ha = 'center')
segments_coordinates = parse_segment_coordinates(args.sg)
reconstructed_structure = ''
all_scores = parse_scores(args.sc)
centro = parse_centro(args.centro)
arm = ''
p_cov, p_cop,p_fcn = parse_rmcap(args.rcmap)
foldbacks_coordinate = parse_foldback_coordinate(args.foldback)
chrom , start , end = detect_start_end(segments_coordinates)
x,y = extract_fcna(p_fcn,chrom,start,end)
if max(x) < min(centro[chrom]):
    arm = 'p'
else:
    arm = 'q'
for index_1, scores in enumerate(all_scores):
    plt.clf()
    fig, ax = plt.subplots()
    fig.set_size_inches(3.5, 3)
    plt.scatter(x, y, c ="#0072b2", s= 0.1,alpha = 0.5)
    plt.stackplot(x, y, color='#d55e00', alpha=0.3)
    plot_segments(segments_coordinates)
    plot_foldbacks(foldbacks_coordinate,max(y),max(x)-min(x), min(x))
    plot_structure(scores, segments_coordinates,arm,max(y),max(x)-min(x),foldbacks_coordinate)
    plt.xlabel('Position')
    plt.ylabel('CopyNumber')
    plot_segments_border(segments_coordinates,ax.get_ylim()[1], max(y))
    ytcik = ax.get_yticks()
    new_ytick = []
    for i in ytcik:
        if i < max(y)*1.1:
            new_ytick.append(i)
    plt.yticks(new_ytick)
    ax.add_patch(plt.Rectangle((ax.get_xlim()[0], max(y)*1.04), ax.get_xlim()[1]-ax.get_xlim()[0],
                                      0.13*max(y), edgecolor='none', facecolor='y',
                                       alpha=0.1))
    plt.subplots_adjust(bottom=0.15, left = 0.175)
    plt.legend(loc ="upper left",title="Score = {score}".format(score = str(scores['Final_score'])[:4]),prop={'size': 4},title_fontsize=6)
    # plt.annotate("Score = {score}".format(score = str(scores['Final_score'])[:4]),xy = (0.1,0.1),xycoords='axes fraction',fontsize=6)
    plt.title(str(chrom)+arm)
    plt.savefig(args.output+'_'+str(index_1+1)+'.png', dpi = 300)
    plt.close()
