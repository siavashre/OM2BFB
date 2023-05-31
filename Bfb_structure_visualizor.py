import matplotlib.pyplot as plt
from collections import defaultdict

def reverse(a):
    if a == 'right':
        return 'left'
    return 'right'

def calculate_cn(bfb, arm):
    d = defaultdict(int)
    for c in bfb:
        d[c] += 1
    current = 'right'
    if arm == 'p':
        current = 'left'
    d_left_right = defaultdict(lambda: defaultdict(int))
    for i in range(len(bfb) - 1):
        if bfb[i] == bfb[i + 1]:
            d_left_right[bfb[i]][current] += 1
            current = reverse(current)
    return d, d_left_right


def bfb_visualizor2(structure, output, segments_cordinate, arm, score, exp_foldback, foldback, fold_back_real_coordinate):
    ax = plt.gca()
    x = 0
    y = 1
    prev = ''
    direction = 'right'
    if arm == 'p':
        direction = 'left'
    min_x = 999999999999999
    max_x = 0
    cn, foldbacknumber = calculate_cn(structure, arm)
    for i, s in enumerate(structure):
        if min_x > min(segments_cordinate[s][0], segments_cordinate[s][1]):
            min_x = min(segments_cordinate[s][0], segments_cordinate[s][1])
        if max_x < max(segments_cordinate[s][0], segments_cordinate[s][1]):
            max_x = max(segments_cordinate[s][0], segments_cordinate[s][1])
    for i, s in enumerate(structure):
        if s != prev:
            x = (segments_cordinate[s][0] + segments_cordinate[s][1]) / 2
            ax.add_patch(
                plt.Rectangle((segments_cordinate[s][0], y), segments_cordinate[s][1] - segments_cordinate[s][0],
                              0.5, edgecolor='r', facecolor='#edc580',
                              clip_on=False, alpha=0.6, linewidth=0.3))
            if y == 1:
                plt.annotate(s, xy=(x - 0.07 * (segments_cordinate[s][1] - segments_cordinate[s][0]),
                                    y + 0.15))
            if direction == 'right' and not fold_back_real_coordinate:
                plt.arrow(x=x + 0.15 * (segments_cordinate[s][1] - segments_cordinate[s][0]), y=y + 0.25,
                          dx=0.3 * (segments_cordinate[s][1] - segments_cordinate[s][0]), dy=0, color='blue', alpha=0.6,
                          length_includes_head=True,
                          head_length=0.03 * (segments_cordinate[s][1] - segments_cordinate[s][0]),
                          head_width=0.2, linewidth=0.5)

            elif not fold_back_real_coordinate:
                plt.arrow(x=x - 0.15 * (segments_cordinate[s][1] - segments_cordinate[s][0]), y=y + 0.25,
                          dx=-0.3 * (segments_cordinate[s][1] - segments_cordinate[s][0]), dy=0, color='blue',
                          alpha=0.6,
                          length_includes_head=True,
                          head_length=0.03 * (segments_cordinate[s][1] - segments_cordinate[s][0]),
                          head_width=0.2, linewidth=0.5)
            prev = s
        else:
            prev = s
            y = y + 1
            x = (segments_cordinate[s][0] + segments_cordinate[s][1]) / 2
            ax.add_patch(
                plt.Rectangle((segments_cordinate[s][0], y), segments_cordinate[s][1] - segments_cordinate[s][0],
                              0.5, edgecolor='r', facecolor='#edc580',
                              clip_on=False, alpha=0.6, linewidth=0.3))
            if direction == 'right':
                if not fold_back_real_coordinate:
                    plt.arrow(x=x - 0.15 * (segments_cordinate[s][1] - segments_cordinate[s][0]), y=y + 0.25,
                          dx=-0.3 * (segments_cordinate[s][1] - segments_cordinate[s][0]), dy=0, color='blue',
                          alpha=0.6,
                          length_includes_head=True,
                          head_length=0.03 * (segments_cordinate[s][1] - segments_cordinate[s][0]),
                          head_width=0.2, linewidth=0.5)
                color = 'blue'
                linestyle = '-'
                if exp_foldback[s]['right'] == 0:
                    linestyle = '--'
                x_point = [segments_cordinate[s][1],
                           segments_cordinate[s][1] + max(0.1 * (segments_cordinate[s][1] - segments_cordinate[s][0]),0.02 * (max_x - min_x)),
                           segments_cordinate[s][1] + max(0.1 * (segments_cordinate[s][1] - segments_cordinate[s][0]),0.02 * (max_x - min_x))]
                arow_x_point = segments_cordinate[s][1] + max(0.1 * (segments_cordinate[s][1] - segments_cordinate[s][0]),0.02 * (max_x - min_x))
                if s in foldback.keys() and fold_back_real_coordinate:
                    if len(foldback[s]['right']) != 0:
                        x_point = [foldback[s]['right'][0][1] - max(0.1 * (segments_cordinate[s][1] - segments_cordinate[s][0]),0.02 * (max_x - min_x)),
                                   foldback[s]['right'][0][1] ,
                                   foldback[s]['right'][0][2]]
                        arow_x_point = foldback[s]['right'][0][2]
                y_point = [y - 1 + 0.25, y - 1 + 0.25, y + 0.25]
                plt.plot(x_point, y_point, color=color, alpha=0.6, linewidth=1, linestyle = linestyle)
                plt.arrow(x=arow_x_point,
                          y=y + 0.25,
                          dx=-1 * max(0.1 * (segments_cordinate[s][1] - segments_cordinate[s][0]),0.02 * (max_x - min_x)), dy=0, color=color,
                          alpha=0.6,
                          length_includes_head=True,
                          head_length=max(0.03 * (segments_cordinate[s][1] - segments_cordinate[s][0]),0.006 * (max_x - min_x)),
                          head_width=0.2, linewidth=0.5,linestyle=linestyle)
                direction = 'left'
            else:
                if not fold_back_real_coordinate:
                    plt.arrow(x=x + 0.15 * (segments_cordinate[s][1] - segments_cordinate[s][0]), y=y + 0.25,
                              dx=0.3 * (segments_cordinate[s][1] - segments_cordinate[s][0]), dy=0, color='blue', alpha=0.6,
                              length_includes_head=True,
                              head_length=0.03 * (segments_cordinate[s][1] - segments_cordinate[s][0]),
                              head_width=0.2, linewidth=0.5)
                color = 'blue'
                linestyle = '-'
                if exp_foldback[s]['left'] == 0:
                    linestyle = '--'
                x_point = [segments_cordinate[s][0],
                           segments_cordinate[s][0] - max(0.1 * (segments_cordinate[s][1] - segments_cordinate[s][0]), 0.02 * (max_x - min_x)),
                           segments_cordinate[s][0] - max(0.1 * (segments_cordinate[s][1] - segments_cordinate[s][0]), 0.02 * (max_x - min_x))]
                arow_x_point = segments_cordinate[s][0] - max(0.1 * (segments_cordinate[s][1] - segments_cordinate[s][0]), 0.02 * (max_x - min_x))
                if s in foldback.keys() and fold_back_real_coordinate:
                    if len(foldback[s]['left']) != 0:
                        x_point = [foldback[s]['left'][0][1]+ max(0.1 * (segments_cordinate[s][1] - segments_cordinate[s][0]), 0.02 * (max_x - min_x)),
                                   foldback[s]['left'][0][1] ,
                                   foldback[s]['left'][0][2]]
                        arow_x_point = foldback[s]['left'][0][2]
                y_point = [y - 1 + 0.25, y - 1 + 0.25, y + 0.25]
                plt.plot(x_point, y_point, color=color, alpha=0.6, linewidth=1, linestyle = linestyle)
                plt.arrow(x=arow_x_point ,
                          y=y + 0.25,
                          dx=max(0.1 * (segments_cordinate[s][1] - segments_cordinate[s][0]), 0.02 * (max_x - min_x)), dy=0, color=color, alpha=0.6,
                          length_includes_head=True,
                          head_length=max(0.03 * (segments_cordinate[s][1] - segments_cordinate[s][0]), 0.006 * (max_x - min_x)),
                          head_width=0.2, linewidth=0.5, linestyle = linestyle)
                direction = 'right'
    y = y + 1
    for k in cn.keys():
        plt.annotate(cn[k], xy=((segments_cordinate[k][0] + segments_cordinate[k][1]) / 2 - 0.07 * (
                    segments_cordinate[k][1] - segments_cordinate[k][0]),
                                0 + 0.15), weight="bold")
    plt.ylim([0, max(y + 4, 10)])
    plt.title('score: ' + str(score[0])[:4] + ' cn: ' + str(score[1])[:4] + ' fn:' + str(score[2])[:4]+ ' sg:'+str(score[3])[:4]+' cosin_sim:'+ str(score[5])[:4]+' euclidean:'+str(score[6])[:4])
    plt.xlim([min_x - 0.105 * (max_x - min_x), max_x + 0.105 * (max_x - min_x)])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.yticks([])
    plt.savefig(output, dpi=300)
    plt.close()
