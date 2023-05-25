from collections import defaultdict
from os.path import exists

#################################### parsing RefAligner Copy Number file ##########################################
def parse_rmcap(cmap_dir):
    cov = {}
    cov = defaultdict(lambda: {}, cov)
    cop = {}
    cop = defaultdict(lambda: {}, cop)
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
                cov['chr' + chrom][pos] = cover
                cop['chr' + chrom][pos] = copynumber
    return cov, cop


#################################### parsing RefAligner smap file ##########################################
#######          Mode = 1 ---> BFB, 3----> detecting integration point
class BP:
    contig_id = ''
    direction1 = ''
    direction2 = ''
    pos1 = ''
    pos2 = ''
    chrom1 = ''
    chrom2 = ''
    line = ''
    type = ''

class SmapEntry:
    smap_id = ''
    q_id = ''
    ref_c_id1 = ''
    ref_c_id2 = ''
    ref_start = 0
    ref_end = 0
    query_start = 0
    query_end = 0
    confidence = 0
    xmap_id1 = ''
    xmap_id2 = ''
    sv_type = ''
    line = ''

def parse_smap(smap_dir):
    with open(smap_dir, 'r') as f:
        bfb_count = {}
        bfb_count = defaultdict(lambda: [], bfb_count)
        breakpoints = []
        translocations = []
        for line in f:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]
            if not line.startswith('#'):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                smap_entry = SmapEntry()
                smap_entry.ref_start = float(fD['RefStartPos'])
                smap_entry.ref_end = float(fD['RefEndPos'])
                smap_entry.xmap_id1 = int(fD['XmapID1'])
                smap_entry.xmap_id2 = int(fD['XmapID2'])
                smap_entry.q_id = int(fD['QryContigID'])
                smap_entry.ref_c_id1 = fD['RefcontigID1']
                smap_entry.ref_c_id2 = fD['RefcontigID2']
                smap_entry.smap_id = int(fD['SmapEntryID'])
                smap_entry.confidence = float(fD['Confidence'])
                smap_entry.query_start = float(fD['QryStartPos'])
                smap_entry.query_end = float(fD['QryEndPos'])
                smap_entry.sv_type = fD['Type']
                smap_entry.line = line
                breakpoints.append(smap_entry)
                if fD['Type'] == 'duplication_inverted':
                    chrom = 'chr' + fD['RefcontigID1']
                    # start_pos = float(fD['RefStartPos'])
                    # end_pos = float(fD['RefEndPos'])
                    # q_id = int(fD['QryContigID'])
                    # xmap_id1 = int(fD['XmapID1'])
                    # xmap_id2 = int(fD['XmapID2'])
                    # bfb_count[chrom].append([start_pos, end_pos, q_id, xmap_id1, xmap_id2])
                    bfb_count[chrom].append(smap_entry)
                # if mode == 2:
                    # line2 = line.strip().split('\t')
                    # chr1 = int(line2[2])
                    # chr2 = int(line2[3])
                    # pos1 = int(float(line2[6]))
                    # pos2 = int(float(line2[7]))
                    # sv_type = line2[9]
                    # item = [chr1, pos1, chr2, pos2, sv_type, line]
                    # breakpoints.append(item)

                # if mode == 3:
                if fD['Type'].startswith('trans'):
                        # t = BP()
                        # t.chrom1 = fD['RefcontigID1']
                        # t.chrom2 = fD['RefcontigID2']
                        # t.pos1 = float(fD['RefStartPos'])
                        # t.pos2 = float(fD['RefEndPos'])
                        # t.contig_id = int(fD['QryContigID'])
                        # t.line = line
                        # translocations.append(t)
                    translocations.append(smap_entry)
    # if mode == 1:
    #     return bfb_count
    # if mode == 2:
    #     return breakpoints
    # if mode == 3:
    #     return translocations
    return bfb_count, breakpoints, translocations


#################################### parsing Alignment xmap file ##########################################
def parse_xmap(xmapf):
    detailFields = ["XmapEntryID", "QryContigID", "RefContigID", "Orientation", "QryLen", "RefLen",
                    "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", "Alignment","Confidence"]
    numeric = ["QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos","Confidence"]
    xmapPair = {}
    with open(xmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                alnstring = ")" + fD["Alignment"] + "("
                # xmapAln[fD["XmapEntryID"]] = alnstring.rsplit(")(")[1:-1]

                xmapPair[fD["XmapEntryID"]] = {x: fD[x] for x in detailFields}
                for keyword in numeric:
                    xmapPair[fD["XmapEntryID"]][keyword] = float(xmapPair[fD["XmapEntryID"]][keyword])

    return xmapPair


#################################### parse centromer ##########################################
def parse_centro(centro):
    r = {}
    r = defaultdict(lambda: [], r)
    with open(centro, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            key = line[0]
            pos1 = int(line[1])
            pos2 = int(line[2])
            r[key].append(pos1)
            r[key].append(pos2)
    return r


####################### parse CNVkit ###########################################
def parse_cnvkit(cnvkit_dir):
    cop = {}
    cop = defaultdict(lambda: [], cop)
    with open(cnvkit_dir, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0] == 'X':
                chrom = 23
            elif line[0] == 'Y':
                chrom = 24
            else:
                chrom = int(line[0])
            start_pos = int(line[1])
            end_pos = int(line[2])
            cn = float(line[4])
            cop[chrom].append([start_pos, end_pos, cn])
            print([chrom, start_pos, end_pos, cn])
    return cop


########################## Parse AA graph file #########################################

def parse_AA_graph(graph_dir):
    discordant = []
    amplicon_region = []
    with open(graph_dir, 'r') as f:
        for line2 in f:
            if line2.startswith('discordant'):
                line = line2.strip().split('\t')
                bp = line[1]
                bp = bp.split('->')
                chr1 = int(bp[0].split(':')[0][3:])
                chr2 = int(bp[1].split(':')[0][3:])
                pos1 = int(bp[0].split(':')[1][:-1])
                pos2 = int(bp[1].split(':')[1][:-1])
                dir1 = bp[0].split(':')[1][-1]
                dir2 = bp[1].split(':')[1][-1]
                discordant.append([chr1, pos1, dir1, chr2, pos2, dir2, line2])
            if line2.startswith('sequence'):
                line = line2.strip().split('\t')
                start = line[1]
                end = line[2]
                chrom = int(start.split(':')[0][3:])
                start_pos = int(start.split(':')[1][:-1])
                end_pos = int(end.split(':')[1][:-1])
                amplicon_region.append([chrom, start_pos, end_pos, line2])
    return discordant, amplicon_region
#######################    Parse FaNDOM SV file ###########################
def parse_SV(SV_dir):
    breakpoints = []
    with open(SV_dir, 'r') as f:
        for line in f:
            if line[0].isdigit():
                line = line.strip().split('\t')
                chr1 = int(line[0])
                chr2 = int(line[3])
                pos1 = int(float(line[1]))
                pos2 = int(float(line[4]))
                dir1 = line[2]
                dir2 = line[5]
                item = [chr1, pos1, dir1, chr2, pos2, dir2, 'translocation', line[7]]
                breakpoints.append(item)
            elif line.startswith('dup') or line.startswith('inv'):
                line = line.strip().split('\t')
                chrom = int(line[1])
                pos1 = int(float(line[2]))
                pos2 = int(float(line[4]))
                dir1 = line[3]
                dir2 = line[5]
                item = [chrom, pos1, dir1, chrom, pos2, dir2, line[0], line[6]]
                breakpoints.append(item)
    return breakpoints
####################### Parse FaNDOM indel SV ##########################
def parse_indel(indel_dir,breakpoints):
    with open(indel_dir, 'r') as f:
        for line in f:
            if line.startswith('deletion') or line.startswith('insert'):
                line = line.strip().split('\t')
                chrom = int(line[1])
                pos1 = int(float(line[2]))
                pos2 = int(float(line[3]))
                dir1 = '+'
                dir2 = '+'
                item = [chrom, pos1, dir1, chrom, pos2, dir2, line[0], line[5]]
                breakpoints.append(item)
    return breakpoints
################################### Parse FaNDOM SV#################
def parse_fandom_sv(fandom_sv_dir):
    with open(fandom_sv_dir, 'r') as f:
        all = []
        for line in f:
            if line.startswith("#H"):
                head = line.rstrip().rsplit()[1:]
            if not line.startswith('#'):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                # if fD['Type'].startswith('Unknown'):
                t = BP()
                t.ref_c_id1 = fD['Chrom1']
                t.ref_c_id2 = fD['Chrom2']
                t.ref_start = float(fD['RefPos1'])
                t.ref_end = float(fD['RefPos2'])
                t.contig_id = fD['Ids']
                t.line = line
                t.direction1 = fD['Direction1']
                t.direction2 = fD['Directio2']
                t.type = fD['Type']
                all.append(t)
    return all
################## Parse molecule to contig alignment ###################
def generate_mol_contig(ids_list, folder_dir):
    contigs = {}
    contigs = defaultdict(lambda: {}, contigs)
    for id in ids_list:
        file_dir = folder_dir + '/EXP_REFINEFINAL1_noOutlier_contig' + str(id) + '.xmap'
        if not exists(file_dir):
            file_dir = folder_dir + '/EXP_REFINEFINAL1_contig' + str(id) + '.xmap'
            if not exists(file_dir):
                file_dir = folder_dir + '/exp_refineFinal1_contig' + str(id) + '.xmap'
        contigs[id] = defaultdict(lambda: [], contigs[id])
        with open(file_dir, 'r') as f:
            for line in f:
                if line.startswith("#h"):
                    head = line.rstrip().rsplit()[1:]
                if not line.startswith('#'):
                    fields = line.rstrip().rsplit()
                    fD = dict(zip(head, fields))
                    # if not line.startswith('#'):
                    # line = line.strip().split('\t')
                    molecule_id = int(fD['QryContigID'])
                    pair_alignment = fD['Alignment']
                    pair_alignment = pair_alignment.split('(')[1:]
                    for i in pair_alignment:
                        i = int(i.split(',')[0])
                        contigs[id][i].append(molecule_id)
    for key in contigs:
        for key2 in contigs[key]:
            contigs[key][key2] = list(set(contigs[key][key2]))
    return contigs


########################return molecule labels in alignment pair like (a,b)(c,d)(e,f) return [b ,d f]#########################
def query_index(aln):
    aln = aln.split(')')[:-1]
    ans = []
    for i in aln:
        ans.append(int(i.split(',')[1]))
    return ans

################### Return molecules supports intersection of two alignments################
def molecule_support_bp(xmap_id1, xmap_id2, contigs, xmap):
    xmap_id2 = str(xmap_id2)
    xmap_id1 = str(xmap_id1)
    q_id = int(xmap[xmap_id1]['QryContigID'])######### baraye jadide int bardashte shod
    aln1 = query_index(xmap[xmap_id1]['Alignment'])
    aln2 = query_index(xmap[xmap_id2]['Alignment'])
    aln1_mol = []
    for i in aln1:
        aln1_mol.extend(contigs[q_id][i])
    aln2_mol = []
    for i in aln2:
        aln2_mol.extend(contigs[q_id][i])
    aln1_mol = set(aln1_mol)
    aln2_mol = set(aln2_mol)
    return list(aln1_mol.intersection(aln2_mol))


################################ Parse bed file ###########################
def parse_bed(bed_dir):
    l = []
    with open(bed_dir, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            chrom = int(line[0][3:])
            start = int(float(line[1]))
            end = int(float(line[2]))
            l.append([chrom, start, end])
    return l 
################### return contig occ ##########
def parse_occ(cmap_dir):
    d = defaultdict(lambda: defaultdict(float))
    with open(cmap_dir,'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                contig_id = line[0]
                contig_label = line[3]
                occurence = float(line[8])
                d[contig_id][contig_label] = occurence
    return d
##################### parse segmentation output from R #############
def parse_segmentation(file_dir):
    ans = []
    with open(file_dir,'r') as f :
        for line in f:
            if not line.startswith('"output.ID"'):
                line = line.strip().split('\t')
                CN = float(line[-1])
                start = int(line[3])
                end = int(line[4])
                ans.append([start,end,CN])
    return ans
