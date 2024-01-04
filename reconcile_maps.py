import sys
from scipy.stats import linregress

print("This script reconciles Hi-C and genetic map. Required arguments:"
      "    1. Genetic map in csv format, (MSTmap output reformatted by reformat_map.py; "
      "       likely KonxElb_SPLITCONTIGS_binned_map.csv)"
      "    2. Hi-C map in original bed format (PGA_assembly.fasta.bed.bed)"
      "    3. Name (and relative path) of output file, without extension--will make bed, agp, and map with both genetic and physical positions")

# load genetic map
gm = open(sys.argv[1], 'r')
# load Hi-C map file
hic = open(sys.argv[2], 'r')
# make bed and agp file for reconciled output
bedout = open("{}.bed".format(sys.argv[3]), 'w')
agpout = open("{}.agp".format(sys.argv[3]), 'w')
mapout = open("{}_marker_physical_pos.tsv".format(sys.argv[3]), 'w')
mapout.write("{}\t{}\t{}\t{}\n".format("marker", "linkage_group", "genetic_pos_cM", "physical_pos_bp"))

# parse genetic map
colnames = gm.readline().strip().split(',')
line1 = gm.readline()
gmdict = {}
while line1:
    scaf, pos, lg, cM = line1.strip().split(',')
    if lg in gmdict:
        if scaf in gmdict[lg]:
            gmdict[lg][scaf]["pos"] += [float(pos)]
            gmdict[lg][scaf]["cM"] += [float(cM)]
        else:
            gmdict[lg][scaf] = {
                "pos": [float(pos)],
                "cM": [float(cM)]
            }
    else:
        gmdict[lg] = {
            scaf: {
                "pos": [float(pos)],
                "cM": [float(cM)]
            }
        }

    line1 = gm.readline()

# for each scaffold, get orientation in genetic map
for lg in gmdict:
    for scaf in gmdict[lg]:
        if len(gmdict[lg][scaf]["pos"]) == 1:
            gmdict[lg][scaf]["orientation"] = ''
        else:
            slope, intercept, r_value, p_value, std_err = \
                linregress(gmdict[lg][scaf]["pos"], gmdict[lg][scaf]["cM"])
            if slope > 0:
                gmdict[lg][scaf]["orientation"] = '+'
            elif slope < 0:
                gmdict[lg][scaf]["orientation"] = '-'
            else:
                gmdict[lg][scaf]["orientation"] = ''

# parse hic map
line2 = hic.readline()
hicdict = {}
while line2:
    line2split = line2.strip().split('\t')
    if len(line2split) == 4:
        next
    else:
        lg, start, end, scaf, v1, orientation = line2split
        # rename contigs to match genetic map file
        scafname = "S{}".format(scaf.strip('__unscaffolded').strip('|arrow|pilon'))
    # Contig S000001F was split into 1.1 and 1.2--only 1.1 is in Hi-C map #
    if scafname == "S000001F":
        scafname = "S000001.1F"
        start = 0
        end = 3053566
    if lg in hicdict:
        hicdict[lg][scafname] = {
            "pos": int(start),
            "end": int(end),
            "orientation": orientation
        }
    else:
        hicdict[lg] = {
            scafname: {
                "pos": int(start),
                "end": int(end),
                "orientation": orientation
            }
        }

    line2 = hic.readline()

gap = 0
# Loop through linkage groups of genetic map, generate lists:
# 1. Contig names in order
# 2. Contig orientations
for lg_to_parse in gmdict:
    curr_lg = gmdict[lg_to_parse]
    curr_lg_cM = [min(curr_lg[x]["cM"]) for x in curr_lg]
    curr_lg_scaf = [x for x in curr_lg]
    curr_lg_orient = [curr_lg[x]["orientation"] for x in curr_lg]
    # get correct order based on location in cM
    sort = sorted(range(len(curr_lg_cM)), key=lambda k: curr_lg_cM[k])
    gm_cM = [curr_lg_cM[x] for x in sort]
    gm_scaf = [curr_lg_scaf[x] for x in sort]
    gm_orient = [curr_lg_orient[x] for x in sort]

    # get corresponding Hi-C chromosome
    for chrom in hicdict:
        if gm_scaf[0] in hicdict[chrom]:
            curr_chrom = chrom
    hic_scaf_name = curr_chrom.split('_')[1]
    # make same lists: contig names in order, contig orientations
    hic_pos_unsorted = [hicdict[curr_chrom][x]["pos"] for x in hicdict[curr_chrom]]
    hic_scaf_unsorted = [x for x in hicdict[curr_chrom]]
    hic_orient_unsorted = [hicdict[curr_chrom][x]["orientation"] for x in hicdict[curr_chrom]]
    # get correct order based on hic_pos
    sort = sorted(range(len(hic_pos_unsorted)), key=lambda k: hic_pos_unsorted[k])
    hic_pos = [hic_pos_unsorted[x] for x in sort]
    hic_scaf = [hic_scaf_unsorted[x] for x in sort]
    hic_orient = [hic_orient_unsorted[x] for x in sort]

    ### Start reconciling maps ###
    # 1. Get indices of Genetic map contigs in Hi-C map
    gm_in_hic = []
    for x in gm_scaf:
        if x in hic_scaf:
            gm_in_hic += [hic_scaf.index(x)]
        else:
            gm_in_hic += [-9999]

    # 2. Give each contig NOT in genetic map a 'parent' scaffold to follow in final assembly
    #    Creating dictionary with genetic map scaffolds as keys, with each entry containing 2 lists:
    #    contigs = ordered list of contigs associated with genetic map scaffold
    #    orientations = ordered list with orientation of contigs
    hic_parent_scafs = {}
    for i,x in enumerate(hic_scaf):
        if x in gm_scaf:
            closest_value = i
            parent_scaf = x
        else:
            closest_value = min(gm_in_hic, key=lambda list_value: abs(list_value - i))
            parent_scaf = hic_scaf[closest_value]
        if parent_scaf in hic_parent_scafs:
            hic_parent_scafs[parent_scaf]["contigs"] += [x]
            hic_parent_scafs[parent_scaf]["orientations"] += [hic_orient[i]]
        else:
            hic_parent_scafs[parent_scaf] = {
                "contigs": [x],
                "orientations": [hic_orient[i]]
            }

    # 3. Assemble contig based on order and orientation in genetic map
    merged_contigs = []
    merged_orientations = []
    gm_hic_agree = True
    for i,curr_scaf in enumerate(gm_scaf):
        orient_in_gm = gm_orient[i]
        # S000001.2F isn't in genetic map and thus won't be in hic_parent_scaf--just add in place to merged contig
        if curr_scaf not in hic_parent_scafs:
            merged_contigs += [curr_scaf]
            merged_orientations += [orient_in_gm]
            continue
        orient_in_hic = hic_parent_scafs[curr_scaf]["orientations"][hic_parent_scafs[curr_scaf]["contigs"].index(curr_scaf)]
        # for contigs with only one marker in genetic map, use orientation from HiC map, (or opposite, if previous contig was in opposite orientation #
        if orient_in_gm == '':
            if gm_hic_agree:
                orient_in_gm = orient_in_hic
            else:
                orient_in_gm = '+' if orient_in_hic == '-' else '-'
        if orient_in_gm == orient_in_hic:
            gm_hic_agree = True
            merged_contigs += hic_parent_scafs[curr_scaf]["contigs"]
            merged_orientations += hic_parent_scafs[curr_scaf]["orientations"]
        else:
            gm_hic_agree = False
            merged_contigs += hic_parent_scafs[curr_scaf]["contigs"][::-1]
            temp = hic_parent_scafs[curr_scaf]["orientations"][::-1]
            merged_orientations += ['+' if x=='-' else '-' for x in temp]

    # 4. Write to output .bed and agp file
    #    in tab-separated format: LGname, contig_start, contig_end, contig_name, 25, contig_orientation
    #    Each contig separated by 100 Ns
    contig_start = 0
    agp_id = 0
    for i,contig in enumerate(merged_contigs):
        # write contig to bed file
        orientation = merged_orientations[i]
        if contig == 'S000001.2F':
            old_contig_start = 0
            old_contig_end = 1549008
        else:
            old_contig_start = hicdict[curr_chrom][contig]["pos"]
            old_contig_end = hicdict[curr_chrom][contig]["end"]
        contig_length = old_contig_end - old_contig_start
        new_contig_end = contig_start + contig_length
        bedout.write("chr{}\t{}\t{}\t{}\t{}\t{}\n".format(lg_to_parse, contig_start, new_contig_end, contig, 1, orientation))
        agp_id += 1
        agpout.write("chr{}\t{}\t{}\t{}\tW\t{}\t1\t{}\t{}\n".format(lg_to_parse, contig_start+1, new_contig_end, agp_id,
                                                                    contig, contig_length, orientation))
        # write genetic and physical positions of markers to file
        if contig in gmdict[lg_to_parse]:
            for j, marker_pos in enumerate(gmdict[lg_to_parse][contig]["pos"]):
                if orientation == '+':
                    physical_pos = int(contig_start + marker_pos + 1)
                else:
                    physical_pos = int(contig_start + (contig_length - marker_pos))
                genetic_pos = gmdict[lg_to_parse][contig]["cM"][j]
                mapout.write("bins_{}_{}\t{}\t{}\t{}\n".format(contig, int(marker_pos), lg_to_parse, genetic_pos, physical_pos))

        contig_start = new_contig_end
        # write gap to bed file (unless last contig on chromosome)
        if i == len(merged_contigs) - 1:
            continue
        else:
            gap += 1
            new_contig_end = contig_start + 100
            bedout.write("chr{}\t{}\t{}\tgap_{}\t{}\t{}\n".format(lg_to_parse, contig_start, new_contig_end, gap,'',''))
            agp_id += 1
            agpout.write("chr{}\t{}\t{}\t{}\tU\t100\tmap\tyes\t.\n".format(lg_to_parse, contig_start+1, new_contig_end,
                                                                           agp_id))
            contig_start = new_contig_end
    # also print names of contigs associated with scaffold but not ordered #
    unordered_id = "unordered_in_{}".format(hic_scaf_name)
    for hic_contig in hicdict:
        if unordered_id in hic_contig:
            contig_id = "S{}".format(hic_contig.split('|')[0])
            lg_unordered_id = "{}_unordered_in_chr{}".format(contig_id, lg_to_parse)
            old_contig_start = hicdict[hic_contig][contig_id]["pos"]
            old_contig_end = hicdict[hic_contig][contig_id]["end"]
            bedout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(lg_unordered_id, old_contig_start, old_contig_end, contig_id, 1, "+"))
            agpout.write("chr{}\t{}\t{}\t1\tW\t{}\t1\t{}\t+\n".format(lg_unordered_id, old_contig_start+1, old_contig_end,
                                                                        contig_id, old_contig_end-old_contig_start))

### Finally, add unscaffolded contigs to bed file
for hic_contig in hicdict:
    if "unscaffolded" in hic_contig:
        contig_id = "S{}".format(hic_contig.split('|')[0])
        lg_unscaffolded_id = "{}_unscaffolded".format(contig_id)
        old_contig_start = hicdict[hic_contig][contig_id]["pos"]
        old_contig_end = hicdict[hic_contig][contig_id]["end"]
        bedout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(lg_unscaffolded_id, old_contig_start, old_contig_end, contig_id, 25,
                                                     "+"))
        agpout.write("chr{}\t{}\t{}\t1\tW\t{}\t1\t{}\t+\n".format(lg_unscaffolded_id, old_contig_start + 1, old_contig_end,
                                                                  contig_id, old_contig_end-old_contig_start))

gm.close()
hic.close()
bedout.close()
agpout.close()
sys.exit()