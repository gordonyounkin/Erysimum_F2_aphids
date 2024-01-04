import sys

print("This script parses genetic map created by MSTmap and reformats for ALLMAPS."
      "Requires two inputs:"
      "    1. Genetic map created in MSTmap"
      "    2. Filename for reformated output output")

map1 = open(sys.argv[1], 'r')
out1 = open(sys.argv[2], 'w')
out1.write('{},{},{},{}\n'.format('Scaffold ID', 'scaffold position', 'LG', 'genetic position'))

lg = 0
line1 = map1.readline()
for line1 in map1:
    if line1.startswith(';BEGINOFGROUP'):
        lg += 1
        line1 = map1.readline()
        while not line1.startswith(';ENDOFGROUP'):
            contig_pos, dist = line1.strip().split('\t')
            contig, pos = contig_pos.split('_')[-2:]
            out1.write('{},{},{},{}\n'.format(contig, pos, lg, dist))
            line1 = map1.readline()
    else:
        line1 = map1.readline()
