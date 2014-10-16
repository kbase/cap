import sys

d = {}

for line in open(sys.argv[1]):
    data = line.rstrip().split('\t')
    query = data[0]
    hit = data[1]
    if d.has_key(query):
        continue
    else:
        d[query] = hit
        print line, 
