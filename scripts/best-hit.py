import sys

d = {}

for line in open(sys.argv[1]):
    if line.startswith("Download"):
        continue
    elif line.startswith("query"):
        print line,
    else:
        data = line.rstrip().split('\t')
        query = '_'.join(data[0].split('|')[1].split('_')[0:6])
        data[0] = query
        hit = data[1]
        if d.has_key(query):
            continue
        else:
            d[query] = hit
            print '\t'.join(data) 
