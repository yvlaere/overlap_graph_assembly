# paths
paf_file = ''
reads_fasta_file = ''
filtered_paf = ''

# parameters
cov_len = 2000
min_overlap = 2000
cov_cutoff = 3

# load reads in memory
read_dict = {}
cov = {}

with open(reads_fasta_file, 'r') as reads:
    lines = iter(reads)
    for key in lines:
        value = next(lines, None)
        key = key.split()[0].strip()[1:]
        value = value.strip()
        read_dict[key] = value
        cov[key] = 0
print("reads loaded in memory")

nbases = {}

with open(paf_file, 'r') as alignments:
    for line in alignments:
        # parse paf
        qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq, *rest = line.split('\t')
        qlen, qstart, qend, tlen, tstart, tend, nmatch, alen, mapq = int(qlen), int(qstart), int(qend), int(tlen), int(tstart), int(tend), int(nmatch), int(alen), int(mapq)

        # avoid self-mappings
        if qname != tname:
            # increase coverage if mapping > 2000 bp
            if ((qend - qstart >= cov_len) and (tend - tstart >= cov_len)):
                cov[qname] += 1
                cov[tname] += 1

            # count number of matched bases
            alignment_name = " ".join(sorted([qname, tname]))
            if alignment_name in nbases.keys():
                if nmatch > nbases[alignment_name]:
                    nbases[alignment_name] = nmatch
            else:
                nbases[alignment_name] = nmatch

print("Calculated coverage")

with open(paf_file, 'r') as alignments:
    for line in alignments:
        # parse paf
        qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq, *rest = line.split('\t')
        qlen, qstart, qend, tlen, tstart, tend, nmatch, alen, mapq = int(qlen), int(qstart), int(qend), int(tlen), int(tstart), int(tend), int(nmatch), int(alen), int(mapq)

        # avoid self-mapping
        if qname != tname:
            # check if we have the best overlap between the current reads
            alignment_name = " ".join(sorted([qname, tname]))
            if nbases[alignment_name] == nmatch:
                # check if overlap length is high enough
                if ((qend - qstart >= min_overlap) and (tend - tstart >= min_overlap)):
                    # check if coverage is high enough for the involved reads
                    if ((cov[qname] >= cov_cutoff) and (cov[tname] >= cov_cutoff)):
                        with open(filtered_paf, 'a') as filtered_alignments:
                            filtered_alignments.write(line)
                        filtered_alignments.close()
