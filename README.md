# overlap_graph_assembly
Tool for de novo genome assembly of long reads using an overlap graph

# Manual
1. Perform self-alignment of reads to create a paf file
   - Example using minimap2 for ONT reads:
   - `minimap2 -x ava-ont ont-reads.fq > reads.paf`
2. Filter aligned reads with Preprocessing.py
3. Create a fasta file for the assembled genome with Assembly.py

# Methodology
## Filtering PAF file
As a first step, the paf file is filtered to remove unwanted alignments. Self-alignments, an alignment between a read and itself, are removed. All alignments associated with reads that have a low coverage (counted as number of alignments, default cutoff is 3) are also removed. Finally, alignments that are too short (2000 bp by default) are removed.

## Creating the overlap graph
The filtered paf file is parsed and all alignments are classified as an internal match (the overhang next to the overlap is too large), contained (one read entirely contains the other read), or an overlap. The overlaps are used to create an overlap graph. Each read is represented by two nodes, one for its normal orientation and one for the reverse complement. Nodes can have directed edges, whose weight represents the length of the nonoverlapping part of the read. To simplify the graph, tips are trimmed, transitive edges are reduced and bubbles are popped.


# Inspiration and further explanation
https://doi.org/10.1093/bioinformatics/btw152

https://doi.org/10.1093/bioinformatics/bti1114
