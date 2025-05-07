## Overview

This repository provides a learning implementation of an Overlap­‐Layout­‐Consensus (OLC) assembler for long reads, inspired by string-graph methods in Miniasm (https://github.com/lh3/miniasm) (https://doi.org/10.1093/bioinformatics/btw152) and the foundational string-graph formulation by Myers (https://doi.org/10.1093/bioinformatics/bti1114). It is **not** intended for production use but rather to illustrate the core data structures and algorithms behind overlap-graph assemblers.
## Prerequisites
- Python 3.7+
- A set of long reads in FASTA/FASTQ
- Minimap2 to compute all-vs-all read overlaps in PAF format

## Getting Started

1. **Compute pairwise overlaps**
    
    - `minimap2 -x ava-ont reads.fq > overlaps.paf`
    
    - This produces a PAF (Pairwise mApping Format) file: a tab-delimited list of read-to-read alignments, where columns 1–12 describe query/target names, positions, strand, match counts, alignment lengths, and mapping quality (https://lh3.github.io/minimap2/minimap2.html).
    
2. **Filter overlaps with `Preprocessing.py`**    
    - Removes self-hits, low-coverage reads (< 3), and short overlaps (< 2 kb)
    
3. **Build and simplify overlap graph with `Assembly.py`**
    - Parses overlaps into a bidirected string graph, then applies tip trimming, transitive-edge reduction, and bubble popping.
## Methodology and Algorithms

### 1. PAF Filtering

- **Self‐alignment removal**  
    Discards any alignment where a read maps to itself. This prevents trivial loops in the graph.
    
- **Coverage‐based read filtering**  
    Reads with fewer than _k_ overlaps (default _k_ = 3) are removed, eliminating low-coverage fragments unlikely to contribute to contigs or caused by sequencing errors.
    
- **Length thresholding**  
    Alignments shorter than a minimum (default 2 000 bp) are discarded to focus on high-confidence overlaps and reduce noise.
    
### 2. Overlap Classification & Graph Construction

Each alignment contains a query read and a target read. Overlaps are very informative to build a de novo assembly, but not every alignment is an overlap, so each remaining alignment is classified into one of four types by comparing:

- **b₁, e₁** = start/end position of the alignment in the query read; **b₂, e₂** = start/end position of the alignment in the target read (or its reverse complement); **ℓ₁, ℓ₂** = length of query and target reads
- **Overhang** = min(b₁,b₂) + min(ℓ₁−e₁, ℓ₂−e₂) is the part next to the overlap where the reads don't align, but should in case of perfect overlap
- **maplen** = max(e₁−b₁, e₂−b₂) is the effective overlap length
#### a. Internal Matches

If `overhang > min(max_overhang, maplen * overhang_ratio)`, the alignment is deemed an internal match (reads overlap mostly in the interior) and is discarded.

#### b. Containment

- Query contained in target if `b₁ ≤ b₂  and  (ℓ₁−e₁) ≤ (ℓ₂−e₂)`
    
- Target contained in query if `b₁ ≥ b₂  and  (ℓ₁−e₁) ≥ (ℓ₂−e₂)`
    
Contained reads are reads that are entirely found in another read and add no new sequence to the graph, so they are discarded.
#### c. Proper Overlaps

Remaining alignments are proper overlaps. For each proper overlap, two orientations of both reads are added as nodes in a bidirected graph: read + and read – (reverse complement). Directed edges link one node to another, weighted by the length of the non-overlapping tail. 

The nodes contain the read names and orientation, the edges contain the distance between overlapping reads. This is done in two directions because an overlap works both ways (query → target and target → query). This means that every node has a reverse complement node where the incoming edges of the original node correspond to the outgoing edges of the reverse complement node.

### 3. Graph Cleaning

To improve robustness to sequencing errors and repeats, the graph is simplified via four standard operations:

#### a. Tip Trimming

Short “dead-end” paths (tips) that diverge from main contigs, often caused by chimeric reads or errors, are removed. Paths of 4 or less nodes where the out-degree = 1, preceded by a node where the out-degree > 1 are removed. This is done using the algorithm described in https://doi.org/10.1093/bioinformatics/bti1114.

#### b. Transitive Edge Reduction

If read A→B and B→C suffice to explain A→C, then the direct A→C edge is called a transitive edge. It is redundant and removed. The collapse of redundant connections yields a minimal string graph.  This is done using the algorithm described in https://doi.org/10.1093/bioinformatics/bti1114.

#### c. Bubble Popping

A bubble arises when two or more alternative paths diverge from a node and reconverge, often caused by errors, heterozygosity, or repeats. Bubbles are detected with a bounded DFS search. If all incoming edges of a node are explored, it gets added to a stack. Nodes get popped of the stack and their neighbours get explored. This is done using a slightly modified version of Algorithm 6 in https://doi.org/10.1093/bioinformatics/btw152). After bubbles are found, they are popped by removing all but the best (highest depth) path through each bubble.

#### d. Heuristic edge removal

The above methods aren't enough to remove low quality edges, so an additional heuristics-based cleanup step is needed. Very long edges (indicating small overlap length) are removed from nodes with multiple edges. This can simplify complex structures in the graph by removing weaker (shorter) overlaps. 

### 4. Contig reconstruction

After cleaning, contigs are generated by simple path traversal:

1. Identify starting nodes: Nodes with in-degree ≠ 1 or out-degree ≠ 1 mark path ends.
    
2. Traverse non-branching paths: Follow edges through nodes where in-degree = out-degree = 1, concatenating read sequences minus overlap to assemble contiguous sequences.
    
3. Report contigs: Output each assembled path as a FASTA record, including sequence and path metadata (length, read IDs).
    
This greedy traversal yields contigs representing maximal unambiguous sequences in the overlap graph.
## Further Reading

- Miniasm: Li H. Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences. Bioinformatics. 2016 Jul 15;32(14):2103-10. doi: https://doi.org/10.1093/bioinformatics/btw152. Epub 2016 Mar 19. PMID: 27153593; PMCID: PMC4937194.
    
- String Graphs: Eugene W. Myers, The fragment assembly string graph, _Bioinformatics_, Volume 21, Issue suppl_2, September 2005, Pages ii79–ii85, [https://doi.org/10.1093/bioinformatics/bti1114](https://doi.org/10.1093/bioinformatics/bti1114)

- Overlap graphs and de Bruijn graphs: Rizzi, R., Beretta, S., Patterson, M. _et al._ Overlap graphs and _de Bruijn_ graphs: data structures for _de novo_ genome assembly in the big data era. _Quant Biol_ **7**, 278–292 (2019). https://doi.org/10.1007/s40484-019-0181-x

- Visual overview: https://bio.libretexts.org/Bookshelves/Computational_Biology/Book%3A_Computational_Biology_-_Genomes_Networks_and_Evolution_(Kellis_et_al.)/05%3A_Genome_Assembly_and_Whole-Genome_Alignment/5.03%3A_Genome_Assembly_II-_String_graph_methods

## License

This code is released under the GPL-3.0 license ([https://www.gnu.org/licenses/quick-guide-gplv3.html](https://www.gnu.org/licenses/quick-guide-gplv3.html))