class Node:
    """
    Class for nodes in the graph. Each read is represented by two nodes, one for the normal orientation and one for the reverse complement.
    A node can have directed edges pointing to other nodes.
    """

    def __init__(self, node_id):
        """
        Initializes a new node with the given node_id.
        :param node_id: the read name plus information on its orientation (+/- appended to the read name)
        """
        self.node_id = node_id
        # the directed edges are stored in a dictionary, each key is a node_id and each value is the length of the edge
        # the length of the edge is the length of the part of the first read that doesn't overlap with the second
        self.edges = {}

    def add_edge(self, target_node, edge_len):
        """
        Adds a directed edge to the node.
        :param target_node: node_id of the target node
        :param edge_len: length of the directed edge, length of the non-overlapping part of the first read
        :return:
        """
        assert target_node not in self.edges, "Cannot add more than one edge to the same target node"
        self.edges[target_node] = edge_len

    def remove_edge(self, target_node):
        """
        Removes a directed edge from the node.
        :param target_node: node_id of the target node
        :return:
        """
        assert target_node in self.edges, "Can't remove edge that doesn't exist"
        self.edges.pop(target_node)

    def sort_edges(self):
        """
        Sorts the directed edges dictionary from shortest ot longest according to the edge length.
        :return:
        """
        self.edges = dict(sorted(self.edges.items(), key=lambda x: x[1]))

class AssemblyGraph:
    """
    Class for the assembly graph. The graph is a dictionary of Node objects
    """
    def __init__(self):
        """
        Initializes a new empty graph.
        """
        self.nodes = {}

    def add_node(self, node_id):
        """
        Adds a new node to the graph if it is not already present.
        :param node_id: the read name plus information on its orientation (+/- appended to the read name)
        :return:
        """
        # add node if it doesn't already exist, if it already exists, don't do anything
        if node_id not in self.nodes:
            self.nodes[node_id] = Node(node_id)

    def add_edge(self, from_id, to_id, edge_len):
        """
        Adds a directed edge between the from_id and to_id nodes in the graph.
        :param from_id: node_id of the source node
        :param to_id: node_id of the target node
        :param edge_len: length of the directed edge, length of the non-overlapping part of the first read
        :return:
        """
        assert from_id in self.nodes and to_id in self.nodes, "from_id and/or to_id aren't present in the graph"
        self.nodes[from_id].add_edge(to_id, edge_len)

    def remove_node(self, node_id):
        """
        Removes a node from the graph.
        :param node_id: node_id of the node to remove
        :return:
        """
        # remove node if it exists, otherwise do nothing
        if node_id in self.nodes:
            self.nodes.pop(node_id)

def load_reads(fasta_file):
    """
    Load reads in a dictionary
    :param fasta_file: file path to fasta file containing reads
    :return:
    """

    read_dict = {}
    with open(reads_fasta_file, 'r') as reads:
        lines = iter(reads)
        for key in lines:
            value = next(lines, None)
            key = key.split()[0].strip()[1:]
            value = value.strip()
            read_dict[key] = value
    return read_dict

def create_string_graph(paf_file, max_overhang, overhang_ratio):
    """
    Create the string graph from filtered all-to-all alignment.
    :param paf_file: filtered all-to-all alignment file in paf format
    :param max_overhang: maximum overhang length to consider an alignment an overlap
    :param overhang_ratio: ratio of read length to consider an alignment an overlap
    :return:
    """

    # initialize
    nr_internal_match = 0
    nr_first_contained = 0
    nr_second_contained = 0
    nr_proper_overlaps = 0
    g = AssemblyGraph()

    # classify alignments
    with open(paf_file, 'r') as alignments:
        for line in alignments:
            # parse alignment
            qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq, *rest = line.split('\t')
            qlen, qstart, qend, tlen, tstart, tend, nmatch, alen, mapq = int(qlen), int(qstart), int(qend), int(tlen), int(tstart), int(tend), int(nmatch), int(alen), int(mapq)

            # define beginning and end of overlap based on read orientation
            # naming corresponds to that in https://doi.org/10.1093/bioinformatics/btw152
            b1, e1, l1 = qstart, qend, qlen
            if strand == "+":
                b2, e2, l2 = tstart, tend, tlen
            else:
                b2, e2, l2 = tlen - tend, tlen - tstart, tlen

            # overhang is the part next to the overlap where the reads don't align, but should in case of perfect overlap
            # it corresponds with the grey region in Figure 1 in https://doi.org/10.1093/bioinformatics/btw152
            overhang = min(b1, b2) + min(l1 - e1, l2 - e2)

            # longest overlap length (could be from read a or read b)
            maplen = max(e1 - b1, e2 - b2)

            # do multiple checks to see if this alignment can be considered an overlap
            # overlaps are a subset of alignments where (in theory) two read edges, one from each read, are part of the alignment
            if overhang > min(max_overhang, maplen*overhang_ratio):
                # overhang is too big, so it's considered an internal match
                nr_internal_match += 1
            elif ((b1 <= b2) and ((l1 - e1) <= (l2 - e2))):
                # first read is contained in the second
                nr_first_contained += 1
            elif ((b1 >= b2) and ((l1 - e1) >= (l2 - e2))):
                # second read is contained in the first
                nr_second_contained += 1
            elif b1 > b2:
                # first to second overlap
                nr_proper_overlaps += 1

                # add nodes and edges, once for each orientation of the reads
                # default orientation
                g.add_node(qname + "+")
                g.add_node(tname + strand)
                edge1_len = b1 - b2
                g.add_edge(qname + "+", tname + strand, edge1_len)

                # reverse complement
                g.add_node(qname + "-")
                rc_strand = "-" if strand == "+" else "+"
                g.add_node(tname + rc_strand)
                edge2_len = (l2 - e2) - (l1 - e1)
                g.add_edge(tname + rc_strand, qname + "-", edge2_len)
            else:
                # second to first overlap
                nr_proper_overlaps += 1

                # add nodes and edges, once for each orientation of the reads
                # default orientation
                g.add_node(qname + "+")
                g.add_node(tname + strand)
                edge1_len = b2 - b1
                g.add_edge(tname + strand, qname + "+", edge1_len)

                # reverse complement
                g.add_node(qname + "-")
                rc_strand = "-" if strand == "+" else "+"
                g.add_node(tname + rc_strand)
                edge2_len = (l1 - e1) - (l2 - e2)
                g.add_edge(qname + "-", tname + rc_strand, edge2_len)

    return g, nr_internal_match, nr_first_contained, nr_second_contained, nr_proper_overlaps

def rc_node(n):
    """
    Create a node_id for the reverse complement of a node.
    :param n: node_id of the node of which the reverse complement is wanted
    :return:
    """

    rc_strand = "-" if n[-1] == "+" else "+"
    return n[:-1] + rc_strand

def check_synchronization(g):
    """
    Check if the bigraph is synchronized:
        1. Every node has a reverse complement.
        2. Ingoing edges of every node correspond to outgoing edges of its reverse complement.
    :param g: graph object
    :return:
    """

    for n in g.nodes:
        # check if reverse complement node exists
        n_rc = rc_node(n)
        assert n_rc in g.nodes, "Reverse complement not found for " + str(n) + ". The bigraph is not synchronized."

        # check if every edge has a counterpart
        for t in g.nodes[n].edges:
            t_rc = rc_node(t)
            assert n_rc in g.nodes[t_rc].edges, "Corresponding edge not found for " + str(n) + ". The bigraph is not synchronized."

def reduce_transitive_edges(g, fuzz):
    """
    Reduce transitive edges. transitive edges are redundant edges that don't add any information to the graph.
    Say read 1 overlaps with read 2 and read 2 overlaps with read 3 and read 1 also overlaps with read 3, then this last overlap is redundant, represented by a transitive edge.
    Algorithm based on https://doi.org/10.1093/bioinformatics/bti1114
    :param g: graph object
    :param fuzz: integer to take overlap length discrepancies into account
    :return:
    """

    # initialize
    mark = {}
    reduce = {}
    nr_reduced = 0

    # mark all nodes as vacant (not in play) and all edges as not to be reduced
    for n1 in g.nodes:
        mark[n1] = "vacant"
        reduce[n1] = {}

        # sort edges, the algorithm assumes all edges are sorted in increasing order
        g.nodes[n1].sort_edges()

        for n2 in g.nodes[n1].edges.keys():
            reduce[n1][n2] = False

    # Mark transitive edges
    # For every node compare nodes encountered two steps into the future with those encountered one step into the future
    for n1 in g.nodes:

        # only calculate for nodes that have outgoing edges
        if len(g.nodes[n1].edges) > 0:

            # mark all edge target nodes of n1 as in play
            for n2 in g.nodes[n1].edges.keys():
                mark[n2] = "inplay"

            # get the longest edge of n1
            longest = list(g.nodes[n1].edges.values())[-1] + fuzz

            for n2 in g.nodes[n1].edges.keys():
                if mark[n2] == "inplay":
                    for n3 in g.nodes[n2].edges.keys():
                        if g.nodes[n2].edges[n3] + g.nodes[n1].edges[n2] <= longest:
                            # check if any of the edge target nodes of n2 are in play (meaning they are also an edge target node of n1)
                            # if the path from n1 to n2 to n3 is smaller than the longest edge out of n1, the n1 to n3 edge gets reduced
                            if mark[n3] == "inplay":
                                mark[n3] = "eliminated"

            for n2 in g.nodes[n1].edges.keys():
                for n3 in g.nodes[n2].edges.keys():
                    if ((g.nodes[n2].edges[n3] < fuzz) or (g.nodes[n2].edges[n3] ==  min(g.nodes[n2].edges.values()))):
                        # check if any of the edge target nodes of n2 are in play (meaning they are also an edge target node of n1)
                        # if they are in play and the n2 to n3 edge is very small or the smallest of the outgoing n2 edges, the n1 to n3 edge gets reduced
                        if mark[n3] == "inplay":
                            mark[n3] = "eliminated"

            # mark edges going to eliminated nodes as reduced
            for n2 in g.nodes[n1].edges.keys():
                if mark[n2] == "eliminated":
                    reduce[n1][n2] = True
                mark[n2] = "vacant"

    # reduce transitive edges
    for n1 in g.nodes:
        edges = list(g.nodes[n1].edges.keys())
        for n2 in edges:
            if reduce[n1][n2]:
                # remove edge
                g.nodes[n1].remove_edge(n2)
                nr_reduced += 1

                # remove reverse complement edge if it exists
                n1_rc = rc_node(n1)
                n2_rc = rc_node(n2)
                if n1_rc in g.nodes[n2_rc].edges:
                    g.nodes[n2_rc].remove_edge(n1_rc)
                    nr_reduced += 1

    return g, nr_reduced

def detect_bubbles(g, max_bubble_dist):
    """
    Find bubbles using a depth first search like approach. If all incoming edges of a node are explored, it gets added to a stack.
    Nodes get popped of the stack and their neighbors get explored.
    This algorithm fails when a tip is encountered, resulting in incorrect outputs. Depth first search between the provided source and sink are needed to verify the validity of the provided bubbles.
    This algorithm is slightly modified from algorithm 6 in https://doi.org/10.1093/bioinformatics/btw152, it was modified to allow tips as bubble sinks
    :param g: graph object
    :param max_bubble_dist: maximum distance to look for bubbles
    :return:
    """
    # Find small bubbles
    # both source and sink
    # iterate over all nodes and see if multiple paths have common nodes
    # dict that stores sources and sinks of bubbles
    bubbles = {}
    nr_bubbles = 0

    for n in g.nodes:
        dist = {}
        # if outdegree < 2, it can't be a bubble
        if len(g.nodes[n].edges) > 1:

            # initalize all nodes with an infinite distance to the starting node
            for n1 in g.nodes:
                dist[n1] = float('inf')

            dist[n] = 0
            # stack contains nodes who's neighbours need to be visited
            stack = []
            nr_incoming = {}
            stack.append(n)
            # visited nodes that haven't been added to the stack yet
            p = 0

            while len(stack) != 0:

                n1 = stack.pop(-1)
                for n2 in g.nodes[n1].edges:

                    # circle found
                    if n2 == n:
                        break

                    # check if the path to this node exceeds the maximum search distance
                    if dist[n1] + g.nodes[n1].edges[n2] > max_bubble_dist:
                        break

                    # not visited before
                    if dist[n2] == float('inf'):

                        # get the indegree of n2
                        n2_rc = rc_node(n2)
                        nr_incoming[n2] = len(g.nodes[n2_rc].edges)
                        # p is the number of visited edges not added to the stack
                        # increase p, because it was not visited before
                        p += 1

                    # update distance if the current path is shorter than any other potential path or the inital infinite value
                    if dist[n1] + g.nodes[n1].edges[n2] < dist[n2]:
                        dist[n2] = dist[n1] + g.nodes[n1].edges[n2]

                    # n2 has been visited (maybe not for the first time), so the amount of unvisited incoming edges can be reduced by 1
                    nr_incoming[n2] = nr_incoming[n2] - 1

                    # if all incoming edges have been visited, n2 can be added to the stack
                    if nr_incoming[n2] == 0:
                        # in the original algorithm, we verify that n2 is not a tip, no longer the case
                        #if len(g.nodes[n2].edges) != 0:
                        stack.append(n2)
                        # n2 has been added to the stack, so p can be decreased
                        p = p - 1

                # found a sink
                if (len(stack) == 1) and (p == 0):
                    bubbles[n] = stack.pop()
                    nr_bubbles += 1

    return bubbles, nr_bubbles

def depth_first_search(g, source_node, target_node, path, all_paths):
    """
    Depth first search algorithm (recursive) to find all possible paths between a source node and a target node.
    :param g: graph object
    :param source_node: starting node of the path
    :param target_node: ending node of the path
    :param path: list of nodes that go from source to target
    :param all_paths: list of all possible paths
    :return:
    """

    # Append the current node to the path
    path.append(source_node)

    # If we've reached the target node, add the current path to the list of all paths
    if source_node == target_node:
        all_paths.append(list(path))

    else:
        # Recur for all the neighbors of the current node
        for t in g.nodes[source_node].edges:

            # Avoid circles
            if t not in path:
                depth_first_search(g, t, target_node, path, all_paths)

    # Backtrack: remove the current node from the path
    path.pop()

def remove_nodes_and_associated_edges(g, n):
    """
    Removes node n and all its associated edges from g.
    :param g: graph object
    :param n: node to be removed
    :return:
    """

    # initialize
    n_rc = rc_node(n)
    nr_edges_removed = 0
    nr_nodes_removed = 0

    # remove all edges that point to these nodes
    for n1 in g.nodes:
        if n in g.nodes[n1].edges:
            g.nodes[n1].remove_edge(n)
            nr_edges_removed += 1
        if n_rc in g.nodes[n1].edges:
            g.nodes[n1].remove_edge(n_rc)
            nr_edges_removed += 1

    # remove nodes
    g.remove_node(n)
    g.remove_node(n_rc)
    nr_nodes_removed += 2

    return g, nr_nodes_removed, nr_edges_removed

def pop_bubbles(g, bubbles):
    """
    Evaluate the bubbles given using depth first search. Pop all but the best paths through the bubble.
    :param g: graph object
    :param bubbles: dictionary containing sources and sinks of bubbles
    :return:
    """

    # initialize
    nr_nodes_removed_total = 0
    nr_edges_removed_total = 0
    nr_bubbles_popped = 0

    for bubble in bubbles:

        # check if source and sink still exist
        if (bubble in g.nodes) and (bubbles[bubble] in g.nodes):

            # find all paths between source and sink
            all_paths = []
            depth_first_search(g, bubble, bubbles[bubble], [], all_paths)

            # check if there is more than one path
            # the bubble detection algorithm sometimes gives incorrect bubbles, so this has to be verified
            if len(all_paths) > 1:
                # calculate node/edge length ratio, highest one is best
                nd_ratio = []
                for path in all_paths:

                    nr_nodes = len(path)
                    dist = 0

                    traverse_len = len(path) - 1
                    for i in range(traverse_len):
                        dist += g.nodes[path[i]].edges[path[i + 1]]

                    nd_ratio.append(nr_nodes/dist)

                # get the best path
                best_path = all_paths[nd_ratio.index(max(nd_ratio))]

                # get rid of all other paths
                for path in all_paths:
                    if path != best_path:

                        for n in path:
                            g, nr_nodes_removed, nr_edges_removed = remove_nodes_and_associated_edges(g, n)
                            nr_nodes_removed_total += nr_nodes_removed
                            nr_edges_removed_total += nr_edges_removed

                        nr_bubbles_popped += 1

    return g, nr_bubbles_popped, nr_nodes_removed_total, nr_edges_removed_total

def trim_tips(g, trim_depth):
    """
    Remove tips from the graphs.
    :param g: graph object
    :param trim_depth: size of tips that can be trimmed
    :return:
    """

    # initialize
    nr_nodes_trimmed_total = 0
    nr_edges_trimmed_total = 0

    # trim edges sequentially, start with the smallest tips, then move to bigger ones
    for current_trim_depth in range(trim_depth):

        # initialize vertices
        to_trim = {}
        for n in g.nodes:
            to_trim[n] = False

        # iterate over all vertices, if one is a tip, go back trim_depth steps. If the outdegree for any of the encountered nodes > 1 -> mark as to_trim
        for n in g.nodes:
            if not to_trim[n]:

                # check if the node is a tip potentially eligible for removal (outdegree = 0, indegree = 1)
                n_rc = rc_node(n)
                if ((len(g.nodes[n].edges.keys()) == 0) and (len(g.nodes[n_rc].edges.keys()) == 1)):
                    visited = [n]
                    reversing = 0
                    search_depth = 0

                    # check the previous nodes
                    while reversing == 0:
                        search_depth += 1

                        # check if trim depth is reached
                        if search_depth >= current_trim_depth - 1:
                            reversing = 1

                        # search one node further, the list of nodes has length 1, this was checked earlier
                        t_rc = list(g.nodes[n_rc].edges.keys())[0]
                        t = rc_node(t_rc)

                        # check if indegree and outdegree are both 1, if so, add to the visited nodes
                        if ((len(g.nodes[t].edges.keys()) == 1) and (len(g.nodes[t_rc].edges.keys()) == 1)):
                            visited.append(t)

                        # if the outdegree is > 1, nodes get trimmed
                        elif len(g.nodes[t].edges.keys()) > 1:
                            reversing = 2

                        # if the indegree is > 1, nodes don't get trimmed
                        elif (len(g.nodes[t_rc].edges.keys()) > 1):
                            reversing = 3

                    # mark nodes to trim
                    if reversing == 2:
                        for tt in visited:
                            to_trim[tt] = True

        # remove tip nodes and associated edges
        for n in to_trim.keys():
            if to_trim[n]:
                g, nr_nodes_trimmed, nr_edges_trimmed = remove_nodes_and_associated_edges(g, n)
                nr_nodes_trimmed_total += nr_nodes_trimmed
                nr_edges_trimmed_total += nr_edges_trimmed

    return g, nr_nodes_trimmed_total, nr_edges_trimmed_total

def heuristic_edge_removal(g, overlap_ratio):
    """
    Remove very long edges (indicating small overlap length) from nodes with multiple edges. This can simplify complex structures in the graph by removing weaker (shorter) overlaps.
    Overlap ratio is calculated as (len(read 1) - len(edge 1)) / (len(read 1) - len(edge 2))
    :param g: graph object
    :param overlap_ratio: overlap length ratio cutoff, edges with overlap ratio's below this one will be removed
    :return:
    """

    nr_heuristic_removed = 0
    for n in g.nodes:

        # node has 2 or more edges
        if len(g.nodes[n].edges) > 1:

            # get shortest edge, it will provide the lowest overlap ratio
            g.nodes[n].sort_edges()
            shortest_edge = list(g.nodes[n].edges.keys())[0]
            shortest_edge_length = g.nodes[n].edges[shortest_edge]

            # check overlap ratio
            read_length = len(read_dict[n[:-1]])
            i = 1
            to_remove = []
            while ((read_length - g.nodes[n].edges[list(g.nodes[n].edges.keys())[-i]]) / (read_length - shortest_edge_length)) < overlap_ratio:
                to_remove.append(list(g.nodes[n].edges.keys())[-i])
                i +=1

            # remove edges with bad overlap ratios
            for t in to_remove:
                g.nodes[n].remove_edge(t)
                nr_heuristic_removed += 1

                t_rc = rc_node(t)
                n_rc = rc_node(n)
                g.nodes[t_rc].remove_edge(n_rc)
                nr_heuristic_removed += 1

    return g, nr_heuristic_removed

def find_unitigs(g):
    """
    Creates a dict from the string graph. Iterate over all nodes, check if they're already in a unitig and create a new unitig if they're not.
    :param g: graph object
    :return:
    """

    # initialize
    unitigs = {}
    nr_unitigs = 0
    in_unitig = []

    for n in g.nodes:
        n_rc = rc_node(n)

        # check if the read is already in a unitig
        if ((n not in in_unitig) and (n_rc not in in_unitig)):

            # check if indegree and outdegree are 1
            if (len(g.nodes[n].edges) == 1) and (len(g.nodes[n_rc].edges) == 1):

                # create new unitig
                unitigs[nr_unitigs] = []
                forward_nodes = []
                reverse_nodes = []

                # keep adding reads in the forward direction until we reach the end of the unitig or until we reach a node we've already encountered
                f = list(g.nodes[n].edges.keys())[0]
                forward_adding = True
                while forward_adding:
                    targets = list(g.nodes[f].edges.keys())
                    # check if this is a tip or node with outdegree > 1
                    if len(targets) == 1:
                        f_rc = rc_node(f)
                        # check if we've already encountered this read
                        if ((f not in in_unitig) and (f_rc not in in_unitig)):
                            forward_nodes.append((f, int(g.nodes[f].edges[targets[0]])))
                            in_unitig.append(f)
                            in_unitig.append(f_rc)
                            f = targets[0]
                        else:
                            forward_adding = False
                    else:
                        forward_adding = False

                # add final node: tip or node with outdegree > 1
                # edge length can be 0, the whole read gets added in the function create_unitigs()
                if len(forward_nodes) != 0:
                    last_forward, _ = forward_nodes[-1]
                    final_node = list(g.nodes[last_forward].edges.keys())[0]
                    forward_nodes.append((final_node, 0))
                else:
                    forward_nodes.append((f, 0))

                # keep adding reads in the reverse complement direction until we reach the end of the unitig or until we reach a node we've already encountered
                r_rc = list(g.nodes[n_rc].edges.keys())[0]
                r = rc_node(r_rc)
                reverse_adding = True
                while reverse_adding:
                    targets_rc = list(g.nodes[r_rc].edges.keys())
                    # check if this is a tip or node with indegree > 1
                    if len(targets_rc) == 1:
                        # check if we've already encountered this read
                        if ((r not in in_unitig) and (r_rc not in in_unitig)):
                            reverse_nodes.append((r, int(g.nodes[r].edges[list(g.nodes[r].edges.keys())[0]])))
                            in_unitig.append(r)
                            in_unitig.append(r_rc)
                            r_rc = targets_rc[0]
                            r = rc_node(r_rc)
                        else:
                            reverse_adding = False
                    else:
                        reverse_adding = False

                # add final node: tip or node with outdegree > 1
                # this edge length has to be correct
                if len(reverse_nodes) != 0:
                    last_reverse, _ = reverse_nodes[-1]
                    final_node = rc_node(list(g.nodes[rc_node(last_reverse)].edges.keys())[0])
                    reverse_nodes.append((final_node, g.nodes[final_node].edges[last_reverse]))
                else:
                    reverse_nodes.append((r, g.nodes[r].edges[n]))

                # combine found nodes into one unitig
                edge_len_n = g.nodes[n].edges[list(g.nodes[n].edges.keys())[0]]
                in_unitig.append(n)
                in_unitig.append(n_rc)
                if len(reverse_nodes) + len(forward_nodes) + 1 > 2:
                    unitigs[nr_unitigs] = reverse_nodes[::-1] + [(n, edge_len_n)] + forward_nodes
                    nr_unitigs += 1

    return unitigs

def reverse_complement(seq):
    """
    Returns the reverse complement of a sequence
    :param seq: string, sequence to get reverse complement of
    :return:
    """

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    rev_seq = [complement[base] for base in seq]
    rev_seq.reverse()
    return ''.join(rev_seq)

def create_unitigs(unitigs, fasta_file):
    """
    Turn the unitig dictionary into a fasta file
    :param unitigs: dictionary that contains the unitig information
    :param fasta_file: file path to fasta file
    :return:
    """

    # iterate over unitigs
    for u in unitigs.keys():

        # read unitig information
        unitig_seq = ""
        for node in unitigs[u]:
            (n, edge_len) = node
            read = n[:-1]
            orientation = n[-1]
            read_seq = read_dict[read]

            # the last read can be fully added to the sequence
            if node == unitigs[u][-1]:
                edge_len = len(read_seq)

            if orientation != "+":
                read_seq = reverse_complement(read_seq)

            read_seq = read_seq[:edge_len]
            unitig_seq += read_seq

        # write to file
        if len(unitig_seq) > 1000:
            with open(fasta_file, 'a') as f:
                f.write('>' + str(u) + '\n')
                f.write(unitig_seq + '\n')

#################
# CONFIGURATION #
#################

paf_file = ''
reads_fasta_file = ''
fasta_file = ''

max_overhang = 1000
overhang_ratio = 0.8
# makes no difference between 10 and 100
fuzz = 10
trim_depth = 4
# has to be <= 1
overlap_ratio = 0.7
max_bubble_dist = 50000

###################
# TESTING GROUNDS #
###################

# load reads in memory
print("## Loading reads into memory ##")
read_dict = load_reads(reads_fasta_file)
print("## Reads loaded in memory ##\n")

# create graph
print("## Sort alignments and create string graph ##")
g, nr_internal_match, nr_first_contained, nr_second_contained, nr_proper_overlaps = create_string_graph(paf_file, max_overhang, overhang_ratio)
check_synchronization(g)
print("internal alignments: " + str(nr_internal_match) + ", contained reads: " + str(nr_first_contained + nr_second_contained) + ", proper overlaps: " + str(nr_proper_overlaps))
print('## String graph created ##\n')

# detect bubbles
print("## Detect bubbles ##")
bubbles, nr_bubbles = detect_bubbles(g, max_bubble_dist)
g, nr_bubbles_popped, nr_nodes_removed, nr_edges_removed = pop_bubbles(g, bubbles)
print("Popped " + str(nr_bubbles_popped) + " bubbles, removed " + str(nr_nodes_removed) + " nodes and " + str(nr_edges_removed) + " edges")
check_synchronization(g)
print("## Bubbles detected ##\n")

# remove transitive edges
print("## Remove transitive edges ##")
reducing = True
while reducing:
    g, nr_reduced = reduce_transitive_edges(g, fuzz)
    print("Reduced " + str(nr_reduced) + " edges")
    if nr_reduced == 0:
        reducing = False
check_synchronization(g)
print("## Transitive edges removed ##\n")

# trim tips
print("## Trimming tips")
trimming = True
while trimming:
    g, nr_nodes_trimmed, nr_edges_trimmed = trim_tips(g, trim_depth)
    print("Trimmed " + str(nr_nodes_trimmed) + " nodes and " + str(nr_edges_trimmed) + " edges")
    if nr_nodes_trimmed == 0 and nr_edges_trimmed == 0:
        trimming = False
check_synchronization(g)
print("## Tips trimmed ##\n")

# Heuristic edge removal
print('## Removing weak edges ##')
heuristic_removal = True
while heuristic_removal:
    g, nr_heuristic_removed = heuristic_edge_removal(g, overlap_ratio)
    print("Removed " + str(nr_heuristic_removed) + " edges")
    if nr_heuristic_removed == 0 :
        heuristic_removal = False
check_synchronization(g)
print("## Weak edges removed ##\n")

# detect bubbles
print("## Detect bubbles ##")
bubbles, nr_bubbles = detect_bubbles(g, max_bubble_dist)
g, nr_bubbles_popped, nr_nodes_removed, nr_edges_removed = pop_bubbles(g, bubbles)
print("Popped " + str(nr_bubbles_popped) + " bubbles, removed " + str(nr_nodes_removed) + " nodes and " + str(nr_edges_removed) + " edges")
check_synchronization(g)
print("## Bubbles detected ##\n")

# tip trimming and heuristic edge removal to simplify the graph
print("## Simplifying graph with tip trimming, transitive edge reduction and weak edge removal ##")
total_removed = float("inf")
nr_edges_trimmed = float("inf")
nr_reduced = float("inf")
nr_heuristic_removed = float("inf")

while total_removed != 0:
    if nr_edges_trimmed != 0:
        print("Tip trimming")
        g, nr_nodes_trimmed, nr_edges_trimmed = trim_tips(g, trim_depth)
        print("Trimmed " + str(nr_nodes_trimmed) + " nodes and " + str(nr_edges_trimmed) + " edges")

    if nr_reduced != 0:
        print("Transitive edge reduction")
        g, nr_reduced = reduce_transitive_edges(g, fuzz)
        print("Reduced " + str(nr_reduced) + " edges")

    if nr_heuristic_removed != 0:
        print("Weak edge removal")
        g, nr_heuristic_removed = heuristic_edge_removal(g, overlap_ratio)
        print("Removed " + str(nr_heuristic_removed) + " edges")

    total_removed = nr_edges_trimmed + nr_reduced + nr_heuristic_removed
print("## Graph simplification finished ##\n")

print("## Detect bubbles ##")
nr_edges_removed = float("inf")
while nr_edges_removed != 0:
    bubbles, nr_bubbles = detect_bubbles(g, max_bubble_dist)
    g, nr_bubbles_popped, nr_nodes_removed, nr_edges_removed = pop_bubbles(g, bubbles)
    print("Popped " + str(nr_bubbles_popped) + " bubbles, removed " + str(nr_nodes_removed) + " nodes and " + str(nr_edges_removed) + " edges")
check_synchronization(g)
print("## Bubbles detected ##\n")

# find unitigs
print("## Collect unitigs ##")
unitigs = find_unitigs(g)

# look for duplicate contigs
final_unitigs = []
for u1 in range(len(unitigs)):
    if unitigs[u1] not in final_unitigs:
        final_unitigs.append(unitigs[u1])
    else:
        print("found duplicate")

print(str(len(unitigs)) + " unitigs")

create_unitigs(unitigs, fasta_file)
print("## Unitigs collected ##\n")

#########
# STATS #
#########

print("## Gathering graph information ##")

print("nr of nodes: " + str(len(g.nodes)))

nr_edges = 0
for n in g.nodes:
    nr_edges += len(g.nodes[n].edges)

print("nr of edges: " + str(nr_edges))
print("edge/node ratio: " + str(nr_edges / len(g.nodes)))

max_len = 0

for u in unitigs.keys():
    if len(unitigs[u]) > max_len:
        max_len = len(unitigs[u])

print("longest unitig: " + str(max_len) + " reads")
print("## Graph information gathered ##")
