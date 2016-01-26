import build_edge_matrix as bem
import sys
import named_narray as na
import os
#import faulthandler
#faulthandler.enable()
import dendropy as dp


class build_edge_matrix_likelihood(bem.build_edge_matrix):

    def __init__(self,tree_f, format,  phyn_f, phym_f, event_f, use_gain_events=True, use_likelihood=True):
        t = dp.Tree()
        self.tree = t
        t.read_from_path(tree_f, format, suppress_internal_node_taxa=False)
        self.phy = na.named_matrix(phym_f, phyn_f)
        #print self.phy.array
        self.label2node = self.get_label_map()
        self.use_gain_events = use_gain_events
        self.char2ev = self.read_events(event_f)
        self.missing_characters = ("?", "-", "-1.0")
        self.event_f=event_f
        self.use_likelihood = use_likelihood


    def get_label_map(self):
        """create a mapping from taxa to each node to speed up the lookup in read events"""
        l2n = {}
        for n in self.tree:
            l2n[n.taxon.label] = n
        return l2n

    def read_events(self, event_f):
        """process max likelihood events"""
        char2ev = {}
        f = open(event_f, 'r')
        #skip header
        f.readline()
        f.readline()
        for l in f:
            etype, char, node, _, _, _, _, prob, exp = l.strip().split('\t')
            #get predecessor of node
            #TODO the dendropy method has to iterate over all nodes in the tree, speed this up by using a mapping from the tree to each node via its label
            if etype == "gain" and not self.use_gain_events:
                continue
            elif etype == "loss" and self.use_gain_events:
                continue
            ancestor = self.label2node[node].parent_node.taxon.label
            #map characters to their events, -1 to begin counting by zero as gainLoss output starts counting at 1!
            if not int(char) - 1 in char2ev:
                char2ev[int(char) - 1] = [(node,ancestor,etype, prob, exp)]
            else: 
                char2ev[int(char) - 1].append((node,ancestor,etype, prob, exp))
        return char2ev


    def map_events(self, node2edge, gt_start, gt_end, char2ev, pt, edges):
        """map likelihood events on to the pt edges"""
        edge2char2val = {}
        for e in edges:
            edge2char2val[tuple(e)] = {}
        #map/match genotype to phenotype edges
        for gt in range(gt_start, gt_end):
            #TODO account for what happen if there are no events
            if not gt in char2ev:
                continue
            for ev in char2ev[gt]:
                #map gt event to pt event and discard if not possible
                if ev[0] in node2edge:
                    edge1 = node2edge[ev[0]]
                else: continue
                if ev[1] in node2edge:
                    edge2 = node2edge[ev[1]]
                else: continue
                isn = edge1.intersection(edge2)
                #if len(isn) == 0:
                    #print "event",ev
                    #print "edge1",edge1
                    #print "edge2",edge2
                isn = isn.pop()
                #TODO adjust summing up when using the likelihoods
                if isn in edge2char2val and gt in edge2char2val[isn]:
                    if self.use_likelihood:
                        edge2char2val[isn][gt] += float(ev[3]) * (1-edge2char2val[isn][gt])
                    else: 
                        edge2char2val[isn][gt] += float(ev[4]) 
                elif isn in edge2char2val and gt not in edge2char2val[isn]:
                    if self.use_likelihood:
                        edge2char2val[isn][gt] = float(ev[3])
                    else:
                        edge2char2val[isn][gt] = float(ev[4]) 
        #account for phenotype edges but only if this is an actual phenotype (no missing values involved)
        if pt > gt_end - 1:
            for ev in char2ev[pt]:
                if ev[0] in node2edge:
                    edge1 = node2edge[ev[0]]
                else: continue
                if ev[1] in node2edge:
                    edge2 = node2edge[ev[1]]
                else: continue
                isn = edge1.intersection(edge2).pop()
                if isn in edge2char2val and pt not in edge2char2val[isn]:
                    if self.use_likelihood:
                        edge2char2val[isn][pt] = float(ev[3])
                    else:
                        edge2char2val[isn][pt] = float(ev[4])

                elif isn in edge2char2val and pt  in edge2char2val[isn]:
                    if self.use_likelihood:
                        edge2char2val[isn][pt] += float(ev[3]) * (1-edge2char2val[isn][pt])
                    else:
                        edge2char2val[isn][pt] += float(ev[4])


        return edge2char2val

if __name__ == '__main__':
    import getopt
    if len(sys.argv) == 1:
        print """USAGE: python %s
-t <tree> in a dendropy compliant format
-f <format> above format either one of nexus, newick, nexml, phylip
-p <phyletic pattern> phyletic_patterns in fasta format
-e <events> input file with max likelihood events
-l <loss_events> only consider the loss events (default only gain events)
-x <use exepected number of events> (default probability)
-g <range of genotypes> e.g. 1-8400
-h <range of phenotypes> to consider e.g 8550-8560
-o <out_dir> for the output matrices
        """ % (sys.argv[0])
        #testing only: uncomment and run without cmdargs
        #t = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/gainLoss_input_v2_candidatus_sample30/stol_2_bioprojects_20140115_RefSeq_genome_NCBI20140115_gideon.tre.internal_nodes_labeled.newick"
        #e = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/gainLoss_input_v2_candidatus_sample30/RESULTS/gainLossProbExpPerPosPerBranch.txt"
        #phy_m = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/pfams_pts_counts.tsv"
        #phy_n = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/stol_2_bioprojects_20140115_RefSeq_genome_NCBI20140115_gideon_sp+st.taxid_only.txt"
        #f = "newick"
        #bem = build_edge_matrix_likelihood(t,f,phy_n, phy_m,e)
        #bem.get_all_edge_m(0,10,8477,8477, "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/likelihood/gain")
        sys.exit(1)
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "t:f:p:n:e:xlg:h:o:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    t=None
    f=None
    p=None
    e = None
    g1 = g2 = None
    pt1 = pt2 = None
    n = None
    out = None
    use_likelihood = True 
    use_gain_events = True 
    for o, a in optlist:
        if o == "-t":
            t = a
        if o == "-f":
            f = a
        if o == "-p":
            p = a
        if o == "-n":
            n = a
        if o == "-e":
            e = a
        if o == "-x":
            use_likelihood = False 
        if o == "-l":
            use_gain_events = False
        if o == "-g":
            g1, g2 = [int(i) for i in a.split("-")]
        if o == "-h":
            pt1, pt2 = [int(i) for i in a.split("-")]
        if o == "-o":
            out = a
            #check if the directory already exists
            if os.path.exists(out):
                sys.stderr.write("output directory %s already exists; delete and rerun\n"%a)
                sys.exit(1)
            else:
                os.mkdir(out)
    bem = build_edge_matrix_likelihood(t,f,n,p,e, use_likelihood=use_likelihood, use_gain_events=use_gain_events)
    bem.get_all_edge_m(g1,g2,pt1,pt2, out)
    print "after get all edges"

