import build_edge_matrix as bem
import sys
import named_narray as na
import os
import pandas as pd
#import faulthandler
#faulthandler.enable()
import dendropy as dp


class build_edge_matrix_likelihood(bem.build_edge_matrix):

    def __init__(self,tree_f, format, phy, event_f, pf_ids, pt_ids, use_gain_events=True, use_likelihood=True):
        #t = dp.Tree()
        t = dp.Tree.get(path = tree_f, schema = 'newick', suppress_internal_node_taxa=False)
        self.tree = t
        self.label2node = self.get_label_map()
        self.phy = phy
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
            if not (pf_ids.tolist() + pt_ids.tolist())[int(char) - 1] in char2ev:
                char2ev[(pf_ids.tolist() + pt_ids.tolist())[int(char) - 1]] = [(node,ancestor,etype, prob, exp)]
            else: 
                char2ev[(pf_ids.tolist() + pt_ids.tolist())[int(char) - 1]].append((node,ancestor,etype, prob, exp))
        return char2ev


    def map_events(self, node2edge, feats, char2ev, pt, edges):
        """map likelihood events on to the pt edges"""
        edge2char2val = {}
        for e in edges:
            edge2char2val[tuple(e)] = {}
        #map/match genotype to phenotype edges
        for gt in feats:
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
    import argparse 
    parser = argparse.ArgumentParser("reconstruct likelihood matrix from gainLoss output")
    parser.add_argument("tree", help='phylogenetic tree')
    parser.add_argument("phypat_pt_f", help='phylogenetic tree')
    parser.add_argument("--format", default = "newick", help='phylogenetic tree format')
    parser.add_argument("event_f", help='gainLoss tabular output file')
    parser.add_argument("--consider_loss_events", action = 'store_true',  help='only parse and reconstruct losses')
    parser.add_argument("--use_expected", action = "store_true", help = "--use expected number of events (default probabilities")
    parser.add_argument("feature_f", help = "list of features used")
    parser.add_argument("phenotype_f", help = "list of phenotypes")
    parser.add_argument("outdir", help = "output directory")
    a = parser.parse_args()
    #check if the directory already exists
    if os.path.exists(a.outdir):
        #sys.stderr.write("output directory %s already exists; delete and rerun\n"%a.outdir)
        #sys.exit(1)
        pass
    else:
        os.mkdir(a.outdir)
    #read pf and pt file
    pf_ids = pd.read_csv(a.feature_f, sep = "\t", index_col = 0).index
    pt_ids = pd.read_csv(a.phenotype_f, sep = "\t", index_col = 0).index.astype('string')
    phy = pd.read_csv(a.phypat_pt_f, index_col = 0, sep = "\t")
    phy.index = phy.index.astype('string')
    bem = build_edge_matrix_likelihood(a.tree, a.format, phy, a.event_f, pf_ids, pt_ids, use_likelihood = not a.use_expected, use_gain_events = not a.consider_loss_events)
    bem.get_all_edge_m(pf_ids, pt_ids, a.outdir)
    print "after get all edges"

