import build_edge_matrix as bem
import sys
import os
import pandas as pd
#import faulthandler
#faulthandler.enable()
import dendropy as dp


class build_edge_matrix_continuous(bem.build_edge_matrix):

    def __init__(self,tree_f, format, phy, recon_f, event_f, pf_ids, pt_ids = [], use_gain_events=True):
        #t = dp.Tree()
        t = dp.Tree.get(path = tree_f, schema = 'newick', suppress_internal_node_taxa=False)
        self.pf_ids = pf_ids
        self.pt_ids = pt_ids
        self.tree = t
        self.label2node = self.get_label_map()
        self.phy = phy
        self.use_gain_events = use_gain_events
        self.char2ev = self.read_events(recon_f, event_f)
        self.missing_characters = ("?", "-", "-1.0")
        self.event_f = event_f

    def read_events(self, recon_f, event_f):
        """process max likelhood phenotype events and get branch continuous feature value difference """
        char2ev = {}
        node_m = pd.read_csv(recon_f, index_col = 0, sep = "\t")
        parents = [] 
        nodes = []
        for n in node_m.index:
            if not n == "N1" and not self.label2node[n].parent_node.taxon.label == "N1":
                parents.append(self.label2node[n].parent_node.taxon.label)
                nodes.append(n)
        for char in node_m.columns:
            differences = node_m.loc[nodes, char] - node_m.loc[parents, char].values
            char2ev[char] = [(n, pn, difference) for n, pn, difference in zip(nodes, parents, differences)]
        #process max likelihood events"""
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
            #map phenotypes
            #map characters to their events, -1 to begin counting at zero as gainLoss output starts counting at 1!
            if not self.pt_ids[int(char) - 1] in char2ev:
                char2ev[self.pt_ids[int(char) - 1]] = [(node,ancestor,etype, prob, exp)]
            else: 
                char2ev[self.pt_ids[int(char) - 1]].append((node,ancestor,etype, prob, exp))
        return char2ev

    def get_label_map(self):
        """create a mapping from taxa to each node to speed up the lookup in read events"""
        l2n = {}
        for n in self.tree:
            l2n[n.taxon.label] = n
        return l2n
    
    def map_events(self, node2edge, feats, char2ev, pt, edges):
        """map feature changes on to the pt edges"""
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
                else: 
                    continue
                if ev[1] in node2edge:
                    edge2 = node2edge[ev[1]]
                else: continue
                isn = edge1.intersection(edge2)
                if len(isn) == 0:
                    print "event",ev
                    print "edge1",edge1
                    print "edge2",edge2
                    continue
                isn = isn.pop()
                if isn in edge2char2val and gt in edge2char2val[isn]:
                    edge2char2val[isn][gt] += ev[2]
                elif isn in edge2char2val and gt not in edge2char2val[isn]:
                    edge2char2val[isn][gt] = ev[2]
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
                edge2char2val[isn][pt] = float(ev[3])
            elif isn in edge2char2val and pt  in edge2char2val[isn]:
                edge2char2val[isn][pt] += float(ev[3]) * (1-edge2char2val[isn][pt])


        return edge2char2val


if __name__ == '__main__':
    import argparse 
    parser = argparse.ArgumentParser("reconstruct likelihood matrix from gainLoss output")
    parser.add_argument("tree", help='phylogenetic tree')
    parser.add_argument("phypat_pt_f", help='phylogenetic tree')
    parser.add_argument("--format", default = "newick", help='phylogenetic tree format')
    parser.add_argument("recon_f", help='table of reconstructed node feature values')
    parser.add_argument("--consider_loss_events", action = 'store_true',  help='only parse and reconstruct losses')
    parser.add_argument("event_f", help='likelihood GLOOME output file')
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
    pf_ids = pd.read_csv(a.feature_f, sep = "\t", index_col = 0).index.tolist()
    pt_ids = pd.read_csv(a.phenotype_f, sep = "\t", index_col = 0).index.astype('string').tolist()
    phy = pd.read_csv(a.phypat_pt_f, index_col = 0, sep = "\t")
    phy.index = phy.index.astype('string')
    bem = build_edge_matrix_continuous(a.tree, a.format, phy, a. recon_f, a.event_f, pf_ids, pt_ids, use_gain_events = not a.consider_loss_events)
    bem.get_all_edge_m(pf_ids, pt_ids, a.outdir)
    print "after get all edges"

