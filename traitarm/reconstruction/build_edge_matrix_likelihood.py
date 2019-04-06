import build_edge_matrix as bem
import sys
import os
import pandas as pd
#import faulthandler
#faulthandler.enable()
import dendropy as dp


class build_edge_matrix_likelihood(bem.build_edge_matrix):

    def __init__(self,tree_f, format, phy, event_f, pf_ids, pt_ids = [], use_gain_events=True, use_likelihood=True):
        #t = dp.Tree()
        #preserve underscores: make sure dendropy doesn't convert things like isol_56b to isol 56b to avoid downstream problems
        t = dp.Tree.get(path = tree_f, schema = 'newick', suppress_internal_node_taxa=False, preserve_underscores = True)
        self.pf_ids = pf_ids
        self.pt_ids = pt_ids
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
            if not (self.pf_ids + self.pt_ids)[int(char) - 1] in char2ev:
                char2ev[(self.pf_ids + self.pt_ids)[int(char) - 1]] = [(node,ancestor,etype, prob, exp)]
            else: 
                char2ev[(self.pf_ids + self.pt_ids)[int(char) - 1]].append((node,ancestor,etype, prob, exp))
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
        if pt in char2ev:
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
