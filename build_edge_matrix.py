import sys
import dendropy as dp
class build_edge_matrix:

    def __init__(self,tree_f, format, phyletic_f, event_f):
        t = dp.Tree()
        self.tree = t
        t.read_from_path(tree_f, format, suppress_internal_node_taxa=False)
        self.phy = dp.StandardCharacterMatrix.get_from_path(phyletic_f, "fasta")
        self.char2ev = self.read_events(event_f)
        self.missing_characters = ("?", "-")
        self.event_f=event_f

    def get_edges(self, pt):
        edges =[]
        #recurse through tree to get edges
        tl, tr = self.tree.seed_node.child_nodes()
        self.get_edges_h( tl, tr,self.tree.seed_node, edges, pt)
        return edges

    def get_edges_h(self, l, r,node, edges, pt):
        print "getting edges for phenotype", pt
        if l.is_leaf() and r.is_leaf():
            return self.update(l, r, node, edges, pt)
        elif l.is_leaf():
            rl, rr = r.child_nodes()
            return self.update(l, self.get_edges_h(rl, rr, r, edges, pt), node, edges, pt)
        elif r.is_leaf():
            ll, lr = l.child_nodes()
            return self.update(self.get_edges_h(ll, lr, l, edges, pt), r, node, edges, pt)
        else:
            ll, lr = l.child_nodes()
            rl, rr = r.child_nodes()
            return self.update(self.get_edges_h(ll, lr, l, edges, pt), self.get_edges_h(rl, rr,r , edges, pt), node, edges, pt)


    def update(self, l, r,node, edges, pt):
        #print type(node.taxon), type(l), type(r)
        #both children stem from missing values (NAs)
        if node.taxon.label == 'N45':
            print l, r
        l1, r1 = node.child_nodes()
        if  l1.taxon.label == 'N45':
            print "links", l,r
        if  r1.taxon.label == 'N45':
            print "rechts", l,r
        if l is None and r is None:
            return None
        #node l is ancestor of NA node return the right node
        elif l is None:
            if type(r) == type([]):
                r.append(r1.taxon.label)
                return r
            elif type(r) == type(True):
                return [r1.taxon.label]
            elif type(r) == type(node):
                if str(self.phy[r.taxon.label][pt]) not in self.missing_characters:
                    return [r.taxon.label]
                else:
                    return None
        #node r is ancestor of NA node return the left node
        elif r is None:
            if type(l) == type([]):
                l.append(l1.taxon.label)
                return l
            elif type(l) == type(True):
                return [l1.taxon.label]
            elif type(l) == type(node):
                if str(self.phy[l.taxon.label][pt]) not in self.missing_characters:
                    return [l.taxon.label]
                else:
                    return None
        #one node is list the other leave
        elif type(l) == type(node) and type(r) == type([]):
            if str(self.phy[l.taxon.label][pt]) not in self.missing_characters:
                r.append(node.taxon.label)
                edges.append(r)
                edges.append((l.taxon.label, node.taxon.label))
                return True
            else:
                r.append(r1.taxon.label)
                return r
        elif type(r) == type(node) and type(l) == type([]):
            if str(self.phy[r.taxon.label][pt]) not in self.missing_characters:
                l.append(node.taxon.label)
                edges.append(l)
                edges.append((r.taxon.label, node.taxon.label))
                return True
            else:
                l.append(l1.taxon.label)
                return l
        #both nodes derive from non NA nodes
        #one of the nodes is sane and isn't an ancestor of a NA node
        elif type(l) == type([]) and type(r) == type([]):
            l.append(node.taxon.label)
            r.append(node.taxon.label)
            edges.append(l)
            edges.append(r)
            return True
        elif type(l) == type(True) or type(r) == type(True):
            #check if left node is either type node or list of nodes
            if type(l) != type(True):
                #check if left node is leaf node
                if type(l) == type(l1) and l.is_leaf():
                    if str(self.phy[l.taxon.label][pt]) not in self.missing_characters:
                        edges.append((l.taxon.label, node.taxon.label))
                        edges.append((r1.taxon.label, node.taxon.label))
                        return True
                    else: return [r1.taxon.label]
                #check if left node is type list
                elif type(l) == type([]):
                    l.append(node.taxon.label)
                    edges.append(l)
                    edges.append((r1.taxon.label, node.taxon.label))
                    return True
            #check if right node is either type node or list of nodes
            if type(r) != type(True):
                #check if right node is leave node
                if type(r) == type(r1) and r.is_leaf():
                    if str(self.phy[r.taxon.label][pt]) not in self.missing_characters:
                        edges.append((r.taxon.label, node.taxon.label))
                        edges.append((l1.taxon.label, node.taxon.label))
                        return True
                    else: return [l1.taxon.label]
                #check if right node is type list
                elif type(r) == type([]):
                    r.append(node.taxon.label)
                    edges.append(r)
                    edges.append((l1.taxon.label, node.taxon.label))
                    return True
            else:
                edges.append((l1.taxon.label,node.taxon.label))
                edges.append((r1.taxon.label, node.taxon.label))
                return True
        #both nodes are leave nodes
        elif l.is_leaf() and r.is_leaf():
            if str(self.phy[l.taxon.label][pt]) not in self.missing_characters and str(self.phy[r.taxon.label][pt]) not in self.missing_characters:
                edges.append((l.taxon.label,node.taxon.label))
                edges.append((r.taxon.label,node.taxon.label))
                return True
            elif str(self.phy[l.taxon.label][pt]) not in self.missing_characters:
                return [l.taxon.label, node.taxon.label]
            elif str(self.phy[r.taxon.label][pt]) not in self.missing_characters:
                return [r.taxon.label, node.taxon.label]
            else:
                return None


    def read_events(self,event_f):
        """read parsimony events for all characters"""
        char2ev = {}
        f = open(event_f, 'r')
        #skip header
        f.readline()
        for l in f:
            char, n1,n2,etype, ecount = l.strip().split("\t")
            if not int(char) in char2ev:
                char2ev[int(char)] = [(n1,n2,etype,ecount)]
            else: char2ev[int(char)].append((n1,n2,etype,ecount))
        return char2ev

    def get_pt_dict(self,edges):
        """for a given pt character get a mapping from nodes to edges """
        node2edge = {}
        for e in edges:
            for node in e:
                if not node in node2edge:
                    node2edge[node] = set()
                    node2edge[node].add(tuple(e))
                else:
                    node2edge[node].add(tuple(e))

        return node2edge

    def map_events(self,node2edge, gt_start, gt_end, char2ev, pt, edges):
        """map events from all  genes in gt_range to pt edges and add pt states too"""
        edge2char2val = {}
        for e in edges:
            edge2char2val[tuple(e)] = {}
        #map/match genotype to phenotype edges
        for gt in range(gt_start, gt_end+1):
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
                #print edge1, edge2
                isn = edge1.intersection(edge2)
                if len(isn) == 0:
                    print "event",ev
                    print "edge1",edge1
                    print "edge2",edge2
                isn = isn.pop()
                if isn in edge2char2val and gt in edge2char2val[isn]:
                    if ev[2] == 'loss':
                        edge2char2val[isn][gt] -= int(ev[3])
                    else: edge2char2val[isn][gt] += int(ev[3])
                elif isn in edge2char2val and gt not in edge2char2val[isn]:
                    if ev[2] == 'loss':
                        edge2char2val[isn][gt] = -int(ev[3])
                    else: edge2char2val[isn][gt] = int(ev[3])
                    #account for phenotype edges
        #print "edge2char2val", edge2char2val
        for ev in char2ev[pt]:
            print ev[0], ev[1]
            edge1 = node2edge[ev[0]]
            edge2 = node2edge[ev[1]]
            isn = edge1.intersection(edge2).pop()
            if ev[2] == 'loss':
                edge2char2val[isn][pt]=-int(ev[3])
            else: edge2char2val[isn][pt]=int(ev[3])
        return edge2char2val

    def get_edge_m(self,edge2char2val, edges, gt_start, gt_end,pt, out_f):
        """generate for each edge a vector  of all characters"""
        out_fo = open(out_f, 'w')
        for e in edges:
            s="%s\t"%str("_".join(e))
            for gt in range(gt_start, gt_end+1):
                if gt in edge2char2val[tuple(e)]:
                    s+="%s\t"%edge2char2val[tuple(e)][gt]
                else:
                    s+="0\t"
            if pt in edge2char2val[tuple(e)]:
                s+="%s\n"%edge2char2val[tuple(e)][pt]
            else: s+="0\n"
            out_fo.write(s)
    def get_all_edge_m(self,gt_start,gt_end, pt_start, pt_end, out_dir):
        """for all phenotypes generate a edge based matrix"""
        for pt in range(pt_start, pt_end+1):
            print "current phenotype number being processed", pt
            edges = self.get_edges(pt)
            print edges
            pt_dict = self.get_pt_dict(edges)
            edge2char2val =self.map_events(pt_dict, gt_start,gt_end, self.char2ev,pt, edges)
            #TODO check if outdir exists otherwise create
            self.get_edge_m(edge2char2val, edges, gt_start,gt_end,pt, "%s/pt%s.dat"%(out_dir,pt))

if __name__ == '__main__':
    import getopt
    if len(sys.argv) == 1:
        print """USAGE: python %s
-t <tree> in a dendropy compliant format
-f <format> above format either one of nexus, newick, nexml, phylip
-p <phyletic pattern> phyletic_patterns in fasta format
-e <events> file with parsimony events
-g <range of genotypes> e.g. 1-8400
-pt <range of phenotypes> to consider e.g 8550-8560
> <out file> with gain and loss events
        """ % (sys.argv[0])
    #t = "/net/metagenomics/projects/phenotypes_20130523/code/build_edge_m/sampletree.newick"
    #p = "/net/metagenomics/projects/phenotypes_20130523/code/build_edge_m/charm.fasta"
    #e = "/net/metagenomics/projects/phenotypes_20130523/code/build_edge_m/events.txt"
    #f = "newick"
    #bem = build_edge_matrix(t,f,p,e)
    ##test edge generation
    #e = bem.get_edges(1)
    #print e
    #for i in e:
    #    s=""
    #    for j in i:
    #        s+="\t%s"%j
    #    print s

    #test edge matrix generation
    #char2ev =  bem.char2ev
    #pt_dict = bem.get_pt_dict(e)
    #edge2char2val = bem.map_events(pt_dict, 0,1, char2ev,1, e)
    #print edge2char2val
    #bem.get_edge_m(edge2char2val, e, 0,1,1, "out_matrix.test.txt")

    #test bulk edge m creation
    #bem.get_all_edge_m(0,1,1,1,"out_matrices")

    #test gideon bulk edge m creation

    #t = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/gainLoss_input_v2_candidatus_sample30/stol_2_bioprojects_20140115_RefSeq_genome_NCBI20140115_gideon.tre.internal_nodes_labeled.newick"
    #p = "/net/metagenomics/projects/phenotypes_20130523/code/cordero_parsimony/test_examples/gideon_candidatus_sample30_pfams.fasta"
    #e = "/net/metagenomics/projects/phenotypes_20130523/code/cordero_parsimony/test_examples/gideon_candidatus_sample30_pfams.events"
    #f = "newick"
    #bem = build_edge_matrix(t,f,p,e)
    #bem.get_all_edge_m( 0,10,560,561, "out_matrices")

    t = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/gainLoss_input_v2_candidatus_sample30/stol_2_bioprojects_20140115_RefSeq_genome_NCBI20140115_gideon.tre.internal_nodes_labeled.newick"
    p = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/pfams_pts.fasta"
    e = "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/parsimony/events_g2_l1.txt"
    f = "newick"
    bem = build_edge_matrix(t,f,p,e)
    bem.get_all_edge_m(0,8475,8476,8568, "/net/metagenomics/projects/phenotypes_20130523/gideon/mapping/stol_2_NCBI20140115_candidatus_sample30/parsimony/input_g2_l1/")
    sys.exit(1)
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "t:f:p:e:")
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    t=None
    f=None
    p=None
    e = None
    for o, a in optlist:
        if o == "-t":
            t = a
        if o == "-f":
            f = a
        if o == "-p":
            p = a
        if o == "-e":
            e = a
    bem = build_edge_matrix(t,f,p,e)

