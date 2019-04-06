import pandas as pd
import ete2
from ete2 import faces, Tree, AttrFace, TreeStyle
import pylab
from matplotlib.colors import hex2color, rgb2hex, hsv_to_rgb, rgb_to_hsv

kelly_colors_hex = [
    0xFFB300, # Vivid Yellow
    0x803E75, # Strong Purple
    0xFF6800, # Vivid Orange
    0xA6BDD7, # Very Light Blue
    0xC10020, # Vivid Red
    0xCEA262, # Grayish Yellow
    0x817066, # Medium Gray

    # The following don't work well for people with defective color vision
    0x007D34, # Vivid Green
    0xF6768E, # Strong Purplish Pink
    0x00538A, # Strong Blue
    0xFF7A5C, # Strong Yellowish Pink
    0x53377A, # Strong Violet
    0xFF8E00, # Vivid Orange Yellow
    0xB32851, # Strong Purplish Red
    0xF4C800, # Vivid Greenish Yellow
    0x7F180D, # Strong Reddish Brown
    0x93AA00, # Vivid Yellowish Green
    0x593315, # Deep Yellowish Brown
    0xF13A13, # Vivid Reddish Orange
    0x232C16, # Dark Olive Green
    ]

def my_layout(node):
    if node.is_leaf():
         # If terminal node, draws its name
         name_face = AttrFace("name")
    else:
         # If internal node, draws label with smaller font size
         name_face = AttrFace("name", fsize=10)
    # Adds the name face to the image at the preferred position
    faces.add_face_to_node(name_face, node, column=0, position="branch-right")
    
def adjust_kelly_brightness(hex_color, val, recon_min, recon_max):
    """set brightness according to change in continuous reconstruction value"""
    h, s, v =  rgb_to_hsv(hex2color('#{0:06X}'.format(hex_color)))
    scale_factor = 1 - (recon_max - val) / (recon_max - recon_min)
    v_new = v - (v * (scale_factor)) 
    return rgb2hex(hsv_to_rgb(pd.np.array([h, s, v_new])))

def get_style():
    ts = TreeStyle()
    # Do not add leaf names automatically
    ts.show_leaf_name = False
    ts.show_scale = True 
    ts.force_topology = False 
    # Use my custom layout
    ts.layout_fn = my_layout
    return ts

def plot_tree(pt_tree, target_node, out):
    #pt_tree, feats, pf2color = get_tree(phenotype = phenotype, feat_list = "top_cor", is_ml_plus_phypat = True, target_node = target_node)
    pt_tree.dist = 0
    target = pt_tree.search_nodes(name = target_node)[0]
    target.render(out + '_tree.pdf',  tree_style = get_style())
    #target.render(out + '_tree.png', tree_style = get_style())
    return target, feats, pf2color

def plot_legend(feats, out, pf2color,  pf_desc = False, pf_acc = True, include_class = False):
    fig = pylab.figure()
    figlegend = pylab.figure(figsize = (9, 6))
    ax = fig.add_subplot(111)
    x = [0,1]
    lines = [ax.plot(x, pd.np.ones(len(x)), 'o', color = "#%06x" % (pf2color[feats.index[i]]))[0] for i in range(len(pf2color))]
    labels= [i for i in feats.index]
    #labels= ["%s" %(feats.loc[:,"Pfam_acc"].iloc[i]) for i in range(feats.shape[0])]
    #if include_class:
    #    labels= ["%s %s" %(labels[i], feats.loc[:, "class"].iloc[i]) for i in range(len(labels))]
    #if pf_desc:
    #    labels = ["%s %s" % (labels[i], pf2short_desc.loc[feats.loc[:,"Pfam_acc"].iloc[i], 1]) for i in range(len(labels))]
    #if pf_acc:
    #    labels = ["%s %s" % (labels[i], pf2acc.loc[feats.loc[:,"Pfam_acc"].iloc[i], 1]) for i in range(len(labels))]

    figlegend.legend(lines, labels, markerscale = 2.5, numpoints = 1, frameon = False)
    #fig.show()
    fig.tight_layout()
    figlegend.savefig(out + "_legend.svg")
    figlegend.savefig(out + "_legend.png")
    return figlegend

def get_tree(phenotype, tree, gain_recon, loss_recon, node_recon, pfam_mapping, feat_list, sample_mapping, threshold = 0.5, target_node = None, are_continuous_features_with_discrete_phenotype = False, max_feats = 10, miscl = None, node_annotation = None):
    #read target feats
    feats = pd.read_csv(feat_list, index_col = 0, sep = "\t")
    pt_tree = ete2.Tree(tree, format = 1)
    pt_tree.ladderize()
    if not node_annotation is None:
        node_table = pd.read_csv(node_annotation, sep = "\t", index_col = 0)
    sample_mapping = pd.read_csv(sample_mapping, index_col = 0, sep = "\t")
    #read node and edge reconstruction matrices
    node_recon = pd.read_csv(node_recon, sep = "\t", index_col = 0)
    gain_recon = pd.read_csv(gain_recon, sep = "\t", index_col = 0)
    gain_recon.index = ["_".join(("_".join(i.split("_")[:-1]), i.split("_")[-1])) for i in gain_recon.index.values]
    loss_recon = pd.read_csv(loss_recon, sep = "\t", index_col = 0)
    loss_recon.index = ["_".join(("_".join(i.split("_")[:-1]), i.split("_")[-1])) for i in loss_recon.index.values]
    #prune to target node
    if target_node is not None:
        pt_tree = pt_tree.search_nodes(name = target_node)[0]
    node2name = dict((i.name, i.name) for i in pt_tree.traverse(strategy = 'preorder'))
    pfams_with_event = set()
    pfam2color = {}
    #set the style of the branches and nodes according to the posterior probability
    top10_feats = feats.iloc[:max_feats,]
    #for visualization of continuous feature get the range of values for each feature
    if are_continuous_features_with_discrete_phenotype:
        recon_min = gain_recon.abs().apply(pd.np.min)
        recon_max = gain_recon.abs().apply(pd.np.max)
    if not miscl is None:
        miscl_m = pd.read_csv(miscl, sep = "\t", index_col = 0)
    for n in pt_tree.traverse():
        #ignore the root
        if n.name == "N1":
            continue
        if not node_annotation is None:
            if n.name in node_table.index:
                for attr,i  in zip(node_table.columns, range(len(node_table.columns))):
                    value = node_table.loc[n.name, attr]
                    if not pd.isnull(value):
                        if value == 0:
                            rf = ete2.CircleFace(radius = 8, style = "circle", color = 'red')
                        elif value == 2:
                            rf = faces.CircleFace(radius = 8, style = "circle", color = 'orange')
                        else:
                            rf = faces.CircleFace(radius = 8, style = "circle", color = 'green')
                    else:
                        rf = faces.CircleFace(radius = 8, style = "circle", color = 'grey')
                    n.add_face(rf, column = i, position = "aligned")

        ns = node_recon.loc[n.name, phenotype] 
        style = ete2.NodeStyle()
        style["shape"] = 'square'
        style['size'] = 10
        if pd.isnull(ns):
            style['fgcolor'] = 'grey'
        elif ns < threshold:
            style['fgcolor'] = 'darkred'
        else:
            style['fgcolor'] = 'green'
        if not n.name == "N1":
            branch_id = n.name + "_" + n.up.name
            if gain_recon.loc[branch_id, phenotype] > threshold:
                style["hz_line_type"] = 1
                style["hz_line_color"] = 'green' 
                style["hz_line_width"] = 3
            elif loss_recon.loc[branch_id, phenotype] > threshold:
                style["hz_line_type"] = 1
                style["hz_line_color"] = 'red'
                style["hz_line_width"] = 3
            else:
                style["hz_line_type"] = 0
                style["hz_line_color"] = 'black'
            n.set_style(style)
            #check if sample was misclassified and add misclassified label
            if not miscl is None:
                if node2name[n.name] in miscl_m.index:
                    tf = faces.TextFace("misclassified")
                    n.add_face(tf, column = 0, position = "branch-right")
            #set species name instead of tax id
            if n.name in sample_mapping.index:
                node2name[n.name] = sample_mapping.loc[n.name,][0]
            #add majority feature gains and losses
            events = []
            for i in range(top10_feats.shape[0]): 
                if not are_continuous_features_with_discrete_phenotype:
                    cf = faces.CircleFace(radius = 8, style = "circle", color = kelly_colors_hex[i])
                    #gain events
                    if gain_recon.loc[branch_id, top10_feats.index[i]] > threshold:
                        pfam2color[top10_feats.index[i]] = kelly_colors_hex[i]
                        tf = faces.TextFace("-")
                        events.append(tf)
                        pfams_with_event.add(node_recon.index[i])
                        events.append(cf)
                    #loss events
                    elif loss_recon.loc[branch_id, top10_feats.index[i]] > threshold:
                        pfam2color[top10_feats.index[i]] = kelly_colors_hex[i]
                        tf = faces.TextFace("-")
                        events.append(tf)
                        pfams_with_event.add(node_recon.index[i])
                        events.append(cf)
                #continuous features
                else: 
                    adjusted_color = adjust_kelly_brightness(kelly_colors_hex[i], abs(loss_recon.loc[branch_id, top10_feats.index[i]]), recon_min.loc[top10_feats.index[i]], recon_max.loc[top10_feats.index[i]])
                    #tf = faces.TextFace(gain_recon.loc[branch_id, top10_feats.index[i]])
                    if loss_recon.loc[branch_id, top10_feats.index[i]] < 0:
                        tf = faces.TextFace("-")
                    else:
                        tf = faces.TextFace("+")
                    cf = faces.CircleFace(radius = 8, style = "circle", color = adjusted_color)
                    pfam2color[top10_feats.index[i]] = kelly_colors_hex[i]
                    pfams_with_event.add(node_recon.index[i])
                    events.append(cf)
                    events.append(tf)
            for i in range(len(events)):
                n.add_face(events[i], column = i, position = "branch-top")
    for n in pt_tree.traverse():
        if n.name in node2name:
            n.name = node2name[n.name]
    #filtered_pfams = filter(lambda i: i in list(pfams_with_event), top10_feats.loc[:,"Pfam_acc"].values)
    #print filtered_pfams
    #filtered_ids = pt_gt2id.loc[filtered_pfams, 0] - 1
    #print filtered_ids
    #top10_feats_with_event = top10_feats.loc[filtered_ids,]
    #process node annotation
    return pt_tree, top10_feats, pfam2color

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("""visualize target list of features""")
    parser.add_argument("node_recon", help = "node ancestral character state reconstruction") 
    parser.add_argument("gain_recon", help = "gain events ancestral character state reconstruction") 
    parser.add_argument("loss_recon", help = "loss events ancestral character state reconstruction") 
    parser.add_argument("tree", help = "tree with internal nodes labeled") 
    parser.add_argument("pfam_mapping",  help = "feature mapping/list")
    parser.add_argument("feat_list",  help = "list of features")
    parser.add_argument("--target_node", default = "N1", help = "list of features")
    parser.add_argument("phenotype",  help = "target phenotype")
    parser.add_argument("--are_continuous_features_with_discrete_phenotype", action = 'store_true',  help = "set if using continuous features with a discrete phenotype")
    parser.add_argument("threshold", type = float,  help = "threshold to call genotype/phenotype events")
    parser.add_argument("sample_mapping",  help = "mapping between sample ids and names")
    parser.add_argument("out",  help = "output file")
    parser.add_argument("--max_feats", type = int, default = 10, help = "visualize at most max_feats features")
    parser.add_argument("--miscl", help = "table of misclassified samples")
    parser.add_argument("--node_annotation", help = "table of binary features for labeling the nodes")
    a = parser.parse_args()
    pt_tree, feats, pf2color = get_tree(node_recon = a.node_recon, gain_recon = a.gain_recon, loss_recon = a.loss_recon, pfam_mapping = a.pfam_mapping, tree = a.tree, feat_list = a.feat_list, phenotype = a.phenotype, target_node = a.target_node, threshold = a.threshold, sample_mapping = a.sample_mapping, are_continuous_features_with_discrete_phenotype = a.are_continuous_features_with_discrete_phenotype, max_feats = a.max_feats, miscl = a.miscl, node_annotation = a.node_annotation)
    plot_tree(pt_tree, a.target_node, a.out)
    plot_legend(feats, a.out, pf2color) 
