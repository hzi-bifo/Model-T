import pandas as ps
import numpy as np
import ete2
from ete2 import TreeStyle, faces, AttrFace, CircleFace, TextFace
import os
import colorsys

def map_bacc2color(val, minval, maxval):
    #color range from colorbrewer2.org
    #range from 247, 252, 185 in RGB to 49,163,84 in hsv 
    #convert val in range minval..maxval to the range hue 0..360 degrees, 0..100 saturation and 0...100 value 
    #if minval is not in the prespecified range set to minvalue; same for max value
    if val < minval:
        val = minval
    if val > maxval:
        val = maxval
    h = (138 - (1 - (float(val-minval) / (maxval-minval))) * 74) / 360
    s = (69.9 - (1 - (float(val-minval) / (maxval-minval))) * 42.3) / 100 
    v = (98.8 - (float(val-minval) / (maxval-minval)) * 34.9) / 100
    # convert hsv color (h,s,v) to its rgb equivalent
    return  '#%02x%02x%02x' %tuple(map(lambda x: int(x*255), colorsys.hsv_to_rgb(h, s, v)))


def miscl_layout(node):
    if node.is_leaf():
        name_face = AttrFace("miscl_full")
        faces.add_face_to_node(name_face, node, column=0, position="branch-top")
    else:
        name_face = AttrFace("miscl_full", fsize = 10)
        faces.add_face_to_node(name_face, node, column=0, position="branch-top")
    if not ps.isnull(miscl_phypat_ml.loc[node.name, "macc"]): 
        #make size of radius depend on number of samples per clade
        radius = 8
        if 'gc' in miscl_phypat_ml.columns and miscl_phypat_ml.loc[node.name, 'gc'] !=  0 :
            import math
            radius = 8 + (math.log(miscl_phypat_ml.loc[node.name, 'gc'])) * (8.0 / math.log(miscl_phypat_ml.loc[:, 'gc'].max()))
            print miscl_phypat_ml.loc[node.name, 'gc']
            print radius
        #make size of radius depend on number of labels per clade
        #radius = 8 + (8 * min(miscl_phypat_ml.loc[node.name, 'positve_samples'] + miscl_phypat_ml.loc[node.name, 'negative_samples'], 4000) / 4000)
        c_mp = CircleFace( radius = radius , color = map_bacc2color(miscl_phypat_ml.loc[node.name, "macc"], 0.75, 1.0))
        c_p = CircleFace( radius = radius, color = map_bacc2color(miscl_phypat.loc[node.name, "macc"], 0.75, 1.0))
    else:
        c_mp = CircleFace( radius = 8, color = 'grey')
        c_p = CircleFace( radius = 8,  color = 'grey')
    c_p.opacity = 0.9
    c_mp.opacity = 0.9
 
    faces.add_face_to_node(c_mp, node, column=1, position="branch-bottom")
    faces.add_face_to_node(c_p, node, column=2, position="branch-bottom")

if __name__ == "__main__":
    import argparse
    FULL = False
    parser = argparse.ArgumentParser("visualize ancestral classification performance")
    parser.add_argument("phypat", help='the input directory with the extant misclassified statistics for phypat ')
    parser.add_argument("phypat_ml", help='the input directory with the extant misclassified statistics for phypat_ml')
    parser.add_argument("tree_f", help = 'tree file')
    parser.add_argument("-f", "--more_info", help = "add statistics about false positives / false negatives to the nodes", action = 'store_true')
    parser.add_argument("-r", "--rank", help = "specify rank at which to cut the tree", default = "genus")
    parser.add_argument("out_f", help='output image path')
    args = parser.parse_args()
    miscl_phypat = ps.read_csv(args.phypat, index_col = 0, sep = "\t")
    miscl_phypat_ml = ps.read_csv(args.phypat_ml, index_col = 0, sep = "\t")
    t = ete2.Tree(args.tree_f, 1)
    t.scientific_name = "Bacteria"
    t.name = "NoName"
    t.rank = "Kingdom"
    #bacc_max = max(miscl_phypat['0'].max(), miscl_phypat_ml['0'].max())
    #bacc_min = min(miscl_phypat['0'].min(), miscl_phypat_ml['0'].min())
    #prune tree to genus level
    #get genus nodes
    genera = []
    for i in t.traverse():
        if i.rank == args.rank:
            genera.append(i)
    t.prune(genera)
    #join scientific name and miscl statistics
    for n in t.traverse():
        if args.more_info:
            n.miscl_full = "%s | %i %i %i %i %s | %i %i %i %i %s" % tuple([n.scientific_name] + [round(miscl_phypat_ml.loc[n.name,].iloc[i], 2) for i in range(5)] + [round(miscl_phypat.loc[n.name,].iloc[i], 2) for i in range(5)])
        else:
            #n.miscl_full = "%s| %s%% | %s%% | %s" % (n.scientific_name , round(miscl_phypat_ml.loc[n.name, "macc"], 3) * 100,  round(miscl_phypat.loc[n.name, "macc"], 3) * 100, miscl_phypat_ml.loc[n.name, "gc"])
            n.miscl_full = "%s" % n.scientific_name

    ts = TreeStyle()
    #ts.rotation = 90
    ts.layout_fn = miscl_layout
    #don't use the default node names
    ts.show_leaf_name = False
    #don't show scale (unit branch lengths)
    ts.show_scale = False
    ts.title.add_face(TextFace("Accuracy\t"), column = 0)
    ts.title.add_face(CircleFace( radius = 8, color = map_bacc2color(0.75, 0.75, 1.0)), column=1)
    ts.title.add_face(TextFace("75%"), column = 1)
    ts.title.add_face(CircleFace( radius = 8, color = map_bacc2color(0.8, 0.75, 1.0)), column=2)
    ts.title.add_face(TextFace("80%"), column = 2)
    ts.title.add_face(CircleFace( radius = 8, color = map_bacc2color(0.85, 0.75,  1)), column=3)
    ts.title.add_face(TextFace("85%"), column = 3)
    ts.title.add_face(CircleFace( radius = 8, color = map_bacc2color(0.90, 0.75,  1.0)), column=4)
    ts.title.add_face(TextFace("90%"), column = 4)
    ts.title.add_face(CircleFace( radius = 8, color = map_bacc2color(0.95, 0.75,  1)), column=5)
    ts.title.add_face(TextFace("95%"), column = 5)
    ts.title.add_face(CircleFace( radius = 8, color = map_bacc2color(1.0, 0.75,  0.99)), column=6)
    ts.title.add_face(TextFace("100%"), column = 6)

    t.render(args.out_f, tree_style = ts, dpi = 300)    



