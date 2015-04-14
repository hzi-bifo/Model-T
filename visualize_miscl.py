import pandas as ps
import numpy as np
import ete2
from ete2 import TreeStyle, faces, AttrFace, CircleFace
import os
import colorsys
miscl_phypat = ps.read_csv("misclassified_extd_phypat.tsv", index_col = 0, sep = "\t")
print miscl_phypat
miscl_phypat_ml = ps.read_csv("misclassified_extd_phypat+ml.tsv", index_col = 0, sep = "\t")
t = ete2.Tree("gideon_tree.nwck")
t.scientific_name = "Bacteria"
t.rank = "Kingdom"
bacc_max = max(miscl_phypat.loc[:, "bacc"].max(), miscl_phypat_ml.loc[:, "bacc"].max())
bacc_min = min(miscl_phypat.loc[:, "bacc"].min(), miscl_phypat_ml.loc[:, "bacc"].min())

#prune tree to genus level
#get genus nodes
genera = []
for i in t.traverse():
    if i.rank == "genus":
        genera.append(i)
t.prune(genera)
#join scientific name and miscl statistics
for n in t.traverse():
    n.miscl_full = "%s | %i %i %i %i %s | %i %i %i %i %s" % tuple([n.scientific_name] + [round(miscl_phypat_ml.loc[n.name,].iloc[i], 2) for i in range(5)] + [round(miscl_phypat.loc[n.name,].iloc[i], 2) for i in range(5)])

def map_bacc2color(val, minval, maxval):
    # convert val in range minval..maxval to the range 0..120 degrees which
    # correspond to the colors red..green in the HSV colorspace
    h = (float(val-minval) / (maxval-minval)) * 120
    # convert hsv color (h,1,1) to its rgb equivalent
    # note: the hsv_to_rgb() function expects h to be in the range 0..1 not 0..360
    return  '#%02x%02x%02x' %tuple(map(lambda x: int(x*255), colorsys.hsv_to_rgb(h/360, 1., 1.)))


def miscl_layout(node):
    if node.is_leaf():
        name_face = AttrFace("miscl_full")
        faces.add_face_to_node(name_face, node, column=0, position="branch-right")
    else:
        name_face = AttrFace("miscl_full", fsize = 10)
        faces.add_face_to_node(name_face, node, column=0, position="branch-right")
    c_mp = CircleFace( radius = 10, style="sphere", color = map_bacc2color(miscl_phypat_ml.loc[node.name, "bacc"], bacc_min, bacc_max))
    c_mp.opacity = 0.7
    c_p = CircleFace( radius = 10, style="sphere", color = map_bacc2color(miscl_phypat.loc[node.name, "bacc"], bacc_min, bacc_max))
    c_p.opacity = 0.7
 
    faces.add_face_to_node(c_mp, node, column=1, position="float")
    faces.add_face_to_node(c_p, node, column=2, position="float")

 

ts = TreeStyle()
#ts.rotation = 90
ts.layout_fn = miscl_layout
ts.show_leaf_name = False
t.render("miscl_tree.png", tree_style = ts)    
