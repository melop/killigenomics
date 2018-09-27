#!/usr/bin/python


import os
import sys
import argparse
import time
from ete3 import Tree, NodeStyle, TreeStyle, TextFace

##########
# inputs #
##########
start_time = time.time()

### Option defining
parser = argparse.ArgumentParser(prog="pcoc_taxa2nodenum.py",
                                 description='')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

##############
requiredOptions = parser.add_argument_group('Required arguments')
requiredOptions.add_argument('-t', "--tree", type=str,
                             help='input tree name', required=True)
requiredOptions.add_argument('-L', "--tiplabels", type=str,
                             help='the tip labels', required=True)

##############

### Option parsing
args = parser.parse_args()


def init_tree(nf):
    t = Tree(nf)

    #Alternatively we could read a tree from a file into a string "line", and then use:
    # t =  Tree( line )

    nodeId = 0
    for n in t.traverse("postorder"):
        n.add_features(ND=nodeId)
        nodeId = nodeId + 1

    return t

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))


if not os.path.isfile(args.tree):
    print ("Error %s does not exist" %args.tree)
    sys.exit(1)

oTree = init_tree(args.tree)
arrEvents = args.tiplabels.split('/');
arrOutAll = [];

for sEvent in arrEvents:
    arrTips = sEvent.split(',');
    arrTreeLeaves = oTree.get_leaf_names();

    arrTipsFound = intersection(arrTips ,arrTreeLeaves );

    if len(arrTipsFound) < 2:
        print ("Error no tips found" )
        sys.exit(1)

    #print(arrTips);
    #print(intersection(arrTips ,arrTreeLeaves ));

    oMonophyly = oTree.check_monophyly(values=arrTipsFound, target_attr="name")
    if not oMonophyly[0]:
        print ("Error taxa not monophyletic\n" )
        sys.exit(1)

    oAnc = oTree.get_common_ancestor(arrTipsFound);
    #print(oAnc.ND);

    arrOut = [str(oAnc.ND)];

    for oNode in oAnc.iter_descendants("postorder"):
      # Do some analysis on node
      arrOut.append(str(oNode.ND))

    
    arrOutAll.append(','.join(arrOut))

print('/'.join(arrOutAll));


