

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from colors import default_dbf_color_map
from annotation_utils import get_orf_data, get_macisaac_data

########### some internal configurations ##############
annotate_ymin = -0.2 # the orf and other annotation y position
annotate_ymax = -0.05 # the orf and other annotation y position
ymax = 1.1 # the y limit 
#######################################################

def visualize_dbf_color_map(dbf_color_map):
    plt.figure(figsize=(10,10))
    dbfs = sorted(dbf_color_map.keys())
    ndbfs = len(dbfs)
    ys = np.arange(0, ndbfs * 20, 20);

    plt.hlines(y = ys, xmin=0, xmax=10, colors = [dbf_color_map[_] for _ in dbfs], linewidth=10);
    plt.ylim(-15, ys[-1] + 15)
    plt.yticks(ys, dbfs);

def annotate_orf(orf_start, orf_end, orf_strand, orf_name, start, end):

    arrow_length = (end - start) / 100.0
    # if both gene start and gene end are in range, draw a rectange like annotation
    if orf_start > start and orf_start < end and orf_end > start and orf_end < end:
        plt.vlines(orf_start, annotate_ymin, annotate_ymax)
        plt.vlines(orf_end, annotate_ymin, annotate_ymax)
        
        if orf_strand == '+':
            plt.arrow(orf_start, annotate_ymin, orf_end-orf_start-arrow_length, 0, head_width=0.05, 
                head_length=arrow_length, fc='black', ec='black')
        else:
            plt.arrow(orf_end, annotate_ymin, orf_start-orf_end+arrow_length, 0, head_width=0.05, 
                head_length=arrow_length, fc='black',ec='black')
        
        plt.annotate(orf_name, ((orf_start+orf_end)/2-10, annotate_ymin - 0.06))
    elif orf_start > start and orf_start < end:
        plt.vlines(orf_start, annotate_ymin, annotate_ymax)
        if orf_strand == '+':
            plt.arrow(orf_start, annotate_ymin, end - orf_start - arrow_length,0, head_width=0.05, 
                head_length=arrow_length, fc='black', ec='black')
        else:
            plt.arrow(end, annotate_ymin, orf_start - end + arrow_length,0, head_width=0.05, 
                head_length=arrow_length, fc='black', ec='black')
        
        plt.annotate(orf_name, ((orf_start+end)/2-10, annotate_ymin - 0.06))
    elif orf_end > start and orf_end < end:
        plt.vlines(orf_end, annotate_ymin, annotate_ymax)
        if orf_strand == '+':
            plt.arrow(start, annotate_ymin, orf_end - start - arrow_length, 0, head_width=0.05, 
                head_length=arrow_length, fc='black', ec='black' )
        else:
            plt.arrow(orf_end, annotate_ymin, start - orf_end + arrow_length, 0, head_width=0.05, 
                head_length=arrow_length, fc='black', ec='black')
        plt.annotate(orf_name, ((start+orf_end)/2-10, annotate_ymin - 0.06))


def plot_orf_annotation(chromo, start, end, orf_annotation):

    #orf data
    orf = get_orf_data(chromo, start, end, orf_annotation)
    for i in range(orf.shape[0]):
        orf_start = orf.start.iloc[i]
        orf_end = orf.end.iloc[i]
        orf_name = orf.name.iloc[i]
        orf_strand = orf.strand.iloc[i]

        annotate_orf(orf_start, orf_end, orf_strand, orf_name, start, end)

def plot_macisaac_annotation(chromo, start, end, dbf_color_map, macisaac_annotation):

    macisaac = get_macisaac_data(chromo, start, end, macisaac_annotation)
    for i in range(macisaac.shape[0]):
        x = (macisaac.start.iloc[i]+macisaac.end.iloc[i])/2.0
        marker = '>' if macisaac.strand.iloc[i] == '+' else '<'
        name = macisaac.name.iloc[i]
        
        if name in dbf_color_map:
            color = dbf_color_map[name]  
        elif name in default_dbf_color_map:
            color = default_dbf_color_map[name]
        else:
            color = 'black'

        plt.plot(x, annotate_ymax, marker=marker, color=color)
        plt.text(x, (annotate_ymax + annotate_ymin)/2.0, name, color=color,rotation='vertical')

def preprocess_occupancy_profile(op, drop_threshold):
    '''process occupancy profile before plotting'''

    #drop columns that do not need plotting
    drop_columns = set(['background', 'nuc_padding', 
        'nuc_center', 'nuc_start', 'nuc_end', 'Tata'])

    #drop columns that have low binding prob
    drop_columns = drop_columns.union( set(op.columns[op.apply(max) < drop_threshold]))
    
    for col in drop_columns:
        if col in op:
            op = op.drop(col, axis = 1)

    return op

def plot_dbf_binding(op, dbf_color_map):
    
    #separate nucleosome from other DBFs
    dbfs = list(op.columns.values)
    dbfs.remove('coordinate')
    if 'nucleosome' in dbfs:
        dbfs.remove('nucleosome')
        nuc_present = True
    else:
        nuc_present = False

    #plot nucleosome first
    if nuc_present:
        plt.plot(op.coordinate, op.loc[:, 'nucleosome'], color = dbf_color_map['nucleosome'], label = 'nuc')
        plt.fill_between(op.coordinate, op.loc[:, 'nucleosome'], color = dbf_color_map['nucleosome'])

    #plot all other dbfs
    for dbf in dbfs:
        plt.plot(op.coordinate, op.loc[:, dbf], color = dbf_color_map[dbf], label = dbf)
        plt.fill_between(op.coordinate, op.loc[:, dbf], color = dbf_color_map[dbf])

def read_occupancy_profile(file_name):

    return pd.read_csv(file_name, sep = '\t')

def plot_occupancy_profile(op, chromo, coordinate_start, padding = 0, threshold = 0.1, figsize=(18,6), orf_annotation = None, macisaac_annotation = None, file_name = None, dbf_color_map = default_dbf_color_map):
    
    plt.figure(figsize=figsize)

    op['coordinate'] = np.arange(coordinate_start, coordinate_start + op.shape[0])
    op = op.iloc[padding:(-padding-1), :]
    
    op = preprocess_occupancy_profile(op, drop_threshold=threshold)

    #plot DBFs using area plot
    plot_dbf_binding(op, dbf_color_map)
    
    if orf_annotation is not None:
        plot_orf_annotation(chromo, op.coordinate.iloc[0], op.coordinate.iloc[-1], orf_annotation)

    if macisaac_annotation is not None:
        plot_macisaac_annotation(chromo, op.coordinate.iloc[0], 
            op.coordinate.iloc[-1], dbf_color_map=dbf_color_map,
             macisaac_annotation=macisaac_annotation)

    #######################  set axis properties #######################
    plt.xlim(op.coordinate.iloc[0]-op.shape[0]/100.0, op.coordinate.iloc[-1] + op.shape[0]/100.0)

    if orf_annotation is not None or\
        macisaac_annotation is not None:
        # draw an annotation line
        plt.hlines(annotate_ymax, plt.xlim()[0], plt.xlim()[1], linewidth=2, color='orange')

        # set y axis limits so that the annotation can be shown properly
        plt.ylim(annotate_ymin - 0.15, 1)

    else:
        plt.ylim(0, 1)


    # legend
    leg = plt.legend(loc='lower left', ncol = 12, bbox_to_anchor = (0., 1.),
        borderaxespad=0, frameon=False, framealpha=0)
    
    for legobj in leg.legendHandles:
        legobj.set_linewidth(10.0)

    #set the y axis bounds so only the 0-1.0 part is shown
    plt.axes().spines['left'].set_bounds(0, 1.0)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1]);

    # delete the x axis, but the ticks will remain
    plt.axes().spines['bottom'].set_linewidth(0)

    # remove top and right borders   
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    plt.xlabel('Chr %d' % chromo)

    if file_name is None:
        plt.show()
    else:
        plt.savefig(file_name)


