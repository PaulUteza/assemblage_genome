#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import re
from operator import itemgetter
import statistics
from random import randint
import random
import networkx as nx
import matplotlib
random.seed(9001)

__author__ = "UTEZA Paul"
__copyright__ = "CY TECH"
__credits__ = ["UTEZA Paul"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "UTEZA Paul"
__email__ = "utezapaul@eisti.eu"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file, 'rt') as file:
        for _ in file:
            yield next(file).strip()
            next(file)
            next(file)


def cut_kmer(read, kmer_size):
    for i in range(len(read)-kmer_size+1):
        yield read[i:kmer_size+i]


def build_kmer_dict(fastq_file, kmer_size):
    reads = read_fastq(fastq_file)
    kmer_dict = {}
    for read in reads:
        kmers = cut_kmer(read, kmer_size)
        for kmer in kmers:
            if kmer in kmer_dict.keys():
                continue
            else:
                # Find and count all overlapping matches
                kmer_dict[kmer] = len(re.findall('(?={})'.format(re.escape(kmer)), read))
    return kmer_dict


def build_graph(kmer_dict):
    kmer_cut_dict = {}
    for kmer in kmer_dict.keys():
        kmer_cut_dict[kmer[:-1]] = {kmer[1:]: {'weight': kmer_dict[kmer]}}
    graph = nx.DiGraph(kmer_cut_dict)
    return graph

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass

def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    pass

def get_sink_nodes(graph):
    pass

def get_contigs(graph, starting_nodes, ending_nodes):
    pass

def save_contigs(contigs_list, output_file):
    pass

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()




if __name__ == '__main__':
    main()
