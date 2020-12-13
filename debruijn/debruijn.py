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
from itertools import chain
import statistics
from collections import Counter
from random import randint
import random
import networkx as nx
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
    """
    Read a fastq file
    :param fastq_file: fastq file
    :return:
    """
    with open(fastq_file, 'rt') as file:
        for _ in file:
            yield next(file).strip()
            next(file)
            next(file)


def cut_kmer(read, kmer_size):
    """
    Take a read and return kmer cut in given size
    :param read: sequence
    :param kmer_size: size the kmer needs to be cut
    :return:
    """
    for i in range(len(read)-kmer_size+1):
        yield read[i:kmer_size+i]


def build_kmer_dict(fastq_file, kmer_size):
    """
    Read sequence from a fastq file then create a dictionary with each kmer and its occurences
    :param fastq_file: fastq file
    :param kmer_size: size to cut the kmer
    :return: Dict of kmer and number of occurences
    """
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
    """
    Create graph of from given kmer dict
    :param kmer_dict: Dict of kmer and number of occurences
    :return: Directed Graph
    """
    graph = nx.DiGraph()
    for kmer in kmer_dict.keys():
        graph.add_edge(kmer[:-1], kmer[1:], weight=kmer_dict[kmer])
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """
    Remove given paths from the graph
    :param graph: Directed graph
    :param path_list: Paths to remove from graph
    :param delete_entry_node: If true will remove entry node as well
    :param delete_sink_node: If true will remove sink node as well
    :return: cleaned graph
    """
    for path in path_list:
        path = list(path)
        for index, node in enumerate(path):
            if index == 0:
                if delete_entry_node:
                    graph.remove_node(node)
                else:
                    continue
            elif index == len(path)-1:
                if delete_sink_node:
                    graph.remove_node(node)
                else:
                    continue
            else:
                graph.remove_node(node)
    return graph


def std(data):
    """

    :param data: List of values
    :return: Standard deviation of data
    """
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """
    Determine best path from a list of paths depending on their weights, length
     and if still tied, randomness
    :param graph: Directed graph
    :param path_list: List of paths
    :param path_length: Corresponding length of the paths
    :param weight_avg_list: Corresponding weights of the paths
    :param delete_entry_node: If true will remove entry node as well
    :param delete_sink_node: If true will remove sink node as well
    :return: Best path
    """
    max_weight = max(weight_avg_list)
    max_weight_idx = [i for i, j in enumerate(weight_avg_list) if j == max_weight]
    if len(max_weight_idx) == 1:
        path_list.pop(max_weight_idx[0])
        graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    else:
        max_length = max(path_length)
        max_length_idx = [i for i, j in enumerate(path_length) if j == max_length]
        if len(max_length_idx) == 1:
            path_list.pop(max_length_idx[0])
            graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
        else:
            random_index = randint(0, len(path_list)-1)
            path_list.pop(random_index)
            graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph, path):
    """
    Give average weight of a path
    :param graph: Directed graph
    :param path: A path in the graph
    :return: Average weight of the path
    """
    weight_list = []
    for i in range(len(path)-1):
        weight_list.append(graph[path[i]][path[i+1]]['weight'])
    return statistics.mean(weight_list)


def solve_bubble(graph, ancestor_node, descendant_node):
    """
    Remove bubbles between two nodes
    :param graph: Directed graph
    :param ancestor_node: Ancestor node
    :param descendant_node: Descendant node
    :return: Graph without the bubbles
    """
    paths = nx.all_simple_paths(graph, ancestor_node, descendant_node)
    path_list = []
    average_weight_list = []
    path_length = []
    for path in paths:
        path_list.append(tuple(path))
        average_weight_list.append(path_average_weight(graph, path))
        path_length.append(len(path))
    graph = select_best_path(graph, path_list, path_length, average_weight_list,
                             delete_entry_node=False, delete_sink_node=False)
    return graph


def simplify_bubbles(graph):
    """
    Clean graph of all bubbles
    :param graph: Directed Graph
    :return: Cleaned graph
    """
    nodes = graph.nodes
    bubbles = []
    for node in nodes:
        predecessors = graph.predecessors(node)
        pred = []
        for predecessor in predecessors:
            predecessors_of_predecessors = graph.predecessors(predecessor)
            for predecessor_of_predecessors in predecessors_of_predecessors:
                pred.append(predecessor_of_predecessors)
            counter = Counter(pred)
            for node_pred in counter.keys():
                if counter[node_pred] > 1:
                    bubbles.append((node_pred, node))
    for bubble in bubbles:
        graph = solve_bubble(graph, bubble[0], bubble[1])
    return graph

def solve_entry_tips(graph, starting_nodes):
    """
    Clear Graph of bad entry nodes
    :param graph: Directed Graph
    :param starting_nodes: Entry nodes
    :return: Graph with no bad entry nodes
    """
    paths = []
    for node in starting_nodes:
        for sink in get_sink_nodes(graph):
            paths.append(nx.all_simple_paths(graph, node, sink))
    all_paths = []
    for path_gen in paths:
        for path in path_gen:
            all_paths.append(path)
    counts = Counter(chain(*map(set, all_paths)))
    tips = [[i for i in sublist if counts[i] == 1] for sublist in all_paths]
    longest_path = [x for x in tips if len(x) == max(map(len, tips))]
    if len(longest_path) == 1:
        graph.remove_nodes_from([x for x in tips if x not in longest_path][0])
    else:
        weight_list = []
        for index, path in enumerate(tips):
            path.append(all_paths[index][len(path)])
            weight_list.append(path_average_weight(graph, path))
        min_index, _ = min(enumerate(weight_list), key=itemgetter(1))
        graph.remove_nodes_from(tips[min_index][:-1])
    return graph


def solve_out_tips(graph, ending_nodes):
    """
     Clear Graph of bad ending nodes
     :param graph: Directed Graph
     :param ending_nodes: Ending nodes
     :return: Graph with no bad ending nodes
     """
    paths = []
    for node in ending_nodes:
        for start in get_starting_nodes(graph):
            paths.append(nx.all_simple_paths(graph, start, node))
    all_paths = []
    for path_gen in paths:
        for path in path_gen:
            all_paths.append(path)
    counts = Counter(chain(*map(set, all_paths)))
    tips = [[i for i in sublist if counts[i] == 1] for sublist in all_paths]
    longest_path = [x for x in tips if len(x) == max(map(len, tips))]
    if len(longest_path) == 1:
        graph.remove_nodes_from([x for x in tips if x not in longest_path][0])
    else:
        weight_list = []
        for index, path in enumerate(tips):
            path.append(all_paths[index][-len(path)-1])
            path.reverse()
            weight_list.append(path_average_weight(graph, path))
        min_index, _ = min(enumerate(weight_list), key=itemgetter(1))
        graph.remove_nodes_from(tips[min_index][1:])
    return graph

def get_starting_nodes(graph):
    """
    Get all starting nodes from a graph
    :param graph: Directed graph
    :return: Starting nodes
    """
    nodes = graph.nodes()
    starting_nodes = []
    for node in nodes:
        if list(nx.edge_dfs(graph, node, orientation='reverse')):
            continue
        else:
            starting_nodes.append(node)
    return starting_nodes


def get_sink_nodes(graph):
    """
    Get all ending nodes
    :param graph: Directed graph
    :return: Ending nodes
    """
    nodes = graph.nodes()
    sink_nodes = []
    for node in nodes:
        if list(nx.dfs_edges(graph, node)):
            continue
        else:
            sink_nodes.append(node)
    return sink_nodes


def get_contigs(graph, starting_nodes, ending_nodes):
    """
    Get all contigs and their length from a set of starting nodes to ending nodes
    :param graph: Directed graph
    :param starting_nodes: Starting nodes
    :param ending_nodes: Ending nodes
    :return: List of tuples of contigs and their length
    """
    contigs = []
    for starting_node in starting_nodes:
        for ending_node in ending_nodes:
            for path in nx.all_simple_paths(graph, source=starting_node, target=ending_node):
                full_path = path[0]
                for node in path[1:]:
                    full_path += node[1:]
                contigs.append((full_path, len(full_path)))
    return contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i + width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):
    """
    Save contigs in a file
    :param contigs_list: Contigs and their length
    :param output_file: file name
    :return: File with contigs with fasta format
    """
    with open(output_file, 'wt') as file:
        for index, contig in enumerate(contigs_list):
            file.write('>contig_{0} len={1}\n'.format(index, contig[1]))
            seq = fill(contig[0])
            file.write(seq+'\n')


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Create graph
    print('Creating the graph ...')
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)

    # Solve bubbles
    print('Solving bubbles ...')
    graph = simplify_bubbles(graph)

    # Solve tips
    print('Solving tips ...')
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)

    # Get and save contigs
    print('Saving contigs ...')
    contigs_list = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs_list, args.output_file)
    print('Done')




if __name__ == '__main__':
    main()
