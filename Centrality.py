#!/usr/bin/python

from __future__ import division
import argparse
import numpy
import csv
import re
import brandes

__author__ = "Edward Kong"
__copyright__ = "Copyright 2016"
__credits__ = ["ELK"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Edward Kong"
__email__ = "edward.kong@yale.edu"

# Usage:      python Centrality.py -i <input .txt file in MI TAB 2.5 format>
# Example:    python Centrality.py -i exampleMITAB.txt -o centrality_out.csv
# Note:       Outputs a csv. with normalized centrality and betweenness for each gene
#               (Genes identified by their DIP identifiers in each case)


# Set up argument parser
parser = argparse.ArgumentParser(description='CentralityCalculator')
parser.add_argument('-i', '--input', help='input file, tab delimited in a .txt file, MI TAB 2.5 format', required=True)
parser.add_argument('-o', '--output', help='output file, csv format', required=False, default='centrality_output.csv')
args = parser.parse_args()


def runCentral(inputFile, outputFile):

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Degree Centrality
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # create a dictionary to hold the unique genes (keys) as well as the degree of the node (values)
    # this will be used to calculate degree centrality
    gene_degree = {}
    # this variable holds the interaction data (pairs of gene names)
    interactions = []

    # Open the input file
    with open(inputFile,'r') as f:
        # skip headings
        next(f)
        reader = csv.reader(f, delimiter='\t')
        # set up regular expression to capture gene names in DIP format
        reg = re.compile('DIP-[\w]+')

        # loop over the rows, and the two genes in each row
        for i, row in enumerate(reader):
            # append dummy names to the list of interactions
            interactions.append(['geneA', 'geneB'])
            for r in range(0, 2):

                # Extract DIP name from each gene
                gene = reg.findall(row[r])[0]

                # Add the gene to the interactions list
                interactions[i][r] = gene

                # only proceed if the interaction is not the same gene with itself
                if reg.findall(row[0])[0] != reg.findall(row[1])[0]:
                    # Add gene to dictionary if not already added
                    if gene not in gene_degree:
                        # add the gene to dictionary and give it 1 connection
                        gene_degree[gene] = 1
                    else:
                        # if gene exists in dictionary, increment num of connections -
                        gene_degree[gene] += 1
                # if the interaction is the same gene with itself, add to dictionary but with 0
                else:
                    # Add gene to dictionary if not already added
                    if gene not in gene_degree:
                        # add the gene to dictionary and give it 0 connections
                        gene_degree[gene] = 0
                    else:
                        # if gene exists in dictionary, do nothing
                        pass

    # n := total number of genes
    n = len(gene_degree)

    # Normalize the degree centrality (divide by n-1)
    for gene, degree in gene_degree.iteritems():
        # normalized degree centrality divides the degree of a node by the total number of nodes - 1
        gene_degree[gene] = degree / (n - 1)

    # Other considerations (which we do not explicitly calculate):
    # We can consider the amount of variation in degree centrality in the network
    #   To assess this, we can calculate various measures of variation. E.g. Variance/stdev or Freeman's general formula
    #   for centralization: C = [(n-1)*(n-2)]^(-1) * sum(over nodes) of (max_centrality - centrality(node i))

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Betweenness Centrality
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # create a dictionary A that contains the genes as keys, and a list of interaction partners as values,
    # then use brandes.py to generate the betweenness centrality
    # geneset := list of the unique genes in the "genes" dictionary. This list will be used to construct the adjacency
    #            matrix and implement brandes
    geneset = gene_degree.keys()

    # initialize the dictionary (make a copy)
    interaction_dict = gene_degree.copy()
    for key in interaction_dict.iterkeys():
        interaction_dict[key] = []

    # loop through the interactions list, filling in the dictionary for each interaction
    for pair in interactions:
        # index of genes A and B in the interaction pair
        geneA = pair[0]
        geneB = pair[1]
        # fill in the dictionary (symmetric)
        interaction_dict[geneA].append(geneB)
        interaction_dict[geneB].append(geneA)

    # use the brandes function in the brandes module to get the betweenness centrality
    gene_betweenness = brandes.brandes(geneset, interaction_dict)

    # normalize the betweenness (for undirected graphs, this amounts to dividing by (n-1)(n-2)*0.5. However, since
    # Brandes counts each pair of genes twice, we only divide by (n-1)(n-2)
    for key in gene_betweenness.iterkeys():
        gene_betweenness[key] /= ((n-1)*(n-2))

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # OUTPUT
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # prepare output file (write csv headers):
    with open(outputFile, 'wb') as out:
        headings = ['GeneID', 'DegreeCentrality', 'Betweenness']
        writer = csv.writer(out, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(headings)

    # write output to file
    with open(outputFile, 'ab') as out:
        writer = csv.writer(out, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        # loop through geneset, which gives keys for both degree and betweenness dictionaries
        for key in geneset:
            writer.writerow([key, gene_degree[key], gene_betweenness[key]])

# Run the centrality algorithm
runCentral(args.input, args.output)
