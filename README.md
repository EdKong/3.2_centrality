# CBB752_3.2_centrality
Given a MITAB gene network interaction file, find the degree and betweenness centralities of each node

Usage:      python Centrality.py -i <input .txt file in MI TAB 2.5 format> -o <output .csv file>

Example:    python Centrality.py -i exampleMITAB.txt -o centrality_out.csv

Note:       Outputs a csv. with normalized centrality and betweenness for each gene as columns
(Genes identified by their DIP identifiers in each case)
