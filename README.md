# ShNeigh
ShNeigh code

[posMatrix,retVals]=ShNeigh(heatmapMatrix,method_type)

Input:
heatmapMatrix: heatmap matrix
method_type: 1 for ShNeigh1; 2 for ShNeigh2

Output:
posMatrix: N*4 matrix, the first column is a vector that marks valid genome loci, columns 2-4 are coordinates
retVal: a structure that contains some resulting information
