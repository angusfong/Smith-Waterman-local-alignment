#!/usr/bin/python
__author__ = "Angus Fong"
__email__ = "hoching.fong@yale.edu"
__copyright__ = "Copyright 2019"
__license__ = "GPL"
__version__ = "1.0.0"
### Usage: python hw1.py -i <input file> -s <score file>
### Example: python hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm
### Scripting must be done from scratch, without the use of any pre-existing packages.
### Python standard library (I/O) and numpy are allowed.
import argparse
import numpy as np
### This is one way to read in arguments in Python.

### We need to read input file and score file.
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
args = parser.parse_args()

### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap):	

	### Print input and score file names. You can comment these out.

	# create an output file
	f = open('output.txt', 'w')

	#read seqs
	seqs = [line.rstrip() for line in open(inputFile)]

	f.write('Sequences:\n')
	f.write(seqs[0] + '\n')
	f.write(seqs[1] + '\n\n')

	#read sim
	sim_file = [line.rstrip().split() for line in open(scoreFile)]
	sim_mat = [s[1:] for s in sim_file if len(s) != 0][1:]
	keys = sim_file[0]

	#store sim as dictionary to make lookup quick
	sim_dict = {keys[i]:{keys[j]:float(sim_mat[i][j]) for j in range(len(keys))} for i in range(len(keys))}

	#affine gap penalty function
	def g(k): #g(k) = openGap + k * extGap
		return openGap + (k-1) * extGap

	#init DP matrices
	#A[i,j] has score for optimal alignment ending with no gaps at the end
	A = [[0 for j in range(len(seqs[0])+1)] for i in range(len(seqs[1])+1)]
	bk = [[None for j in range(len(seqs[0])+1)] for i in range(len(seqs[1])+1)]
	
	#write DP matrices
	for i in range(1,len(seqs[1])+1):
		for j in range(1,len(seqs[0])+1):
			# cell one diagonal above
			sim = sim_dict[seqs[1][i-1]][seqs[0][j-1]]
			diag = A[i-1][j-1] + sim
			diag_coord = (i-1, j-1)
			
			# iterate through all cells to the left, find max score
			lefts = [A[i][j-l] + g(l) for l in range(1,j+1)]
			left = max(lefts)
			left_ind = lefts.index(left)
			left_coord = (i, j - left_ind - 1)

			# iterate through all cells above, find max score
			ups = [A[i-l][j] + g(l) for l in range(1,i+1)]
			up = max(ups)
			up_ind = ups.index(up) 
			up_coord = (i - up_ind - 1, j)

			# build array of candidate scores and sources
			scores = [diag, up, left, 0]
			sources = [diag_coord, up_coord, left_coord, None]

			# write to DP and backtracing matrix
			A[i][j] = int(max(scores))
			bk[i][j] = sources[scores.index(max(scores))]

	#print(A); print(bk)

	# save A to DP.txt, adding some rownames and column names to A
	colnames = [''] + [c for c in seqs[0]]
	rownames = ['',''] + [c for c in seqs[1]]
	DP = [colnames] + A
	DP = [[rownames[row_ind]] + DP[row_ind] for row_ind in range(len(DP))]
	f.write('DP matrix:\n\n')
	np.savetxt(f, DP, fmt="%s",delimiter='\t')

	# now do backtracking
	best_score = max(map(max, A))
	f.write('\nBest score = ' + str(best_score) + '\n\n')
	best_coords = [(i,j) for i in range(len(A)) for j in range(len(A[0])) if A[i][j] == best_score]

	aligned_seq0s = []; aligned_seq1s = []

	#function takes one starting coordinate and backtracks, in case there are ties between max scores
	def backtrack(start_coord): 
		i, j = start_coord
		seq0 = [seqs[1][i-1]]; seq1 = [seqs[0][j-1]]
		while A[i][j] > 0:
			new_i, new_j = bk[i][j]
			if A[new_i][new_j] == 0: #can break if zero at backtrack location
				break
			#4 possibilities
			if new_i is None: #shouldn't happen
				raise Exception('None where there should be a source')
			elif new_i == i-1 and new_j == j-1: #diagonal
				seq0.append(seqs[1][new_i-1]); seq1.append(seqs[0][new_j-1])
			elif new_i == i: #left
				seq0.extend(['-']*(j-new_j)); seq1.extend([seqs[0][x] for x in range(j-1,new_j-1,-1)])
			else: #up
				seq1.extend(['-']*(i-new_i)); seq0.extend([seqs[1][y] for y in range(i-1,new_i-1,-1)])
			i, j = new_i, new_j
		aligned_seq0s.append(seq0); aligned_seq1s.append(seq1)

	for start_coord in best_coords:
		backtrack(start_coord)

	# render alignment		
	for seq_id in range(len(aligned_seq0s)):
		f.write('Local alignment ' + str(seq_id+1) + ':\n')
		seq0, seq1 = aligned_seq0s[seq_id], aligned_seq1s[seq_id]
		match = [' ']*len(seq0)
		match_inds = [i for i in range(len(seq0)) if seq0[i]==seq1[i]]
		for i in match_inds:
			match[i] = '|'
		seq0.reverse(); match.reverse(); seq1.reverse()
		f.write(' '.join(seq0) + '\n')
		f.write(' '.join(match) + '\n')
		f.write(' '.join(seq1) + '\n')
		f.write('\n')

	f.close()

	print('Smith-Waterman algorithm complete. See output.txt for output.')

### Run your Smith-Waterman Algorithm
runSW(args.input, args.score, args.opengap, args.extgap)