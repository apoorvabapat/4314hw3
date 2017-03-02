#!/usr/bin/env python
from sys import *
import numpy

#################################################
#   Authors: Tyler Glotz and Apoorva Bapat
#################################################
# Algorithm Outline:
#
# Python implementation of Needleman-Wunsch to achieve global alignment
# Algorithm runs in time complexity of O(mn) and a space complexity of O(mn).
#
#
# needlemanW Process:
# 1. Create matrix
# 2. Fill matrix with distance scores based on given values (0 for match, 1 for gap, 1 for mismatch)
# 3. Traceback starting a bottom right of matrix and work way backwards, taking least costly move (min)
# 4. Add gaps if neccessary
# 5. Reverse strings back, and output distance matrix

##################################################


#TODO: Handle FASTA file input
#TODO: For each sequence, run against every other sequence in fasta
#TODO: print each global alignment and corresponding distance matrix

seq1='GAGACTAGA'
seq2='AGACATCGA'

def needlemanW(seq1, seq2):
   
 m = len(seq1) + 1
 n = len(seq2) + 1  # length of two sequences
    
    # Generate DP table of zeros using numpy
 distance = numpy.zeros(shape=(m,n))  # matrix
   
    #Dynamic Programming calculations
 for i in range(0, m):
  distance[i][0] = 1 * i
 for j in range(0, n):
  distance[0][j] = 1 * j
 for i in range(1, m):
  for j in range(1, n):
   #Here we check to see if indices match, have gap, or mismatch
   if seq1[i-1] == seq2[j-1]:
    match = distance[i - 1][j - 1]  #if matched adds no penalty score
   elif seq1[i-1] == '-' or seq2[j-1] == '-':
    match = distance[i - 1][j - 1] + 1 #penalized for gaps +1
   else:
    match = distance[i - 1][j - 1] + 1 #penalized for mismatch +1
   
   delete = distance[i - 1][j] + 1 # +1 indicates penalty
   insert = distance[i][j - 1] + 1 # +1 indicates penalty
   
   distance[i][j] = min(match, delete, insert) #take min to find least costly move


    # Traceback
 alignment1 = '' 
 alignment2 = '' #initialize alignment strings
 i = m-1
 j = n-1 # position i and j to bottom right index of matrix, so we can work our way up
 
 while i > 0 and j > 0: # while we are not at top left index of matrix
  
  #directional names given indices to make code clearer
  current = distance[i][j]
  diagonal = distance[i-1][j-1]
  up = distance[i][j-1]
  left = distance[i-1][j]
  
  if seq1[i-1] == seq2[j-1]:
   a = 0 #matching score
  elif seq1[i-1] == '-' or seq2[j-1] == '-':
   a = 1 #gap penalty
  else:
   a = 1 #gap penalty
  
  if current == diagonal + a: #a is either nothing in event of match, or 1 in case of gap or mismatch
   alignment1 += seq1[i-1]
   alignment2 += seq2[j-1]
   i -= 1
   j -= 1
  elif current == left + 1:
   alignment1 += seq1[i-1]
   alignment2 += '-'
   i -= 1
  elif current == up + 1:
   alignment1 += '-'
   alignment2 += seq2[j-1]
   j -= 1

  #tracing up to the top left cell, 
 while i > 0:
  alignment1 += seq1[i-1]
  alignment2 += '-'
  i -= 1
 while j > 0:
  alignment1 += '-'
  alignment2 += seq2[j-1]
  j -= 1

    
 if alignment1 == alignment2 or alignment1 > alignment2:
   la = len(alignment1)
 else:
   la = len(alignment2)

 #all the print outs
 print '\n' + alignment1[::-1] #reverse alignment1 back to finalize the alignment
 print alignment2[::-1] #reverse alignment2 back to finalize the alignment
 print '=' * la + ' (',la,')' #outputs '=' * length of longest string + length
 print '\n' , distance , '\n' #outputs distance matrix
        








print(needlemanW(seq1, seq2))