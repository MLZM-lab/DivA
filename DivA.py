#! /usr/bin/env python

##Import libraries
from modules import *

####################
###Get the options

parser = argparse.ArgumentParser(description='Identify very divergent potentially non-homologous windows in a protein multiple sequence alignment.')

parser.add_argument('alnNamesFile', help='A txt file with the file name(s) of the MSA(s) on which to perform the method')
parser.add_argument("--mask", help='Flag for the output of an alignment with the wrong windows masked with XXs [default not set]', action="store_true")
parser.add_argument("--printAllwindows", help='Flag for the output of a file with the parameter values and start and end positions of all the windows in the MSA(s) [default not set]', action="store_true")
parser.add_argument("-w", type=int, default=12, help='The size of the sliding window [default 12]')
parser.add_argument("-g", type=float, default=0.6, help='Maximum gap content in a window to be considered [default 0.6]')
parser.add_argument("-p", type=int, default=1, help='The number of standard deviations from the mean of the alpha parameter to use as threshold [default 1]')
parser.add_argument("-zp", type=int, default=2, help='The number of standard deviations from the mean of the Zalpha parameter to use as threshold [default 2]')
parser.add_argument("-d", type=int, default=2, help='The number of standard deviations from the mean of the beta parameter to use as threshold  [default 2]')
parser.add_argument("-zd", type=int, default=2, help='The number of standard deviations from the mean of the Zbeta parameter to use as threshold [default 2]')
parser.add_argument("-o", default="out", help='Output basename prefix [default "out"]')
parser.add_argument("-m", default="blosum62.txt", help='The amino acid distance matrix [default "blosum62.txt"]')

args = parser.parse_args()

####################
##Initialize variables

alnNamesFile= args.alnNamesFile
window= args.w
maxGapContent= args.g
outprefix= args.o
aaMatrix= args.m
Psd= args.p
Zsd= args.zp
Dsd= args.d
Z62sd= args.zd


if args.printAllwindows:
	out2=open("%s_PredictorVariables.txt"%(outprefix), "w")

alnNames=open(alnNamesFile)
alnNames=alnNames.readlines()
matrixAA=readMatrix(aaMatrix)
step=1

#Initialize dict
aaFreqsDict={}
aminoacids=["X","A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "-"] # amino acid frequencies, gaps and X will be given a zero prob
for element in aminoacids:
	aaFreqsDict[element]=[] 

windowsMaster={}
wrongMaster={} #This is going to be the windows from {windows} that have 1 in the [classification][i]
#################################################### Main body

########## Get and trim the windows and associated parameters values

for filename in alnNames:
	windows= {} #I had to make a dictionary with all the output with the values to then calculate the thresholds. It's: windows[name]=[ [filename], [start],[end],[P],[Z],[D],[Z62],[classificacion] ]
	filename=filename.rstrip()
	print("I'm in %s"%(filename))   
	align = AlignIO.read(filename, "fasta")
        j=0
        while j<len(align[1].seq)-window:
                windowAlign=align[:, j:j+window]
                scoresProb=getAveProbPerWindowPerSeq(windowAlign,maxGapContent,aaFreqsDict)
                scoresMatrixPerCol=getColumnDistScoreToClosest(windowAlign, matrixAA)
                #here is where it goes species per species testing the window values        
                for record in align: 
                        name=record.id
                        P=scoresProb[name][0]
			if P!="NA" and float(P)<0.7: #Now I don't want the windows that would result in NAs which are given if P==NA or float(P)>0.7 (too conserved)
				start=j+1
				end=j+window
                                firstAA,lastAA=getFirstLastAA(windowAlign,aaFreqsDict,name)
				if ( (lastAA!="NA" and lastAA>0.9) or (firstAA!="NA" and firstAA>0.9) ): #Check if I need to trim and recalculate
	                                windowAlign, start, end= trimWindowEdges(windowAlign, align, aaFreqsDict,name, start, end, firstAA,lastAA ) #I turned this part into this function
	                                scoresProb=getAveProbPerWindowPerSeq(windowAlign, maxGapContent,aaFreqsDict)
        	                        scoresMatrixPerCol=getColumnDistScoreToClosest(windowAlign, matrixAA)
					P=scoresProb[name][0]
                                        if P!="NA" and float(P)<0.7: #If after the trimming it's still a useful window
						windows=fillParametersLists(scoresProb, scoresMatrixPerCol, name, P, windows, start, end)						
				else:
					windows=fillParametersLists(scoresProb, scoresMatrixPerCol, name, P, windows, start, end)
		j+=step
	if len(windows)>0: #fill windowsMaster which has the windows of all the filenames
		if not windowsMaster.has_key(filename): 
			windowsMaster[filename]=windows #windows has [SpeciesName]=start, end, P, Z, D, Z62


########## Calculate thresholds
print "Calculating the parameters thresholds"
Pthresh,Zthresh,Dthresh,Z62thresh=defineThresholds(windowsMaster, Psd, Zsd, Dsd, Z62sd)
### Classify windows
print "Classifying the windows"
windowsMaster=classifyWindows(windowsMaster, Pthresh , Zthresh , Dthresh, Z62thresh )
### Merge windows
print "Merging the windows"
wrongMaster=getAndMergeWrongWindows(windowsMaster, wrongMaster)
### Print wrongMaster 

out1=open("%s_WrongWindows.txt"%(outprefix), "w")

#print header
out1.write("alignment\tsequence\tstart\tend\tP\tZP\tD\tZD\n")

if len(wrongMaster)>0:
	OutputDict(out1, wrongMaster) #Print filename per filename
	### Optional alignment masking
	if args.mask:
		print "Masking alignment"
		OptmaskAln( wrongMaster )
### Print optionally the windowsMaster dict
if args.printAllwindows:
	OutputDict(out2, windowsMaster) #Print filename per filename

out1.close()
if args.printAllwindows:
	out2.close()

