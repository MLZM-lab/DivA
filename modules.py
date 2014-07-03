################ Import packages

from numpy import *
from Bio import AlignIO
import re
import os
import sys
import argparse

################ Define the functions

### output the dicts

def OutputDict(out, dictMaster): 
	for filename in dictMaster:
		for name in dictMaster[filename]:
			for i in range(len(dictMaster[filename][name][0])):
				a="\t".join(["%s" % dictMaster[filename][name][el][i] for el in range(len(dictMaster[filename][name])) ])
				out.write("%s\t%s\t%s\n"%(filename,name,a))



##### Optiona alignment masking

def OptmaskAln(wrongMaster): # WrongMaster has wrongs, which have start (human start, no start base 0), end (human start, no start base 0), P, Z, D, Z62
	for filename in wrongMaster:
		aln=open(filename)
		tempAln=aln.readlines()
		aln.close()
		filenameShort=re.sub(".fasta$","", re.sub(".fas$","", re.sub(".fna$","", filename)))
		out1=open("%s_maskedAln.fasta"%(filenameShort), "w")
		for name in wrongMaster[filename]:	
			i=0
			for row in tempAln: #Find the seq of that >name
				if re.search(name, row): 
					pos=i+1
					break
				else:
					i+=1
			start=wrongMaster[filename][name][0]
			end=wrongMaster[filename][name][1]
			i=0
			while i<len(start): #Maske all the windows   ## tempAln[pos][start[i]-1:end[i]]
				tempAln[pos]=  tempAln[pos][:start[i]-1]  +   "X"*(end[i]-start[i]+1)  + tempAln[pos][end[i]:]	
				i+=1
		#Write the masked tempAln
		for row in tempAln:
			out1.write(row)
		out1.close()



######## merge and return the wrong windows

def getAndMergeWrongWindows(windowsMaster, wrongMaster):
	for filename in windowsMaster:
		wrong={}
		for key in windowsMaster[filename]:
			k=0
			i=[]
			i = [index for index, x in enumerate(windowsMaster[filename][key][6]) if x==1 ] #Keep only the ones that are wrong 
			if len(i)>0: #initialize it if it does have wrong windows
				if not wrong.has_key(key): 
					wrong[key]=[[],[],[],[],[],[]] #start, end, P, Z, D, Z62
				##Fill it merging
				while k < (len(i)-1): #While k is not the index of the last cell
					if not windowsMaster[filename][key][1][i[k]]>= windowsMaster[filename][key][0][i[(k+1)]]:  #If they do not overlap
						wrong[key][0].append(windowsMaster[filename][key][0][i[k]])
						wrong[key][1].append(windowsMaster[filename][key][1][i[k]])
						wrong[key][2].append(windowsMaster[filename][key][2][i[k]])
						wrong[key][3].append(windowsMaster[filename][key][3][i[k]])
						wrong[key][4].append(windowsMaster[filename][key][4][i[k]])
						wrong[key][5].append(windowsMaster[filename][key][5][i[k]])
						k+=1
					else:
						count=1
						P= windowsMaster[filename][key][2][i[k]]
						Z= windowsMaster[filename][key][3][i[k]]
						D= windowsMaster[filename][key][4][i[k]]
						Z62= windowsMaster[filename][key][5][i[k]]
						start=windowsMaster[filename][key][0][i[k]]
						while k<(len(i)-1) and windowsMaster[filename][key][1][i[k]]>= windowsMaster[filename][key][0][i[(k+1)]]-1: #If they overlap. 
							end=windowsMaster[filename][key][1][i[(k+1)]]
							P+= windowsMaster[filename][key][2][i[(k+1)]]
							Z+= windowsMaster[filename][key][3][i[(k+1)]]
							D+= windowsMaster[filename][key][4][i[(k+1)]]
							Z62+= windowsMaster[filename][key][5][i[(k+1)]]
							count+=1
							k+=1
						else: #Once I have the values, fill
							k+=1 #I have to add another 1 so that this row is not repeated as is if was not mergeable even though it has been merged
							wrong[key][0].append(start)
							wrong[key][1].append(end)
							wrong[key][2].append(float(P)/count)
							wrong[key][3].append(float(Z)/count)
							wrong[key][4].append(float(D)/count)
							wrong[key][5].append(float(Z62)/count)
				else: #When k is the index of the last cell
					if k == (len(i)-1): #If it was left alone at last not merged with the prev window it's because it's unmergeable. So just fill
						wrong[key][0].append(windowsMaster[filename][key][0][i[k]])
						wrong[key][1].append(windowsMaster[filename][key][1][i[k]])
						wrong[key][2].append(windowsMaster[filename][key][2][i[k]])
						wrong[key][3].append(windowsMaster[filename][key][3][i[k]])
						wrong[key][4].append(windowsMaster[filename][key][4][i[k]])
						wrong[key][5].append(windowsMaster[filename][key][5][i[k]])
		if len(wrong)>0: #Populate wrongMaster
			wrongMaster[filename]=wrong #start, end, P, Z, D, Z62	
	return wrongMaster
	


########

def classifyWindows(windowsMaster, Pthresh , Zthresh , Dthresh, Z62thresh ):
	for filename in windowsMaster:
		for key in windowsMaster[filename]:
			for i in range(len(windowsMaster[filename][key][2])):
				if windowsMaster[filename][key][2][i]<Pthresh and windowsMaster[filename][key][3][i]>=Zthresh and windowsMaster[filename][key][4][i]<Dthresh and windowsMaster[filename][key][5][i]>=Z62thresh:
					windowsMaster[filename][key][6].append(1)
				else:
					windowsMaster[filename][key][6].append(0)
	return windowsMaster	


##########

def defineThresholds(windowsMaster, Psd, Zsd, Dsd, Z62sd):
	#Initialize lists
	Ps=[]
	Zs=[]
	Ds=[]
	Z62s=[]
	for filename in windowsMaster: #Fill in the lists with all the values from all the windows from all the filenames
		for key in windowsMaster[filename]:
			Ps+=windowsMaster[filename][key][2]
			Zs+=windowsMaster[filename][key][3]
			Ds+=windowsMaster[filename][key][4]
			Z62s+=windowsMaster[filename][key][5]
	Pthresh=array(Ps).mean()-(array(Ps).std()*Psd) #These are my parameters thresholds
	Zthresh=array(Zs).mean()+(array(Zs).std()*Zsd)
	Dthresh=array(Ds).mean()-(array(Ds).std()*Dsd)
	Z62thresh=array(Z62s).mean()+(array(Z62s).std()*Z62sd)
	return (Pthresh, Zthresh, Dthresh, Z62thresh)


##### I had to make a dictionary with all the output with the values to then calculate the thresholds

def fillParametersLists(scoresProb, scoresMatrixPerCol, name, P, windows, start, end):
	if not windows.has_key(name): #initialize
		windows[name]=[[],[],[],[],[],[],[]]
	windows[name][0].append(start)
	windows[name][1].append(end)
	windows[name][2].append(P)
	windows[name][3].append(scoresProb[name][1])
	windows[name][4].append(scoresMatrixPerCol[name][0])
	windows[name][5].append(scoresMatrixPerCol[name][1])
	return windows

##### Functionalize the window trimming step

def trimWindowEdges(windowAlign, align, aaFreqsDict,name, start, end, firstAA,lastAA ):
	while lastAA!="NA" and lastAA>0.9 and end>start: #If lastAA is too conserved, I'll trim it
		windowAlign=align[:, start-1:end-1]
		end=end-1
		firstAA,lastAA=getFirstLastAA(windowAlign,aaFreqsDict,name)
	while firstAA!="NA" and firstAA>0.9 and start<end: #If firstAA is too conserved, I'll trim it
		windowAlign=align[:, start:end]
		start=start+1
		firstAA,lastAA=getFirstLastAA(windowAlign,aaFreqsDict,name)
	return (windowAlign, start, end)


def getFirstLastAA(windowAlign,aaFreqsDict,name):
    currentWindowSpeciesProbs={} # [species]=[probsite1, probsite2...]
    currentWindowSpeciesProbs[name]=[]
    window=len(windowAlign[1].seq)
    i=0
    while i<window:
        col= windowAlign[:, i]
        cleancol=col.replace("-", "").replace("X", "")
        # populate the amino acids frequencies per site per column 
        for key in aaFreqsDict:
            if key!="-" and key!="X" and len(cleancol)>0 :
                aaFreqsDict[key]=round(cleancol.count(key)*1.00/len(cleancol), 10)
            else:
                aaFreqsDict[key]="NA"
        # save the prob of a species site aa
        for record in windowAlign:
	    if record.id==name:
                currentAA=record.seq[i]
                siteProb=aaFreqsDict[currentAA]
            	currentWindowSpeciesProbs[name].append(siteProb)
        i+=1
    firstAA=currentWindowSpeciesProbs[name][0]
    lastAA=currentWindowSpeciesProbs[name][-1]
    return(firstAA,lastAA)


def getAveProbPerWindowPerSeq(align, pGaps,aaFreqsDict):
    scorePerSpeciesPerWindow={} # [species]=[Z-score w1, Z-score w2...]
    currentWindowSpeciesProbs={} # [species]=[probsite1, probsite2...]
    aveWindowSpeciesProbs={} # [species]=windowprob
    window=len(align[1].seq)
    i=0
    while i<window:
        col=align[:, i]
        cleancol=col.replace("-", "").replace("X", "")
        # populate the amino acids frequencies per site per column 
        for key in aaFreqsDict:
            if key!="-" and key!="X" and len(cleancol)>0 :
                aaFreqsDict[key]=round(cleancol.count(key)*1.00/len(cleancol), 10)
            else:
                aaFreqsDict[key]="NA"
        # save the prob of a species site aa
        for name in align:
            species=name.id
            currentAA=name.seq[i]
            siteProb=aaFreqsDict[currentAA]
            if not currentWindowSpeciesProbs.has_key(species):
                currentWindowSpeciesProbs[species]=[]        
            currentWindowSpeciesProbs[species].append(siteProb)
        i+=1
    #calculate ave prob per window  per seq
    allAveWindow=[]
    for species in currentWindowSpeciesProbs:
        probs=currentWindowSpeciesProbs[species]
        probsNumeric=[x for x in probs if x!="NA" ]
        if len(probsNumeric)>pGaps*window: #if not mostly-gap seq
            aveProb=array(probsNumeric).mean()
            aveWindowSpeciesProbs[species]=round(aveProb, 2)
            allAveWindow.append(aveProb)
        else:
            aveWindowSpeciesProbs[species]="NA"
    # calculate Z-score
    if len(allAveWindow)>0: #not a gap-only window ### LisNote: Before it did not have the len(), now I added it here. This is the changed good program to run for the tests
        aveWindow=array(allAveWindow).mean()
        stdevWindow=array(allAveWindow).std()
    for species in aveWindowSpeciesProbs:
        probSpecies=aveWindowSpeciesProbs[species]
        Z=0
        if not probSpecies=="NA": #not a gap-only window
            if not stdevWindow<0.001: #very conserved window
                dev=abs(aveWindow-probSpecies)
                Z=round(dev/stdevWindow, 2)
            else:
                Z=0
        else:
            Z="NA"
        if not scorePerSpeciesPerWindow.has_key(species):
            scorePerSpeciesPerWindow[species]=[]
        score=[probSpecies,Z]
        scorePerSpeciesPerWindow[species]=score
    return scorePerSpeciesPerWindow


def readMatrix(nameMatrix):
    #takes an half similarity matrix and turns into a dict
    input=open(nameMatrix)
    lines=input.readlines()
    names=lines[0].rstrip().split()
    aaDict={} #aaDict[(aa1,aa2)]=value
    i=1
    while i<len(lines):
        info=lines[i].rstrip().split()
        aa1=info[0]
        j=1
        while j<len(info):
            aa2=names[j-1]
            value=int(info[j])
            pair=(aa1,aa2)
            aaDict[pair]=value
            pairI=(aa2,aa1)
            aaDict[pairI]=value
            j+=1
        i+=1
    return aaDict


def scorePairAlignedSeqsWithMatrix(seq1,seq2,matrixAA):
    #takes 2 aa seqs as strings with the same size, 
    #and calculates the score of aligned aas ignoring gaps
    totalScore=[]
    gaps=0
    i=0
    while i<len(seq1):
        if seq1[i]!="-" and seq2[i]!="-" and seq1[i]!="X" and seq2[i]!="X":
            totalScore.append(matrixAA[(seq1[i],seq2[i])])
        else:
            gaps+=1
        i+=1
    if gaps==len(seq1):
        total="NA"
    else:
        total=array(totalScore).mean()
    return total


def getPairwiseScore(windowAlign, matrixAA):    
    scoreSlice={} # [species]=lowestScore
    allSeqs={}
    forAverage=[]
    for entry in windowAlign:
        name=entry.id
        allSeqs[name]=str(entry.seq)
    for name in allSeqs:
        allScores=[]
        others=allSeqs.keys()
        others.remove(name)
        for alt in others:
            seq1=allSeqs[name]
            seq2=allSeqs[alt]
            score=scorePairAlignedSeqsWithMatrix(seq1, seq2, matrixAA)            
            if score!="NA":
                allScores.append(score)            
        if len(allScores)>0:
            score=round(allScores[-1], 2)
            forAverage.append(score)
            scoreSlice[name]=[score]        
    aveWindow=array(forAverage).mean()
    stdevWindow=array(forAverage).std()        
    for species in allSeqs:
        if scoreSlice.has_key(species):
            score=round(scoreSlice[species][0], 2)            
            if not stdevWindow<0.001: #very conserved window
                dev=abs(aveWindow-score)
                Z=round(dev/stdevWindow, 2)                
            else:
                Z=0                
            scoreSlice[species].append(Z)
        else:
            scoreSlice[species]=["NA", "NA"]        
    return scoreSlice


def getColumnDistScoreToClosest(windowAlign, matrixAA):    
    scoreSlice={} # [species]=lowestScore
    allSeqs={}
    forAverage=[]
    allScores=[]    
    for entry in windowAlign:
        name=entry.id
        allSeqs[name]=str(entry.seq)
    for name in allSeqs:
        seq1=allSeqs[name]
        others=allSeqs.keys()
        others.remove(name)
        windowScore=[]        
        i=0
        while i<len(seq1):
            scoreCol=[]
            #score each position of the alignment 
            for alt in others:            
                seq2=allSeqs[alt][i]
                score=scorePairAlignedSeqsWithMatrix(seq1[i], seq2, matrixAA)
                if score!="NA":
                    scoreCol.append(score)            
            if len(scoreCol)>0: 
                scoreCol.sort()
                windowScore.append(scoreCol[-1]) #get the highest score for this column                  
            i+=1        
        if len(windowScore)>0:
            scoreAVEwindow=array(windowScore).mean()
            scoreAVEwindow=round(scoreAVEwindow, 2)            
            allScores.append(scoreAVEwindow)
            scoreSlice[name]=[scoreAVEwindow]       
    if len(allScores)>0:
    	aveWindow=array(allScores).mean()
    	stdevWindow=array(allScores).std()    
    for species in allSeqs:
        if scoreSlice.has_key(species):
            score=scoreSlice[species][0]            
            if not stdevWindow<0.001: #very conserved window
                dev=abs(aveWindow-score)
                Z=round(dev/stdevWindow, 2)
            else:
                Z=0                
            scoreSlice[species].append(Z)
        else:
            scoreSlice[species]=["NA", "NA"]        
    return scoreSlice


