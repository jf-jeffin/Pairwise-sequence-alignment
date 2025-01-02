'''
Pair-Wise Sequence Alignment
Local alignment: Smith Waterman Algorithm
Global Alignment: Needleman-Wunsch Algorithm
'''
def sco(seq1: str,seq2: str) -> int:
    #return score
    if seq1==seq2:
        return 3
    else:
        return -1

def tbstartloc(scoreMatrix: list[list],local_ali,maxvalues=1):
    #returns location in matrix to start trace back    
    if local_ali=='l' and maxvalues:#location with max score
        maxs=0
        for r in range(1,ls+1):#r for row, nr for next row
            for c in range(1,gs+1):#c for column, nc for next col
                if scoreMatrix[r][c]>maxs:
                    maxs=scoreMatrix[r][c]
                    mr=r
                    mc=c
        return [mr,mc]
    elif local_ali !='l':#last col row entry
         r=len(scoreMatrix)
         c=len(scoreMatrix[r-1])
         return [r-1,c-1]

def initi(gs: int,ls: int, local_ali)->list[list]:
    #return score matrix
    gap=-2
    if local_ali=='l':
        scorMatr=[[0]]
        scorMatr=scorMatr*(ls+1)
        for i in range(len(scorMatr)):  
            scorMatr[i]=scorMatr[i]*(gs+1)
        return scorMatr
    else:#global alignment
        scorMatr=[[x*gap] for x in range(ls)] 
        scorMatr[0]=[x*gap for x in range(gs)]
        return scorMatr

#load two sequences 
s=input("Enter first sequence: ")
sd=input("Enter second sequence: ")
local_ali=input("Enter l for local alignment or g for global aligment: ")
     
#initialize matrix      
if len(s)<len(sd):#row will be sequence with larger length
    ls=len(s)
    gs=len(sd)  
else:
    ls=len(sd)
    gs=len(s)    

sm=initi(gs,ls,local_ali)
    
#populating matrix with score
gap=-2
r=c=1
sc=sr=0#for sequence index
if local_ali=='l':
    while r<len(sm) and sr<ls:
        while c<len(sm[0]) and sc<gs:
            pr,pc=r-1,c-1
            sm[r][c]=max(0,sm[pr][pc]+sco(s[sc],
                        sd[sr]),sm[pr][c]+gap,
                         sm[r][pc]+gap)
            c,sc=c+1,sc+1
        c,sc=1,0
        r,sr=r+1,sr+1
else:
    while r<len(sm) and sr<ls:
        while c<len(sm[0]) and sc<gs :
            pr,pc=r-1,c-1#previous row(pr) 
            sm[r].append(max(sm[pr][pc]+sco(s[sc],sd[sr]),sm[pr][c]+gap,sm[r][pc]+gap))
            c,sc=c+1,sc+1
        c,sc=1,0
        r,sr=r+1,sr+1
            
#traceback
startl=tbstartloc(sm,local_ali)
mr,mc=startl

cu_sc=0  #cutoff score  
ase=""#alignment sequences
ase1=""
while sm[mr][mc]:

    #score matrix is initialized with first col and first row with zeros
    #so corresponding index in seq will be scorematrix index-1
    nr,nc=mr-1,mc-1  
    if sd[nr]==s[nc]:
        ase=ase+s[nc]
        ase1=ase1+sd[nr]
        mr,mc=mr-1,mc-1
    elif sd[nr]!=s[nc]:#if sequence element is not equal
        ne=[sm[nr][nc]
            ,sm[nr][mc],
            sm[mr][nc]]
        nseq=max(ne)#max score among neighbouring upper triangle entries 
        if nseq==ne[0]:
            ase=ase+s[nc]
            ase1=ase1+sd[nr]
            mr,mc=mr-1,mc-1
        elif  nseq==ne[1] :              
            ase1= ase1+sd[nr]
            ase=ase+"_"
            mr=nr
        elif nseq==ne[2]:               
            ase1=ase1+"_"
            ase=ase+s[nc]          
            mc=nc 
aliseq,aliseq2=ase[::-1],ase1[::-1]

#alignment score
match=unmatch=gap=0
for s1,s2 in zip(aliseq,aliseq2):
    if s1==s2:
        match=match+1
    elif s1=="_" or s2 =="_":
        gap=gap+1
    else:
        unmatch=unmatch+1

print(aliseq)
print(aliseq2)
print(f"{match=},{unmatch=},{gap=}")
