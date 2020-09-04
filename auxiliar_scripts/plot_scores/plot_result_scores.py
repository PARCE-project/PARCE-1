#!/usr/bin/python3

"""
PARCE: Protocol for Amino acid Refinement through Computational Evolution
Script to plot the average scores from the design protocol

From publication "PARCE: Protocol for Amino acid Refinement through Computational Evolution"
Computer Physics Communications 
Authors: Rodrigo Ochoa, Miguel A. Soler, Alessandro Laio, Pilar Cossio
Year: 2020

Explanation:

This script generates plots showing the evolution of the scores, and the mutations accepted and rejected during
the design. The input for that analysis is the mutation_report.txt file obtained after running a design cycle.
To run the script first update it with the path where the mutation_report.txt file is located, and run the
command 'python3 plot_result_scores.py'. Two plots are created for each scoring function: one has only the accepted
mutations whereas the other has also the rejected mutations (red points).
"""

# Import local modules
import matplotlib.pyplot as plt
from statistics import mean
from statistics import stdev
from statistics import variance

# VARIABLES TO MODIFY
path="/home/PARCE-1/auxiliar_scripts/plot_scores"
score_list=["bach","pisa","zrank","irad","bmf-bluues","firedock"]

# Score lists
bach=[]
pisa=[]
zrank=[]
irad=[]
totalbmfbluues=[]
firedock=[]

# Score lists of accepted mutations
bachAccepted=[]
pisaAccepted=[]
zrankAccepted=[]
iradAccepted=[]
totalbmfbluuesAccepted=[]
firedockAccepted=[]

# Additional fields
xvalues=[]
mutations=[]
acceptance=[]
sequence=[]
markers=[]

counter_mutations=0
mutation_data=[x.strip() for x in open(path+"/mutation_report.txt")]

for line in mutation_data:
    fields=line.split()
    if "Iteration" in fields[0]:
        counter_mutations+=1
        xvalues.append(counter_mutations)
        mutations.append(fields[1])
        acceptance.append(fields[3])
        sequence.append(fields[len(fields)-1].split(":")[1])
        
        for f in fields:
            if "bach" in f:
                if "bach" in score_list:
                    bach.append(float(f.split(":")[1]))
                    if fields[3]=="Accepted":
                        bachAccepted.append(float(f.split(":")[1]))
                        markers.append(9)
                    else:
                        bachAccepted.append(bachAccepted[-1])
                        markers.append(6)
            if "pisa" in f:
                if "pisa" in score_list:
                    pisa.append(float(f.split(":")[1]))
                    if fields[3]=="Accepted":
                        pisaAccepted.append(float(f.split(":")[1]))
                    else:
                        pisaAccepted.append(pisaAccepted[-1])
            if "zrank" in f:
                if "zrank" in score_list:
                    zrank.append(float(f.split(":")[1]))
                    if fields[3]=="Accepted":
                        zrankAccepted.append(float(f.split(":")[1]))
                    else:
                        zrankAccepted.append(zrankAccepted[-1])
            if "irad" in f:
                if "irad" in score_list:
                    irad.append(float(f.split(":")[1]))
                    if fields[3]=="Accepted":
                        iradAccepted.append(float(f.split(":")[1]))
                    else:
                        iradAccepted.append(iradAccepted[-1])
            if "totalbmfbluues" in f:
                if "bmf-bluues" in score_list:
                    totalbmfbluues.append(float(f.split(":")[1]))
                    if fields[3]=="Accepted":
                        totalbmfbluuesAccepted.append(float(f.split(":")[1]))
                    else:
                        totalbmfbluuesAccepted.append(totalbmfbluuesAccepted[-1])
            if "firedock" in f:
                if "firedock" in score_list:
                    firedock.append(float(f.split(":")[1]))
                    if fields[3]=="Accepted":
                        firedockAccepted.append(float(f.split(":")[1]))
                    else:
                        firedockAccepted.append(firedockAccepted[-1])

# Generation of plots
if "bach" in score_list:
    # Plot with all the scores
    plt.plot(xvalues,bach,color='b',linewidth=2.0,alpha=1.0,label='BACH')
    for i,m in enumerate(markers):
        if m==9:
            plt.plot(xvalues[i],bach[i],color='b',linestyle='None',marker="o",markersize=m,alpha=1.0)
        if m==6:
            plt.plot(xvalues[i],bach[i],color='r',linestyle='None',marker="o",markersize=m,alpha=1.0)
    
    axes = plt.gca()
    axes.set_xlim([0,xvalues[-1]])
    plt.title('BACH score per mutation')
    plt.xlabel('Attempted mutations',fontsize=18)
    plt.ylabel('Score',fontsize=18)
    plt.legend(loc='upper right',fontsize='large',numpoints=1)
    plt.savefig("bach_score.png")
    plt.close()
    
    # Plot with only the accepted scores
    plt.plot(xvalues,bachAccepted,color='b',linewidth=2.0,alpha=1.0,label='BACH')
    for i,m in enumerate(markers):
        if m==9:
            plt.plot(xvalues[i],bach[i],color='b',linestyle='None',marker="o",markersize=m,alpha=1.0)
            
    axes = plt.gca()
    axes.set_xlim([0,xvalues[-1]])
    plt.title('BACH score per mutation')
    plt.xlabel('Mutations',fontsize=18)
    plt.ylabel('Score',fontsize=18)
    plt.legend(loc='upper right',fontsize='large',numpoints=1)
    plt.savefig("bach_accepted_score.png")
    plt.close()

if "zrank" in score_list:
    # Plot with all the scores
    plt.plot(xvalues,zrank,color='g',linewidth=2.0,alpha=1.0,label='ZRANK')
    for i,m in enumerate(markers):
        if m==9:
            plt.plot(xvalues[i],zrank[i],color='g',linestyle='None',marker="o",markersize=m,alpha=1.0)
        if m==6:
            plt.plot(xvalues[i],zrank[i],color='r',linestyle='None',marker="o",markersize=m,alpha=1.0)
    
    axes = plt.gca()
    axes.set_xlim([0,xvalues[-1]])
    plt.title('ZRANK score per mutation')
    plt.xlabel('Mutations',fontsize=18)
    plt.ylabel('Score',fontsize=18)
    plt.legend(loc='upper right',fontsize='large',numpoints=1)
    plt.savefig("zrank_score.png")
    plt.close()
    
    plt.plot(xvalues,zrankAccepted,color='g',linewidth=2.0,alpha=1.0,label='ZRANK')
    for i,m in enumerate(markers):
        if m==9:
            plt.plot(xvalues[i],zrank[i],color='g',linestyle='None',marker="o",markersize=m,alpha=1.0)
            
    axes = plt.gca()
    axes.set_xlim([0,xvalues[-1]])
    plt.title('ZRANK score per mutation')
    plt.xlabel('Mutations',fontsize=18)
    plt.ylabel('Score',fontsize=18)
    plt.legend(loc='upper right',fontsize='large',numpoints=1)
    plt.savefig("zrank_accepted_score.png")
    plt.close()

if "irad" in score_list:
    # Plot with all the scores
    plt.plot(xvalues,irad,color='m',linewidth=2.0,alpha=1.0,label='IRAD')
    for i,m in enumerate(markers):
        if m==9:
            plt.plot(xvalues[i],irad[i],color='m',linestyle='None',marker="o",markersize=m,alpha=1.0)
        if m==6:
            plt.plot(xvalues[i],irad[i],color='r',linestyle='None',marker="o",markersize=m,alpha=1.0)
    
    axes = plt.gca()
    axes.set_xlim([0,xvalues[-1]])
    plt.title('IRAD score per mutation')
    plt.xlabel('Mutations',fontsize=18)
    plt.ylabel('Score',fontsize=18)
    plt.legend(loc='upper right',fontsize='large',numpoints=1)
    plt.savefig("irad_score.png")
    plt.close()
    
    plt.plot(xvalues,iradAccepted,color='m',linewidth=2.0,alpha=1.0,label='IRAD')
    for i,m in enumerate(markers):
        if m==9:
            plt.plot(xvalues[i],irad[i],color='m',linestyle='None',marker="o",markersize=m,alpha=1.0)
            
    axes = plt.gca()
    axes.set_xlim([0,xvalues[-1]])
    plt.title('IRAD score per mutation')
    plt.xlabel('Mutations',fontsize=18)
    plt.ylabel('Score',fontsize=18)
    plt.legend(loc='upper right',fontsize='large',numpoints=1)
    plt.savefig("irad_accepted_score.png")
    plt.close()

if "firedock" in score_list:
    # Plot with all the scores
    plt.plot(xvalues,firedock,color='gold',linewidth=2.0,alpha=1.0,label='Firedock')
    for i,m in enumerate(markers):
        if m==9:
            plt.plot(xvalues[i],firedock[i],color='gold',linestyle='None',marker="o",markersize=m,alpha=1.0)
        if m==6:
            plt.plot(xvalues[i],firedock[i],color='r',linestyle='None',marker="o",markersize=m,alpha=1.0)
    
    axes = plt.gca()
    axes.set_xlim([0,xvalues[-1]])
    plt.title('Firedock score per mutation')
    plt.xlabel('Mutations',fontsize=18)
    plt.ylabel('Score',fontsize=18)
    plt.legend(loc='upper right',fontsize='large',numpoints=1)
    plt.savefig("firedock_score.png")
    plt.close()
    
    plt.plot(xvalues,firedockAccepted,color='gold',linewidth=2.0,alpha=1.0,label='Firedock')
    for i,m in enumerate(markers):
        if m==9:
            plt.plot(xvalues[i],firedock[i],color='gold',linestyle='None',marker="o",markersize=m,alpha=1.0)
            
    axes = plt.gca()
    axes.set_xlim([0,xvalues[-1]])
    plt.title('Firedock score per mutation')
    plt.xlabel('Mutations',fontsize=18)
    plt.ylabel('Score',fontsize=18)
    plt.legend(loc='upper right',fontsize='large',numpoints=1)
    plt.savefig("firedock_accepted_score.png")
    plt.close()

if "pisa" in score_list:
    # Plot with all the scores
    plt.plot(xvalues,pisa,color='pink',linewidth=2.0,alpha=1.0,label='Pisa')
    for i,m in enumerate(markers):
        if m==9:
            plt.plot(xvalues[i],pisa[i],color='pink',linestyle='None',marker="o",markersize=m,alpha=1.0)
        if m==6:
            plt.plot(xvalues[i],pisa[i],color='r',linestyle='None',marker="o",markersize=m,alpha=1.0)
    
    axes = plt.gca()
    axes.set_xlim([0,xvalues[-1]])
    plt.title('Pisa score per mutation')
    plt.xlabel('Mutations',fontsize=18)
    plt.ylabel('Score',fontsize=18)
    plt.legend(loc='upper right',fontsize='large',numpoints=1)
    plt.savefig("pisa_score.png")
    plt.close()
    
    plt.plot(xvalues,pisaAccepted,color='pink',linewidth=2.0,alpha=1.0,label='Pisa')
    for i,m in enumerate(markers):
        if m==9:
            plt.plot(xvalues[i],pisa[i],color='pink',linestyle='None',marker="o",markersize=m,alpha=1.0)
            
    axes = plt.gca()
    axes.set_xlim([0,xvalues[-1]])
    plt.title('Pisa score per mutation')
    plt.xlabel('Mutations',fontsize=18)
    plt.ylabel('Score',fontsize=18)
    plt.legend(loc='upper right',fontsize='large',numpoints=1)
    plt.savefig("pisa_accepted_score.png")
    plt.close()
    
if "bmf-bluues" in score_list:
    # Plot with all the scores
    plt.plot(xvalues,totalbmfbluues,color='grey',linewidth=2.0,alpha=1.0,label='BMF-BLUUES')
    for i,m in enumerate(markers):
        if m==9:
            plt.plot(xvalues[i],totalbmfbluues[i],color='grey',linestyle='None',marker="o",markersize=m,alpha=1.0)
        if m==6:
            plt.plot(xvalues[i],totalbmfbluues[i],color='r',linestyle='None',marker="o",markersize=m,alpha=1.0)
    
    axes = plt.gca()
    axes.set_xlim([0,xvalues[-1]])
    plt.title('BMF-BLUUES score per mutation')
    plt.xlabel('Mutations',fontsize=18)
    plt.ylabel('Score',fontsize=18)
    plt.legend(loc='upper right',fontsize='large',numpoints=1)
    plt.savefig("bmf-bluues_score.png")
    plt.close()
    
    
    plt.plot(xvalues,totalbmfbluuesAccepted,color='grey',linewidth=2.0,alpha=1.0,label='BMF-BLUUES')
    for i,m in enumerate(markers):
        if m==9:
            plt.plot(xvalues[i],totalbmfbluues[i],color='grey',linestyle='None',marker="o",markersize=m,alpha=1.0)
            
    axes = plt.gca()
    axes.set_xlim([0,xvalues[-1]])
    plt.title('BMF-BLUUES score per mutation')
    plt.xlabel('Mutations',fontsize=18)
    plt.ylabel('Score',fontsize=18)
    plt.legend(loc='upper right',fontsize='large',numpoints=1)
    plt.savefig("bmf-bluues_accepted_score.png")
    plt.close()
