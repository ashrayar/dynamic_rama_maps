# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 16:42:23 2021

@author: ashraya
"""
import matplotlib.pyplot as plt
plt.ioff()
#md_phi=[]
#md_psi=[]
#infile=open("/home/ashraya/RMap_Work/rmap_dynamics/2vb1_opls/2vb1_GLU35_phipsi.csv","r")
#for line in infile:
#    phi,psi=line[0:-1].split(",")
#    md_phi.append(float(phi))
#    md_psi.append(float(psi))
#infile.close()

count=0
for i in range(20001,21001):
    print i
    infile=open("/home/ashraya/RMap_Work/rmap_dynamics/2vb1_opls/frame_maps/2vb1_GLU35_frame_"+str(i)+".csv","r")
    full_all_phi=[]
    full_all_psi=[]
    part_all_phi=[]
    part_all_psi=[]
    for line in infile:
        lineparts=line.split(",")
        phi=int(lineparts[0])
        psi=int(lineparts[1])
        apd=int(lineparts[2][0:-1])
        if apd==1:
            full_all_phi.append(phi)
            full_all_psi.append(psi)
        elif apd==2:
            part_all_phi.append(phi)
            part_all_psi.append(psi)
    infile.close()

    plt.figure(figsize=(7,7))
    if len(full_all_phi)!=0:
        plt.scatter(full_all_phi,full_all_psi,s=5,alpha=1,c="forestgreen",color="forestgreen",zorder=1)
    if len(part_all_phi)!=0:
        plt.scatter(part_all_phi,part_all_psi,s=5,alpha=1,c="burlywood",color="burlywood",zorder=1)
#    plt.scatter(md_phi[count],md_psi[count],s=50,c="blue",zorder=2)
    plt.hlines(0,-180,180)
    plt.vlines(0,-180,180)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel('$\phi$ (in $\degree$)',fontsize=22)
    plt.ylabel('$\psi$ (in $\degree$)',fontsize=22)
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.savefig("/home/ashraya/RMap_Work/rmap_dynamics/2vb1_opls/frame_maps/2vb1_GLU35_frame_"+str(i)+"_nophisi_map.png",dpi=300,bbox_inches="tight")
    count+=1




