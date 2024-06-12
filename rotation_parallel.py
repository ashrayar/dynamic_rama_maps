# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 10:42:36 2018

@author: ashraya
"""

import os
import math
import multiprocessing as mp
 
import numpy as np
      
def cleanPdbFile(atom_lines_array,missing_res=[]):
    cur_res=None
    residue=dict()
    towrite=[]
    chain_id=" "
    for line in atom_lines_array:
        if line[0:4]=="ATOM" and line[21]==chain_id:
#            if line[26]!=" ":
#                continue
            if cur_res==None:
                cur_res=int(line[22:26])
            elif cur_res!=int(line[22:26]) and cur_res not in missing_res: #new residue and not a missing residue, write N, CA and C atom lines, in that order (order important for phi-psi calculation)
                if residue.has_key("N"):
                    towrite.append(residue['N'])
                if residue.has_key("H"):
                    towrite.append(residue["H"])
                if residue.has_key("CD"):
                    towrite.append(residue["CD"])
                if residue.has_key("CA"):
                    towrite.append(residue['CA'])
                if residue.has_key("C"):
                    towrite.append(residue['C'])
                for item in residue.keys():
                    if item not in ['N','CA','C','H','CD']: #write rest of the atom lines of the residue
                        towrite.append(residue[item])
                residue=dict()
                cur_res=int(line[22:26])
                
            if cur_res==int(line[22:26]):
                if line[26]!=" ":
                    residue=dict()
                    cur_res=int(line[22:26])
                else:
                    if not residue.has_key(line[12:16].strip()):
                        if line[12:16].strip() == "HN":
                            residue["H"]=line
                        else:
                            residue[line[12:16].strip()]=line #store atom lines of non-missing residues
        elif line[0:6]=="ANISOU": #ignore anisou lines
            continue
        else:
            towrite.append(line)
    if residue.has_key("N"):
        towrite.append(residue['N'])
    if residue.has_key("H"):
        towrite.append(residue["H"])
    if residue.has_key("CD"):
        towrite.append(residue["CD"])
    if residue.has_key("CA"):
        towrite.append(residue['CA'])
    if residue.has_key("C"):
        towrite.append(residue['C'])
    for item in residue.keys():
        if item not in ['N','CA','C','H','CD']: #write rest of the atom lines of the residue
            towrite.append(residue[item])
    return towrite

def findMissingRes(pdb_id,chain_id):            
    infile=open(pdb_id+".pdb","r")
    cur_res=None
    missing_res=[]
    atoms=dict()
    for line in infile:
        if line[0:4] == "ATOM" and line[21]==chain_id:
            if cur_res==None:
                cur_res=int(line[22:26])
            if cur_res!=int(line[22:26]): #new residue
                if cur_res+1!=int(line[22:26]): #residue number currently read is not +1 of previously read residue
                    added_atoms=atoms.keys()
                    if "N" not in added_atoms or "C" not in added_atoms or "CA" not in added_atoms: #if any of the three backbone atoms not present in a residue, consider it as missing
                        missing_res.append(cur_res)
                    for i in range(cur_res+1,int(line[22:26])):
                        missing_res.append(i) 
                else:
                    added_atoms=atoms.keys()
                    if "N" not in added_atoms or "C" not in added_atoms or "CA" not in added_atoms: #if any of the three backbone atoms not present in a residue, consider it as missing
                        missing_res.append(cur_res)
                cur_res=int(line[22:26])
                atoms=dict()
            if line[21]==chain_id and line[12:16].strip() in ["N","CA","C"] and line[26]==" ":
                atom_type=line[12:16].strip()
                x=float(line[30:38])
                y=float(line[38:46])
                z=float(line[46:54])
                if not atoms.has_key(atom_type): #store first conformer (in case of multiple occupancy)
                    atoms[atom_type]=[x,y,z]
    infile.close()
    print missing_res
    return missing_res

def FindDirCosines(A,B):
    delx=B[0]-A[0]
    dely=B[1]-A[1]
    delz=B[2]-A[2]
    denom=math.sqrt(delx*delx+dely*dely+delz*delz)
    l=delx/denom
    m=dely/denom
    n=delz/denom
    a=[l,m,n]
    return a

def FindDihedralAngle(A,B,C,D):

    BCi=C[0]-B[0]
    BCj=C[1]-B[1]
    BCk=C[2]-B[2]
    
    BAi=A[0]-B[0]
    BAj=A[1]-B[1]
    BAk=A[2]-B[2]
    
    CDi=D[0]-C[0]
    CDj=D[1]-C[1]
    CDk=D[2]-C[2]
    
    Q1i=(BCj*BAk)-(BCk*BAj)
    Q1j=(BCk*BAi)-(BCi*BAk)
    Q1k=(BCi*BAj)-(BCj*BAi)
    
    Q2i=(BCj*CDk)-(BCk*CDj)
    Q2j=(BCk*CDi)-(BCi*CDk)
    Q2k=(BCi*CDj)-(BCj*CDi)
    magQ1=math.sqrt((Q1i*Q1i)+(Q1j*Q1j)+(Q1k*Q1k))
    Q1i=Q1i/magQ1
    Q1j=Q1j/magQ1
    Q1k=Q1k/magQ1
    
    magQ2=math.sqrt((Q2i*Q2i)+(Q2j*Q2j)+(Q2k*Q2k))
    Q2i=Q2i/magQ2
    Q2j=Q2j/magQ2
    Q2k=Q2k/magQ2
    
    Q1dotQ2=(Q1i*Q2i)+(Q1j*Q2j)+(Q1k*Q2k)
    Q1dotQ2=round(Q1dotQ2,10)
    chi=math.acos(Q1dotQ2)
    chinew=math.degrees(chi)
    
    Q1=np.array([Q1i,Q1j,Q1k])
    Q2=np.array([Q2i,Q2j,Q2k])
    Q1crossQ2=np.cross(Q1,Q2)
    magBC=math.sqrt((BCi*BCi)+(BCj*BCj)+(BCk*BCk))
    unitBCi=BCi/magBC
    unitBCj=BCj/magBC
    unitBCk=BCk/magBC
    unitBC=np.array([unitBCi,unitBCj,unitBCk])
    anglesign=np.dot(Q1crossQ2,unitBC)
    if anglesign<0:
        chinew=chinew*-1
        
    return int(round(chinew))

def resultantpoint(l,m,n,a,b,c,x,y,z,angle):
    angle=math.radians(angle)
    term1=(a*(m*m+n*n)-l*(b*m+c*n-l*x-m*y-n*z))*(1-math.cos(angle))
    term2=x*math.cos(angle)
    term3=(b*n-c*m-n*y+m*z)*math.sin(angle)
    x_new=term1+term2+term3
    term1=(b*(l*l+n*n)-m*(a*l+c*n-l*x-m*y-n*z))*(1-math.cos(angle))
    term2=y*math.cos(angle)
    term3=(c*l-a*n+n*x-l*z)*math.sin(angle)
    y_new=term1+term2+term3
    term1=(c*(l*l+m*m)-n*(a*l+b*m-l*x-m*y-n*z))*(1-math.cos(angle))
    term2=z*math.cos(angle)
    term3=(a*m-b*l-m*x+l*y)*math.sin(angle)
    z_new=term1+term2+term3
    newcoords=[0 for z in range(3)]
    newcoords[0]=x_new
    newcoords[1]=y_new
    newcoords[2]=z_new
    return newcoords

def FindDistance(a1,a2):
    return round(math.sqrt(((a2[0]-a1[0])*(a2[0]-a1[0]))+((a2[1]-a1[1])*(a2[1]-a1[1]))+((a2[2]-a1[2])*(a2[2]-a1[2]))),3)

def CheckShortContacts(a1,a2,d1,d2):
    dist=FindDistance(a1,a2)
    if dist<=d1 and dist>=d2:
        return "partially allowed: "+str(dist)
    elif dist>d1:
        return "fully allowed: "+str(dist)
    elif dist<d2:
        return "disallowed: "+str(dist)
    

def BringToZero(tempatoms,prepro,phi,psi):
#    phi=FindDihedralAngle(tempatoms['C1'],tempatoms['N2'],tempatoms['CA2'],tempatoms['C2'])
#    psi=FindDihedralAngle(tempatoms['N2'],tempatoms['CA2'],tempatoms['C2'],tempatoms['N3'])
#    print str(phi)+","+str(psi)
    #bring phi to 0
    dircos=FindDirCosines(tempatoms['N2'], tempatoms['CA2'])
    tempatoms['CB']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CB'][0],tempatoms['CB'][1],tempatoms['CB'][2],-phi)
    tempatoms['CA2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],-phi)
    tempatoms['HA']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['HA'][0],tempatoms['HA'][1],tempatoms['HA'][2],-phi)
    tempatoms['C2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['C2'][0],tempatoms['C2'][1],tempatoms['C2'][2],-phi)
    tempatoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['O2'][0],tempatoms['O2'][1],tempatoms['O2'][2],-phi)
    tempatoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['N3'][0],tempatoms['N3'][1],tempatoms['N3'][2],-phi)
    if prepro=="normal":
        tempatoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['H3'][0],tempatoms['H3'][1],tempatoms['H3'][2],-phi)
    else:
        tempatoms['CD3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CD3'][0],tempatoms['CD3'][1],tempatoms['CD3'][2],-phi)
    tempatoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CA3'][0],tempatoms['CA3'][1],tempatoms['CA3'][2],-phi)
    #bring psi to 0
    dircos=FindDirCosines(tempatoms['CA2'],tempatoms['C2'])
    tempatoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['O2'][0],tempatoms['O2'][1],tempatoms['O2'][2],-psi)
    tempatoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['N3'][0],tempatoms['N3'][1],tempatoms['N3'][2],-psi)
    if prepro=="normal":
        tempatoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['H3'][0],tempatoms['H3'][1],tempatoms['H3'][2],-psi)
    else:
        tempatoms['CD3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['CD3'][0],tempatoms['CD3'][1],tempatoms['CD3'][2],-psi)
    tempatoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['CA3'][0],tempatoms['CA3'][1],tempatoms['CA3'][2],-psi)
#    phi=FindDihedralAngle(tempatoms['C1'],tempatoms['N2'],tempatoms['CA2'],tempatoms['C2'])
#    psi=FindDihedralAngle(tempatoms['N2'],tempatoms['CA2'],tempatoms['C2'],tempatoms['N3'])
#    print phi
#    print psi
def add_missing_phipsi(phipsi):
    keys=phipsi.keys()
    if "-178,-180" not in keys and "-178,180" in keys:
        for i in range(-179,180):
            phipsi[str(i)+","+"-180"]=phipsi[str(i)+","+"180"]
    elif "-178,180" not in keys and "-178,-180" in keys:
        for i in range(-179,180):
            phipsi[str(i)+","+"180"]=phipsi[str(i)+","+"-180"]
    if "-180,-178" not in keys and "180,-178" in keys:
        for j in range(-179,180):
            phipsi["-180,"+str(j)]=phipsi["180,"+str(j)]
    elif "180,-178" not in keys and "-180,-178" in keys:
        for j in range(-179,180):
            phipsi["180,"+str(j)]=phipsi["-180,"+str(j)]

    if "-180,180" in keys:
        phipsi["180,180"]=phipsi["-180,180"]
        phipsi["180,-180"]=phipsi["-180,180"]
        phipsi["-180,-180"]=phipsi["-180,180"]
    elif "-180,-180" in keys:
        phipsi["180,180"]=phipsi["-180,-180"]
        phipsi["180,-180"]=phipsi["-180,-180"]
        phipsi["-180,180"]=phipsi["-180,-180"]
    elif "180,-180" in keys:
        phipsi["180,180"]=phipsi["180,-180"]
        phipsi["-180,180"]=phipsi["180,-180"]
        phipsi["-180,-180"]=phipsi["180,-180"]
    elif "180,180" in keys:
        phipsi["-180,180"]=phipsi["180,180"]
        phipsi["180,-180"]=phipsi["180,180"]
        phipsi["-180,-180"]=phipsi["180,180"]

    
#def Ala_MapsGenerator_parallel(atoms,pdb_id,chain_id,resname,resnumber,prepro,frame=0):
def Ala_MapsGenerator_parallel(params_list):
    #atoms_to_rotate,"GLY",first_res+1,"normal",int(model),pdb_id,chain_id
    atoms=params_list[0]
    pdb_id=params_list[5]
    chain_id=params_list[6]
    resname=params_list[1]
    resnumber=params_list[2]
    prepro=params_list[3]
#    print params_list[4]
    frame=params_list[4]+19999
    #digimapfile=open("/home/ashraya/RMap_Work/2vb1_simulation/2vb1_one_fs/frame_maps/parallel_2vb1_A_"+resname+str(resnumber)+"_frame_"+str(frame)+".csv","w")
    digimapfile=open("/home/ashraya/RMap_Work/rmap_dynamics/2vb1_charmm/frame_maps/2vb1_"+resname+str(resnumber)+"_frame_"+str(frame)+".csv","w")
    framephi=FindDihedralAngle(atoms['C1'],atoms['N2'],atoms['CA2'],atoms['C2'])
    framepsi=FindDihedralAngle(atoms['N2'],atoms['CA2'],atoms['C2'],atoms['N3'])
    BringToZero(atoms,prepro,framephi,framepsi)
    phipsi=dict()
    for i in range(-180,180):
        for j in range(-180,180):
            finalres=''
            #bring psi to 0
            dircos=FindDirCosines(atoms['CA2'],atoms['C2'])
            atoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['O2'][0],atoms['O2'][1],atoms['O2'][2],1)
            atoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['N3'][0],atoms['N3'][1],atoms['N3'][2],1)
            if prepro=="normal":
                atoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['H3'][0],atoms['H3'][1],atoms['H3'][2],1)
            else:
                atoms['CD3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['CD3'][0],atoms['CD3'][1],atoms['CD3'][2],1)
            atoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['CA3'][0],atoms['CA3'][1],atoms['CA3'][2],1)
            phi=FindDihedralAngle(atoms['C1'],atoms['N2'],atoms['CA2'],atoms['C2'])
            psi=FindDihedralAngle(atoms['N2'],atoms['CA2'],atoms['C2'],atoms['N3'])
            #sampleout.writelines("\nphi="+str(phi)+"  psi="+str(psi)+"\n")
            if prepro=="normal":            
                res=CheckShortContacts(atoms['CA1'],atoms['H3'],2.40,2.20)
            else:
                res=CheckShortContacts(atoms['CA1'],atoms['CD3'],3.00,2.90)
            finalres=finalres+res[0]
            #if res[0]=='p' or res[0]=='d':
                #sampleout.write('CA1...H3 : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['HA'],2.40,2.20)
            finalres=finalres+res[0]
            #if res[0]=='p' or res[0]=='d':
                #sampleout.write('C1...HA : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['C2'],3.00,2.90)
            finalres=finalres+res[0]
            #if res[0]=='p' or res[0]=='d':
                #sampleout.write('C1...C2 : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['O2'],2.80,2.70)
            finalres=finalres+res[0]
            #if res[0]=='p' or res[0]=='d':
                #sampleout.write('C1...O2 : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['N3'],2.90,2.80)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('C1...N3 : '+res+'\n')
            if prepro=="normal":  
                res=CheckShortContacts(atoms['C1'],atoms['H3'],2.40,2.20)
            else:
                res=CheckShortContacts(atoms['C1'],atoms['CD3'],3.00,2.90)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('C1...H3 : '+res+'\n')
            res=CheckShortContacts(atoms['C1'],atoms['CB'],3.20,3.00)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('C1...CB : '+res+'\n')
            
            res=CheckShortContacts(atoms['O1'],atoms['CB'],2.80,2.70)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...CB : '+res+'\n')
            res=CheckShortContacts(atoms['O1'],atoms['HA'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...HA : '+res+'\n')
            res=CheckShortContacts(atoms['O1'],atoms['C2'],2.80,2.70)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...C2 : '+res+'\n')    
            res=CheckShortContacts(atoms['O1'],atoms['O2'],2.70,2.60)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...O2 : '+res+'\n')  
            res=CheckShortContacts(atoms['O1'],atoms['N3'],2.70,2.60)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...N3 : '+res+'\n') 
            if prepro=="normal":  
                res=CheckShortContacts(atoms['O1'],atoms['H3'],2.40,2.20)
            else:
                res=CheckShortContacts(atoms['O1'],atoms['CD3'],2.80,2.70)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...H3 : '+res+'\n') 
            res=CheckShortContacts(atoms['O1'],atoms['CA3'],2.80,2.70)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('O1...CA3 : '+res+'\n')  

            res=CheckShortContacts(atoms['N2'],atoms['O2'],2.70,2.60)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('N2...O2 : '+res+'\n')
            res=CheckShortContacts(atoms['N2'],atoms['N3'],2.70,2.60)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('N2...N3 : '+res+'\n')  
            if prepro=="normal":  
                res=CheckShortContacts(atoms['N2'],atoms['H3'],2.40,2.20)
            else:
                res=CheckShortContacts(atoms['N2'],atoms['CD3'],2.90,2.80)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('N2...H3 : '+res+'\n')
            
            res=CheckShortContacts(atoms['H2'],atoms['HA'],2.00,1.90)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('H2...HA : '+res+'\n')   
            res=CheckShortContacts(atoms['H2'],atoms['C2'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('H2...C2 : '+res+'\n')
            res=CheckShortContacts(atoms['H2'],atoms['O2'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('H2...O2 : '+res+'\n')
            res=CheckShortContacts(atoms['H2'],atoms['N3'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('H2...N3 : '+res+'\n')
            if prepro=="normal":  
                res=CheckShortContacts(atoms['H2'],atoms['H3'],2.00,1.90)
            else:
                res=CheckShortContacts(atoms['H2'],atoms['CD3'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('H2...H3 : '+res+'\n')
            res=CheckShortContacts(atoms['H2'],atoms['CB'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('H2...CB : '+res+'\n')
            
            res=CheckShortContacts(atoms['CB'],atoms['O2'],2.80,2.70)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('CB...O2 : '+res+'\n')
            res=CheckShortContacts(atoms['CB'],atoms['N3'],2.90,2.80)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('CB...N3 : '+res+'\n')
            if prepro=="normal":  
                res=CheckShortContacts(atoms['CB'],atoms['H3'],2.40,2.20)
            else:
                res=CheckShortContacts(atoms['CB'],atoms['CD3'],3.20,3.00)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('CB...H3 : '+res+'\n')
            
            res=CheckShortContacts(atoms['HA'],atoms['O2'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('HA...O2 : '+res+'\n')
            
            res=CheckShortContacts(atoms['HA'],atoms['N3'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('HA...N3 : '+res+'\n')
            if prepro=="normal":  
                res=CheckShortContacts(atoms['HA'],atoms['H3'],2.00,1.90)
            else:
                res=CheckShortContacts(atoms['HA'],atoms['CD3'],2.40,2.20)
            finalres=finalres+res[0]
#            if res[0]=='p' or res[0]=='d':
#                sampleout.write('HA...H3 : '+res+'\n')
            
            if 'd' not in finalres and 'p' not in finalres:
#                sampleout.write("This conformation is fully allowed\n--------------------------------------------------------\n")
#                digimapfile.write(str(phi)+","+str(psi)+",1\n")
                phipsi[str(phi)+","+str(psi)]=1
            elif 'd' in finalres:
#                sampleout.write("This conformation is disallowed\n--------------------------------------------------------\n")
#                digimapfile.write(str(phi)+","+str(psi)+",0\n")
                phipsi[str(phi)+","+str(psi)]=0
            else:
#                sampleout.write("This conformation is partially allowed\n--------------------------------------------------------\n")
#                digimapfile.write(str(phi)+","+str(psi)+",2\n") 
                phipsi[str(phi)+","+str(psi)]=2
        
        
        dircos=FindDirCosines(atoms['N2'], atoms['CA2'])
        atoms['CB']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CB'][0],atoms['CB'][1],atoms['CB'][2],1)
        atoms['CA2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],1)
        atoms['HA']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['HA'][0],atoms['HA'][1],atoms['HA'][2],1)
        atoms['C2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['C2'][0],atoms['C2'][1],atoms['C2'][2],1)
        atoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['O2'][0],atoms['O2'][1],atoms['O2'][2],1)
        atoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['N3'][0],atoms['N3'][1],atoms['N3'][2],1)
        if prepro=="normal":
            atoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['H3'][0],atoms['H3'][1],atoms['H3'][2],1)
        else:
            atoms['CD3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CD3'][0],atoms['CD3'][1],atoms['CD3'][2],1)
        atoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CA3'][0],atoms['CA3'][1],atoms['CA3'][2],1)   
            
            
                    
    #sampleout.close()
    #digimapfile.write(str(actphi)+","+str(actpsi)+",3")
    add_missing_phipsi(phipsi)
    for i in range(-180,181):
        for j in range(-180,181):
            digimapfile.write(str(i)+","+str(j)+","+str(phipsi[str(i)+","+str(j)])+"\n")
    digimapfile.close()

def send_to_queue(atom):
#    wanted_residues=[29,103,115]
    wanted_residues=[35]
    global peptide_q
    peptide_q.append(atom)
    #print "atom appended: "+atom[0]+"_"+atom[4]+"_"+str(atom[5])
    if atom[0]=="CA":
        #print "CA atom"
        atom_names=[row[0] for row in peptide_q]
        #print atom_names
        ca_count=atom_names.count("CA")
        #print str(ca_count)
        if ca_count==3:
#        if len(peptide_q)==15:
            #print peptide_q
            first_res=peptide_q[0][5]
            last_res=peptide_q[-1][5]
            if last_res-first_res!=2:
                new_pep=[]
                for item in peptide_q:
                    if item[5]==last_res-1:
                        new_pep.append(item)
                    elif item[5]==last_res:
                        #print "Missing res"
                        if item[0] in ["CA","C","O"]:
                            new_pep.append(item)
                            #print item
                peptide_q=[]
                peptide_q = new_pep[:]
            elif last_res-first_res==2:
                first_res_res=[]
                second_res_res=[]
                third_res_res=[]
                for item in peptide_q:
                    if item[5]==first_res:
                        first_res_res.append(item[0])
                    elif item[5]==first_res+1:
                        second_res_res.append(item[0])
                        second_res_type=item[4]
                    else:
                        third_res_res.append(item[0])
                first_res_flag=second_res_flag=third_res_flag=0
                if "CA" in first_res_res and "C" in first_res_res and "O" in first_res_res:
                    first_res_flag=1
                if "N" in second_res_res and "H" in second_res_res and "CA" in second_res_res and "C" in second_res_res and "O" in second_res_res:
                    if second_res_type!="GLY":
                        if "CB" in second_res_res and "HA" in second_res_res:
                            second_res_flag=1
                    else:
                        if "HA1" in second_res_res and "HA2" in second_res_res:
                            second_res_flag=2
                if "N" in third_res_res and  "CA" in third_res_res:
                    if "H" in third_res_res:
                        third_res_flag=1
                    elif "CD" in third_res_res:
                        third_res_flag=2
                if second_res_flag==0:
                    if "CA" in second_res_res and "C" in second_res_res and "O" in second_res_res:
                        second_res_flag=-1
                if first_res_flag==1 and (second_res_flag==1 or second_res_flag==2) and (third_res_flag==1 or third_res_flag==2):
                    atoms_to_rotate=dict()
                    new_pep=[]
                    for item in peptide_q:
                        if item[5]==first_res:
                            if item[0]=="CA":
                                atoms_to_rotate["CA1"]=[item[1],item[2],item[3]]
                            elif item[0]=="C":
                                atoms_to_rotate["C1"]=[item[1],item[2],item[3]]
                            elif item[0]=="O":
                                atoms_to_rotate["O1"]=[item[1],item[2],item[3]]
                        elif item[5]==first_res+1:
                            if item[0]=="N":
                                atoms_to_rotate["N2"]=[item[1],item[2],item[3]]
                            elif item[0]=="H":
                                atoms_to_rotate["H2"]=[item[1],item[2],item[3]]
                            elif item[0]=="CA":
                                atoms_to_rotate["CA2"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                            elif item[0]=="CB":
                                atoms_to_rotate["CB"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                            elif item[0]=="HA":
                                atoms_to_rotate["HA"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                            elif item[0]=="HA1":
                                atoms_to_rotate["HA1"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                            elif item[0]=="HA2":
                                atoms_to_rotate["HA2"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                            elif item[0]=="C":
                                atoms_to_rotate["C2"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                            elif item[0]=="O":
                                atoms_to_rotate["O2"]=[item[1],item[2],item[3]]
                                new_pep.append(item)
                        elif item[5]==last_res:
                            if item[0]=="N":
                                atoms_to_rotate["N3"]=[item[1],item[2],item[3]]
                            elif item[0]=="H":
                                atoms_to_rotate["H3"]=[item[1],item[2],item[3]]
                            elif item[0]=="CA":
                                atoms_to_rotate["CA3"]=[item[1],item[2],item[3]]
                            elif item[0]=="CD":
                                atoms_to_rotate["CD3"]=[item[1],item[2],item[3]]
                            new_pep.append(item)
                    peptide_q=[]
                    peptide_q=new_pep[:]
                    if second_res_flag==1:
                        if third_res_flag==1 and first_res+1 in wanted_residues:
#                            ala_rotate(atoms_to_rotate,second_res_type,first_res+1,"normal")
                            atom_queue.append([atoms_to_rotate,second_res_type,first_res+1,"normal",int(model),pdb_id,chain_id])
                        elif third_res_flag==2  and first_res+1 in wanted_residues:
#                            ala_rotate(atoms_to_rotate,second_res_type,first_res+1,"prepro")
                            atom_queue.append([atoms_to_rotate,second_res_type,first_res+1,"prepro",int(model),pdb_id,chain_id])
                    else:
                        if third_res_flag==1  and first_res+1 in wanted_residues:
                            atom_queue.append([atoms_to_rotate,"GLY",first_res+1,"normal",int(model),pdb_id,chain_id])
#                            gly_rotate(atoms_to_rotate,"GLY",first_res+1,"normal")
                        elif third_res_flag==2  and first_res+1 in wanted_residues:
                            atom_queue.append([atoms_to_rotate,"GLY",first_res+1,"prepro",int(model),pdb_id,chain_id])
#                            gly_rotate(atoms_to_rotate,"GLY",first_res+1,"prepro")
                elif second_res_flag==-1 and third_res_flag==1:
                    new_pep=[]
                    for item in peptide_q:
                        if item[5]==first_res+1 and item[0] in ["CA","C","O","CB","HA"]:
                            new_pep.append(item)
                        elif item[5]==last_res and item[0] in ["CA","N","H"]:
                            new_pep.append(item)
                    peptide_q=[]
                    peptide_q=new_pep[:]
                elif first_res_flag==0 and second_res_flag==1 and third_res_flag==1:
                    new_pep=[]
                    for item in peptide_q:
                        if item[5]==first_res+1 or item[5]==last_res:
                            new_pep.append(item)
                    peptide_q=[]
                    peptide_q=new_pep[:]
                else:
                    peptide_q=[]
                    peptide_q.append(atom)
                                
                        
    
#pdbfile=open("/home/ashraya/RMap_Work/2vb1_simulation/2vb1_one_fs/2vb1_1_fs_pdb.pdb","r")
flag=0
atom_queue=[]
first_model=True
missing_res=[]
atom_lines=[]
pdb_id="2vb1"
peptide_q=[]
chain_id=" "
first_n=True
first_h=True
#pdbfile=open("/home/ashraya/RMap_Work/2vb1_simulation/2vb1_one_fs/2vb1_1_fs_0_to_12499.pdb","r")
pdbfile=open("/home/ashraya/RMap_Work/rmap_dynamics/2vb1_charmm/2vb1_20000_to_21000.pdb","r")
for line1 in pdbfile:
    if line1.startswith("MODEL"):
        if first_model:
            first_model=False
            continue
        else:
#            print line1[0:-1]
            model=line1.split()[1]
            
            model=int(model)
            if model>10000:
                atom_lines=[]
                new_atom_lines=[]
                atom_queue=[]
                break
#            elif model!=10000:
#                atom_lines=[]
#                new_atom_lines=[]
#                atom_queue=[]
#                continue
            atom_queue=[]
            print model
            new_atom_lines=cleanPdbFile(atom_lines,missing_res)
            for line in new_atom_lines:
                if line.startswith("ATOM"):
                    if chain_id in line[21]:
                        if first_n and line[12:16].strip()=="N":
                            first_n=False
                        elif first_h and line[12:16].strip() in ["H1","H2","H3","H"]:
                            first_h=False
                        elif line[12:16].strip()=="N":
                            atom_to_consider=["N",float(line[30:38]),float(line[38:46]),float(line[46:54]),line[17:20],int(line[22:26])]
                            send_to_queue(atom_to_consider)
                        elif line[12:16].strip() in ["CA","C","O","CB","H","HA","HA1","HA2"]:
                            atom_to_consider=[line[12:16].strip(),float(line[30:38]),float(line[38:46]),float(line[46:54]),line[17:20],int(line[22:26])]
                            send_to_queue(atom_to_consider)
                        elif line[12:16].strip() in ["HN"]:
                            atom_to_consider=["H",float(line[30:38]),float(line[38:46]),float(line[46:54]),line[17:20],int(line[22:26])]
                            send_to_queue(atom_to_consider)
                        if line[17:20]=="PRO" and line[12:16].strip()=="CD":
                            atom_to_consider=[line[12:16].strip(),float(line[30:38]),float(line[38:46]),float(line[46:54]),line[17:20],int(line[22:26])]
                            send_to_queue(atom_to_consider)
            atom_lines=[]
            new_atom_lines=[]
#            print atom_queue
            pool=mp.Pool(2)
            #print atom_queue
            dummy_results=pool.map(Ala_MapsGenerator_parallel,[inputs for inputs in atom_queue])
            pool.close()
#            atom_queue=[]
    elif line1.startswith('ATOM'):
        atom_lines.append(line1)
atom_lines=[]
new_atom_lines=[]
pdbfile.close()

