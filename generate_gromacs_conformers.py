# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 18:17:53 2019

@author: ashraya
"""

import math
import numpy as np
import os
#from CalculatePhiPsi import FindDihedralAngle
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
    chi=math.acos(round(Q1dotQ2,10))
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

def writePDB(atomname):
    writestr=""
    if atomname=="CA1":
        writestr+="ATOM      1  CH3 ACE A   34   "
    elif atomname=="C1":
        writestr+="ATOM      2   C  ACE A   34   "
    elif atomname=="O1":
        writestr+="ATOM      3   O  ACE A   34   "
    elif atomname=="HC11":
        writestr+="ATOM      4 HH31 ACE A   34   "
    elif atomname=="HC12":
        writestr+="ATOM      5 HH32 ACE A   34   "
    elif atomname=="HC13":
        writestr+="ATOM      6 HH33 ACE A   34   "
    elif atomname=="N2":
        writestr+="ATOM      7   N  ALA A   35   "  
    elif atomname=="H2":
        writestr+="ATOM      8   H  ALA A   35   "
    elif atomname=="CB":
        writestr+="ATOM      9   CB ALA A   35   "
    elif atomname=="HB1":
        writestr+="ATOM     10  HB1 ALA A   35   "   
    elif atomname=="HB2":
        writestr+="ATOM     11  HB2 ALA A   35   "   
    elif atomname=="HB3":
        writestr+="ATOM     12  HB3 ALA A   35   "   
    elif atomname=="HA":
        writestr+="ATOM     13   HA ALA A   35   "
    elif atomname=="CA2":
        writestr+="ATOM     14   CA ALA A   35   "
    elif atomname=="C2":
        writestr+="ATOM     15   C  ALA A   35   "
    elif atomname=="O2":
        writestr+="ATOM     16   O  ALA A   35   "
    elif atomname=="N3":
        writestr+="ATOM     17   N  NME A   36   "  
    elif atomname=="H3":
        writestr+="ATOM     18   H  NME A   36   "
    elif atomname=="CA3":
        writestr+="ATOM     19  CH3 NME A   36   "   
    elif atomname=="HC31":
        writestr+="ATOM     20 HH31 NME A   36   "
    elif atomname=="HC32":
        writestr+="ATOM     21 HH32 NME A   36   "
    elif atomname=="HC33":
        writestr+="ATOM     22 HH33 NME A   36   "
    
    for i in range(8-len(str(round(atoms[atomname][0],3)))):
        writestr+=" "
    writestr+=str(round(atoms[atomname][0],3))
    for i in range(8-len(str(round(atoms[atomname][1],3)))):
        writestr+=" "
    writestr+=str(round(atoms[atomname][1],3))
    for i in range(8-len(str(round(atoms[atomname][2],3)))):
        writestr+=" "
    writestr+=str(round(atoms[atomname][2],3))
    writestr+="  1.00"
    writestr+="  0.00"
    writestr+="           "+atomname[0]+" \n"
    return writestr
    

def BringToZero(tempatoms):
    phi=FindDihedralAngle(tempatoms['C1'],tempatoms['N2'],tempatoms['CA2'],tempatoms['C2'])
    psi=FindDihedralAngle(tempatoms['N2'],tempatoms['CA2'],tempatoms['C2'],tempatoms['N3'])
#    print phi
#    print psi
    #bring phi to 0
    dircos=FindDirCosines(tempatoms['N2'], tempatoms['CA2'])
    atoms['CB']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CB'][0],tempatoms['CB'][1],tempatoms['CB'][2],-phi)
    atoms['HB1']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['HB1'][0],tempatoms['HB1'][1],tempatoms['HB1'][2],-phi)    
    atoms['HB2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['HB2'][0],tempatoms['HB2'][1],tempatoms['HB2'][2],-phi)    
    atoms['HB3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['HB3'][0],tempatoms['HB3'][1],tempatoms['HB3'][2],-phi)        
    atoms['CA2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],-phi)
    atoms['HA']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['HA'][0],tempatoms['HA'][1],tempatoms['HA'][2],-phi)
    atoms['C2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['C2'][0],tempatoms['C2'][1],tempatoms['C2'][2],-phi)
    tempatoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['O2'][0],tempatoms['O2'][1],tempatoms['O2'][2],-phi)
    tempatoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['N3'][0],tempatoms['N3'][1],tempatoms['N3'][2],-phi)
    tempatoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['H3'][0],tempatoms['H3'][1],tempatoms['H3'][2],-phi)
    tempatoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['CA3'][0],tempatoms['CA3'][1],tempatoms['CA3'][2],-phi)
    tempatoms['HC31']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['HC31'][0],tempatoms['HC31'][1],tempatoms['HC31'][2],-phi)   
    tempatoms['HC32']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['HC32'][0],tempatoms['HC32'][1],tempatoms['HC32'][2],-phi)   
    tempatoms['HC33']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['N2'][0],tempatoms['N2'][1],tempatoms['N2'][2],tempatoms['HC33'][0],tempatoms['HC33'][1],tempatoms['HC33'][2],-phi)   

    #bring psi to 0
    dircos=FindDirCosines(tempatoms['CA2'],tempatoms['C2'])
    atoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['O2'][0],tempatoms['O2'][1],tempatoms['O2'][2],-psi)
    atoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['N3'][0],tempatoms['N3'][1],tempatoms['N3'][2],-psi)
    atoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['H3'][0],tempatoms['H3'][1],tempatoms['H3'][2],-psi)
    atoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['CA3'][0],tempatoms['CA3'][1],tempatoms['CA3'][2],-psi)
    atoms['HC31']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['HC31'][0],tempatoms['HC31'][1],tempatoms['HC31'][2],-psi)   
    atoms['HC32']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['HC32'][0],tempatoms['HC32'][1],tempatoms['HC32'][2],-psi)   
    atoms['HC33']=resultantpoint(dircos[0],dircos[1],dircos[2],tempatoms['CA2'][0],tempatoms['CA2'][1],tempatoms['CA2'][2],tempatoms['HC33'][0],tempatoms['HC33'][1],tempatoms['HC33'][2],-psi)   
    
def MapsGenerator(peptide):
    BringToZero(atoms)
#    print atoms
    phi=FindDihedralAngle(atoms['C1'],atoms['N2'],atoms['CA2'],atoms['C2'])
    psi=FindDihedralAngle(atoms['N2'],atoms['CA2'],atoms['C2'],atoms['N3'])
#    print phi
#    print psi
    opfile=open("/home/ashraya/RMap_Work/2vb1_simulation/energy_calc/ser72_frame1671/ser72_conformers/phipsi_0_0.pdb","w")
    opfile.write("REMARK phi=0 psi=0\n")
    opfile.write(writePDB("CA1"))
    opfile.write(writePDB("HC11"))
    opfile.write(writePDB("HC12"))
    opfile.write(writePDB("HC13"))
    opfile.write(writePDB("C1"))
    opfile.write(writePDB("O1"))
    opfile.write(writePDB("N2"))
    opfile.write(writePDB("H2"))
    opfile.write(writePDB("CA2"))
    opfile.write(writePDB("CB"))
    opfile.write(writePDB("HA"))
    opfile.write(writePDB("HB1"))
    opfile.write(writePDB("HB2"))
    opfile.write(writePDB("HB3"))
    opfile.write(writePDB("C2"))
    opfile.write(writePDB("O2"))
    opfile.write(writePDB("N3"))
    opfile.write(writePDB("H3"))
    opfile.write(writePDB("CA3"))
    opfile.write(writePDB("HC31"))
    opfile.write(writePDB("HC32"))
    opfile.write(writePDB("HC33"))
    opfile.close()
    count=2
    
    for i in range(-180,180):
        for j in range(-180,180):
            #bring psi to 0
            dircos=FindDirCosines(atoms['CA2'],atoms['C2'])
            atoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['O2'][0],atoms['O2'][1],atoms['O2'][2],1)
            atoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['N3'][0],atoms['N3'][1],atoms['N3'][2],1)
            atoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['H3'][0],atoms['H3'][1],atoms['H3'][2],1)
            atoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['CA3'][0],atoms['CA3'][1],atoms['CA3'][2],1)
            atoms['HC31']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['HC31'][0],atoms['HC31'][1],atoms['HC31'][2],1)   
            atoms['HC32']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['HC32'][0],atoms['HC32'][1],atoms['HC32'][2],1)   
            atoms['HC33']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],atoms['HC33'][0],atoms['HC33'][1],atoms['HC33'][2],1)   
            phi=FindDihedralAngle(atoms['C1'],atoms['N2'],atoms['CA2'],atoms['C2'])
            psi=FindDihedralAngle(atoms['N2'],atoms['CA2'],atoms['C2'],atoms['N3'])
            #print str(phi)+"\t"+str(psi)
            phi=int(round(phi))
            psi=int(round(psi))
            opfile=open("/home/ashraya/RMap_Work/2vb1_simulation/energy_calc/ser72_frame1671/ser72_conformers/phipsi_"+str(phi)+"_"+str(psi)+".pdb","w")
            count+=1
            opfile.write("REMARK phi="+str(phi)+" psi="+str(psi)+"\n")
            opfile.write(writePDB("CA1"))
            opfile.write(writePDB("HC11"))
            opfile.write(writePDB("HC12"))
            opfile.write(writePDB("HC13"))
            opfile.write(writePDB("C1"))
            opfile.write(writePDB("O1"))
            opfile.write(writePDB("N2"))
            opfile.write(writePDB("H2"))
            opfile.write(writePDB("CA2"))
            opfile.write(writePDB("CB"))
            opfile.write(writePDB("HA"))
            opfile.write(writePDB("HB1"))
            opfile.write(writePDB("HB2"))
            opfile.write(writePDB("HB3"))
            opfile.write(writePDB("C2"))
            opfile.write(writePDB("O2"))
            opfile.write(writePDB("N3"))
            opfile.write(writePDB("H3"))
            opfile.write(writePDB("CA3"))
            opfile.write(writePDB("HC31"))
            opfile.write(writePDB("HC32"))
            opfile.write(writePDB("HC33"))
            opfile.close()

        dircos=FindDirCosines(atoms['N2'], atoms['CA2'])
        atoms['CB']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CB'][0],atoms['CB'][1],atoms['CB'][2],1)
        atoms['HB1']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['HB1'][0],atoms['HB1'][1],atoms['HB1'][2],1)
        atoms['HB2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['HB2'][0],atoms['HB2'][1],atoms['HB2'][2],1)
        atoms['HB3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['HB3'][0],atoms['HB3'][1],atoms['HB3'][2],1)
        atoms['CA2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CA2'][0],atoms['CA2'][1],atoms['CA2'][2],1)
        atoms['HA']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['HA'][0],atoms['HA'][1],atoms['HA'][2],1)
        atoms['C2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['C2'][0],atoms['C2'][1],atoms['C2'][2],1)
        atoms['O2']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['O2'][0],atoms['O2'][1],atoms['O2'][2],1)
        atoms['N3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['N3'][0],atoms['N3'][1],atoms['N3'][2],1)
        atoms['H3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['H3'][0],atoms['H3'][1],atoms['H3'][2],1)
        atoms['CA3']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['CA3'][0],atoms['CA3'][1],atoms['CA3'][2],1)   
        atoms['HC31']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['HC31'][0],atoms['HC31'][1],atoms['HC31'][2],1)   
        atoms['HC32']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['HC32'][0],atoms['HC32'][1],atoms['HC32'][2],1)   
        atoms['HC33']=resultantpoint(dircos[0],dircos[1],dircos[2],atoms['N2'][0],atoms['N2'][1],atoms['N2'][2],atoms['HC33'][0],atoms['HC33'][1],atoms['HC33'][2],1)   

        phi=FindDihedralAngle(atoms['C1'],atoms['N2'],atoms['CA2'],atoms['C2'])
        psi=FindDihedralAngle(atoms['N2'],atoms['CA2'],atoms['C2'],atoms['N3'])
        #print str(phi)+"\t"+str(psi)
        phi=int(round(phi))
        psi=int(round(psi))
        opfile=open("/home/ashraya/RMap_Work/2vb1_simulation/energy_calc/ser72_frame1671/ser72_conformers/phipsi_"+str(phi)+"_"+str(psi)+".pdb","w")
        count+=1
        opfile.write("REMARK phi="+str(phi)+" psi="+str(psi)+"\n")
        opfile.write(writePDB("CA1"))
        opfile.write(writePDB("HC11"))
        opfile.write(writePDB("HC12"))
        opfile.write(writePDB("HC13"))
        opfile.write(writePDB("C1"))
        opfile.write(writePDB("O1"))
        opfile.write(writePDB("N2"))
        opfile.write(writePDB("H2"))
        opfile.write(writePDB("CA2"))
        opfile.write(writePDB("CB"))
        opfile.write(writePDB("HA"))
        opfile.write(writePDB("HB1"))
        opfile.write(writePDB("HB2"))
        opfile.write(writePDB("HB3"))
        opfile.write(writePDB("C2"))
        opfile.write(writePDB("O2"))
        opfile.write(writePDB("N3"))
        opfile.write(writePDB("H3"))
        opfile.write(writePDB("CA3"))
        opfile.write(writePDB("HC31"))
        opfile.write(writePDB("HC32"))
        opfile.write(writePDB("HC33"))
        opfile.close()
                    
#    opfile.close()
    #digimapfile.write(str(actphi)+","+str(actpsi)+",3")
    
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


atoms=dict()
ipfile=open("/home/ashraya/RMap_Work/2vb1_simulation/energy_calc/ser72_frame1671/2vb1_ser72_h.pdb","r")
for line in ipfile:
    if "ATOM" in line:
        lineparts=line.split()
        xcoord=float(lineparts[5])
        ycoord=float(lineparts[6])
        zcoord=float(lineparts[7])
        atoms[lineparts[2]]=[xcoord,ycoord,zcoord]  
#        print lineparts[2]
ipfile.close()
MapsGenerator("AAtAA1")

        
        