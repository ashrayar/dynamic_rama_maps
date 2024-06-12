# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 11:37:58 2019

@author: ashraya
"""
def round_down(number,divisor):
    if number>=0:
        low = number-(number%divisor)
        up=low+9
    else:
        low = number+(abs(number)%divisor)
        up=low-9
    return low,up
def round_up(number,divisor):
    return number+(number%divisor)

opfile=open("/home/ashraya/RMap_Work/2vb1_simulation/energy_calc/ser72_frame1670/potential_vaccuum_unrest.csv","w")
conf_file=open("/home/ashraya/RMap_Work/2vb1_simulation/energy_calc/ser72_frame1670/list_of_conformers.txt","r")
for line1 in conf_file:
    conform=line1.split(".")[0]
    phi,psi=conform.split("_")[1:]
    phi=int(phi)
    psi=int(psi)
    if phi==-180:
        folder_name="phi_-180"
    elif phi<=-1 and phi >=-9:
        folder_name="phi_-1_to_-9"
    else:
        lower_bound,upper_bound=round_down(phi,10)
        folder_name="phi"+"_"+str(lower_bound)+"_to_"+str(upper_bound)
    infile=open("/home/ashraya/RMap_Work/2vb1_simulation/energy_calc/ser72_frame1670/ser72_gromacs_outputs/"+folder_name+"/phipsi_"+str(phi)+"_"+str(psi)+"_unrest_vaccuum.log","r")
    flag=0        
    for line in infile:
        if "  Potential" in line:
            flag=1
            continue
        if flag==1:
            lineparts=line.split()
            potential=float(lineparts[-1])
            flag=0
            break
    infile.close()
    opfile.write(str(phi)+","+str(psi)+","+str(potential)+"\n")
opfile.close()
conf_file.close()
'''
mapfile1=open("/home/ashraya/RMap_Work/2vb1_simulation/energy_calc/ser72_frame1670/potential_vaccuum_unrest.csv","r")
phipsi=dict()
for line in mapfile1:
    phipsiparts=line.split(",")
    #print phipsiparts
    phi=int(phipsiparts[0])
    psi=int(phipsiparts[1])
    apd=float(phipsiparts[2])
    phipsi[str(phi)+","+str(psi)]=apd
mapfile1.close()
keys=phipsi.keys()
for phi in range(-180,181):
    if str(phi)+",-180" in keys:
        phipsi[str(phi)+",180"]=phipsi[str(phi)+",-180"]
    elif str(phi)+",180" in keys:
        phipsi[str(phi)+",-180"]=phipsi[str(phi)+",180"]
if "-180,180" in keys:
    phipsi["180,180"]=phipsi["-180,-180"]=phipsi["180,-180"]=phipsi["-180,180"]
elif "-180,-180" in keys:
    phipsi["180,180"]=phipsi["180,-180"]=phipsi["-180,180"]=phipsi["-180,-180"]
elif "180,-180" in keys:
    phipsi["180,180"]=phipsi["-180,180"]=phipsi["-180,-180"]=phipsi["180,-180"]
elif "180,180" in keys:
    phipsi["-180,180"]=phipsi["-180,-180"]=phipsi["180,-180"]=phipsi["180,180"]
for psi in range(-180,181):
    phipsi["180,"+str(psi)]=phipsi["-180,"+str(psi)]
opfile=open("/home/ashraya/RMap_Work/2vb1_simulation/energy_calc/ser72_frame1670/potential_vaccuum_unrest_full.csv","w")
for i in range(-180,181):
    for j in range(-180,181):
        opfile.write(str(i)+","+str(j)+","+str(phipsi[str(i)+","+str(j)])+"\n")
opfile.close()
'''     
        