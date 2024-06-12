# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 15:17:09 2020

@author: ashraya
"""
'''
#infile=open("/home/ashraya/RMap_Work/2vb1_simulation/2vb1_one_fs/2vb1_1_fs_nowater.gro","r")
infile=open("/home/ashraya/RMap_Work/rmap_dynamics/2vb1_charmm/2vb1_no_water.gro","r")
opfile=open("/home/ashraya/RMap_Work/rmap_dynamics/2vb1_charmm/2vb1_vmd_ba_script.tcl","w")
prev_res=1
atoms_to_consider=['C','CA','CB','HN','N','O','HA']
atom_index=dict()
for line in infile:
    resno=int(line[0:5])
    resname=line[5:8].strip()
    if resno==prev_res:
       atom_name=line[10:15].strip()
       if atom_name in atoms_to_consider:
           atom_index[atom_name]=int(line[15:20])-1
       prev_res=resno
       prev_resname=resname
    elif resno!=prev_res:
       print resno
       atom_name=line[10:15].strip()
       if atom_name=="N":
           atom_index["N2"]=int(line[15:20])-1
       opfile.write('set output [open "'+str(prev_res)+'_N_CA_bond.txt" w]\n')
       opfile.write('puts $output [measure bond {'+str(atom_index["N"])+' '+str(atom_index["CA"])+'} first 1]\n')
       opfile.write('close $output\n')
       if prev_res!=1 and prev_resname!="PRO":
           opfile.write('set output [open "'+str(prev_res)+'_N_H_bond.txt" w]\n')
           opfile.write('puts $output [measure bond {'+str(atom_index["N"])+' '+str(atom_index["HN"])+'} first 1]\n')
           opfile.write('close $output\n')
       if prev_resname!="GLY":
           opfile.write('set output [open "'+str(prev_res)+'_CA_CB_bond.txt" w]\n')
           opfile.write('puts $output [measure bond {'+str(atom_index["CA"])+' '+str(atom_index["CB"])+'} first 1]\n')
           opfile.write('close $output\n')
           opfile.write('set output [open "'+str(prev_res)+'_CA_HA_bond.txt" w]\n')
           opfile.write('puts $output [measure bond {'+str(atom_index["CA"])+' '+str(atom_index["HA"])+'} first 1]\n')
           opfile.write('close $output\n')
       opfile.write('set output [open "'+str(prev_res)+'_CA_C_bond.txt" w]\n')
       opfile.write('puts $output [measure bond {'+str(atom_index["CA"])+' '+str(atom_index["C"])+'} first 1]\n')
       opfile.write('close $output\n')
       opfile.write('set output [open "'+str(prev_res)+'_C_O_bond.txt" w]\n')
       opfile.write('puts $output [measure bond {'+str(atom_index["C"])+' '+str(atom_index["O"])+'} first 1]\n')
       opfile.write('close $output\n')
       opfile.write('set output [open "'+str(prev_res)+'_C_N_bond.txt" w]\n')
       opfile.write('puts $output [measure bond {'+str(atom_index["C"])+' '+str(atom_index["N2"])+'} first 1]\n')
       opfile.write('close $output\n')
       if prev_res!=1 and prev_resname!="PRO":
           opfile.write('set output [open "'+str(prev_res)+'_H_N_CA_angle.txt" w]\n')
           opfile.write('puts $output [measure angle {'+str(atom_index["HN"])+' '+str(atom_index["N"])+' '+str(atom_index["CA"])+'} first 1]\n')
           opfile.write('close $output\n')
       opfile.write('set output [open "'+str(prev_res)+'_N_CA_C_angle.txt" w]\n')
       opfile.write('puts $output [measure angle {'+str(atom_index["N"])+' '+str(atom_index["CA"])+' '+str(atom_index["C"])+'} first 1]\n')
       opfile.write('close $output\n')
       if prev_resname!="GLY":
           opfile.write('set output [open "'+str(prev_res)+'_N_CA_CB_angle.txt" w]\n')
           opfile.write('puts $output [measure angle {'+str(atom_index["N"])+' '+str(atom_index["CA"])+' '+str(atom_index["CB"])+'} first 1]\n')
           opfile.write('close $output\n')
           opfile.write('set output [open "'+str(prev_res)+'_CB_CA_HA_angle.txt" w]\n')
           opfile.write('puts $output [measure angle {'+str(atom_index["CB"])+' '+str(atom_index["CA"])+' '+str(atom_index["HA"])+'} first 1]\n')
           opfile.write('close $output\n')
           opfile.write('set output [open "'+str(prev_res)+'_N_CA_HA_angle.txt" w]\n')
           opfile.write('puts $output [measure angle {'+str(atom_index["N"])+' '+str(atom_index["CA"])+' '+str(atom_index["HA"])+'} first 1]\n')
           opfile.write('close $output\n')
           opfile.write('set output [open "'+str(prev_res)+'_CB_CA_C_angle.txt" w]\n')
           opfile.write('puts $output [measure angle {'+str(atom_index["CB"])+' '+str(atom_index["CA"])+' '+str(atom_index["C"])+'} first 1]\n')
           opfile.write('close $output\n')
           opfile.write('set output [open "'+str(prev_res)+'_HA_CA_C_angle.txt" w]\n')
           opfile.write('puts $output [measure angle {'+str(atom_index["HA"])+' '+str(atom_index["CA"])+' '+str(atom_index["C"])+'} first 1]\n')
           opfile.write('close $output\n')
       opfile.write('set output [open "'+str(prev_res)+'_CA_C_O_angle.txt" w]\n')
       opfile.write('puts $output [measure angle {'+str(atom_index["CA"])+' '+str(atom_index["C"])+' '+str(atom_index["O"])+'} first 1]\n')
       opfile.write('close $output\n')
       opfile.write('set output [open "'+str(prev_res)+'_CA_C_N_angle.txt" w]\n')
       opfile.write('puts $output [measure angle {'+str(atom_index["CA"])+' '+str(atom_index["C"])+' '+str(atom_index["N2"])+'} first 1]\n')
       opfile.write('close $output\n')
       opfile.write('set output [open "'+str(prev_res)+'_O_C_N_angle.txt" w]\n')
       opfile.write('puts $output [measure angle {'+str(atom_index["O"])+' '+str(atom_index["C"])+' '+str(atom_index["N2"])+'} first 1]\n')
       opfile.write('close $output\n')
       atom_index=dict()
       if atom_name in atoms_to_consider:
           atom_index[atom_name]=int(line[15:20])-1
       prev_res=resno
       prev_resname=resname
infile.close()

opfile.write('set output [open "'+str(prev_res)+'_N_CA_bond.txt" w]\n')
opfile.write('puts $output [measure bond {'+str(atom_index["N"])+' '+str(atom_index["CA"])+'} first 1]\n')
opfile.write('close $output\n')
opfile.write('set output [open "'+str(prev_res)+'_N_H_bond.txt" w]\n')
opfile.write('puts $output [measure bond {'+str(atom_index["N"])+' '+str(atom_index["HN"])+'} first 1]\n')
opfile.write('close $output\n')
opfile.write('set output [open "'+str(prev_res)+'_CA_CB_bond.txt" w]\n')
opfile.write('puts $output [measure bond {'+str(atom_index["CA"])+' '+str(atom_index["CB"])+'} first 1]\n')
opfile.write('close $output\n')
opfile.write('set output [open "'+str(prev_res)+'_CA_HA_bond.txt" w]\n')
opfile.write('puts $output [measure bond {'+str(atom_index["CA"])+' '+str(atom_index["HA"])+'} first 1]\n')
opfile.write('close $output\n')
opfile.write('set output [open "'+str(prev_res)+'_CA_C_bond.txt" w]\n')
opfile.write('puts $output [measure bond {'+str(atom_index["CA"])+' '+str(atom_index["C"])+'} first 1]\n')
opfile.write('close $output\n')
if prev_res!=1 and prev_resname!="PRO":
   opfile.write('set output [open "'+str(prev_res)+'_H_N_CA_angle.txt" w]\n')
   opfile.write('puts $output [measure angle {'+str(atom_index["HN"])+' '+str(atom_index["N"])+' '+str(atom_index["CA"])+'} first 1]\n')
   opfile.write('close $output\n')
opfile.write('set output [open "'+str(prev_res)+'_N_CA_C_angle.txt" w]\n')
opfile.write('puts $output [measure angle {'+str(atom_index["N"])+' '+str(atom_index["CA"])+' '+str(atom_index["C"])+'} first 1]\n')
opfile.write('close $output\n')
if prev_resname!="GLY":
   opfile.write('set output [open "'+str(prev_res)+'_N_CA_CB_angle.txt" w]\n')
   opfile.write('puts $output [measure angle {'+str(atom_index["N"])+' '+str(atom_index["CA"])+' '+str(atom_index["CB"])+'} first 1]\n')
   opfile.write('close $output\n')
   opfile.write('set output [open "'+str(prev_res)+'_CB_CA_HA_angle.txt" w]\n')
   opfile.write('puts $output [measure angle {'+str(atom_index["CB"])+' '+str(atom_index["CA"])+' '+str(atom_index["HA"])+'} first 1]\n')
   opfile.write('close $output\n')
   opfile.write('set output [open "'+str(prev_res)+'_N_CA_HA_angle.txt" w]\n')
   opfile.write('puts $output [measure angle {'+str(atom_index["N"])+' '+str(atom_index["CA"])+' '+str(atom_index["HA"])+'} first 1]\n')
   opfile.write('close $output\n')
   opfile.write('set output [open "'+str(prev_res)+'_CB_CA_C_angle.txt" w]\n')
   opfile.write('puts $output [measure angle {'+str(atom_index["CB"])+' '+str(atom_index["CA"])+' '+str(atom_index["C"])+'} first 1]\n')
   opfile.write('close $output\n')
   opfile.write('set output [open "'+str(prev_res)+'_HA_CA_C_angle.txt" w]\n')
   opfile.write('puts $output [measure angle {'+str(atom_index["HA"])+' '+str(atom_index["CA"])+' '+str(atom_index["C"])+'} first 1]\n')
   opfile.write('close $output\n')
opfile.close()
'''

#'''      
#infile=open("/home/ashraya/RMap_Work/2vb1_simulation/2vb1_one_fs/2vb1_1_fs_nowater.gro","r")
infile=open("/home/ashraya/RMap_Work/rmap_dynamics/2vb1_charmm/2vb1_no_water.gro","r")
opfile=open("/home/ashraya/RMap_Work/rmap_dynamics/2vb1_charmm/2vb1_vmd_ba_remaining_script.tcl","w")
atom_index=dict()
first=True
resnum=1
resname="LYS"
for line in infile:
    cur_resnum=int(line[0:5])
    cur_resname=line[5:8].strip()       
    atom_name=line[10:15].strip()
    if atom_name=="C":
        if first:
            first=False
            atom_index[atom_name]=int(line[15:20])-1
            continue
        else:
            print resnum
            opfile.write('set output [open "'+str(resnum)+'_C_N_CA_angle.txt" w]\n')
            opfile.write('puts $output [measure angle {'+str(atom_index["C"])+' '+str(atom_index["N"])+' '+str(atom_index["CA"])+'} first 1]\n')
            opfile.write('close $output\n')
            if cur_resname!="PRO":
                opfile.write('set output [open "'+str(resnum)+'_C_N_H_angle.txt" w]\n')
                opfile.write('puts $output [measure angle {'+str(atom_index["C"])+' '+str(atom_index["N"])+' '+str(atom_index["HN"])+'} first 1]\n')
                opfile.write('close $output\n')
            atom_index=dict()
            resnum=cur_resnum
            resname=cur_resname
            atom_index[atom_name]=int(line[15:20])-1
    elif atom_name in ["N","CA","HN"]:
        if cur_resnum==resnum+1:
            atom_index[atom_name]=int(line[15:20])-1
infile.close()
opfile.close()
#'''
            