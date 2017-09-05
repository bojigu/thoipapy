import re
import os
import math
###in homotypic_data/daba/Bind/pdb  ZhelixPair gives interacting TMP homodimers

def homodimer_residue_closedist_calculate_from_complex(pathdict,set_,logging):
    helixpair_file=set_["homodimer_structure_helix_pair_infor"]
    #helixpair_file=r"/scratch2/zeng/homotypic_data/data/Bind/pdb/zHelixPair"
    helixpair_file_handle=open(helixpair_file,"r")
    for row in helixpair_file_handle:
        row_array=row.strip().split("\t")
        pdb=row_array[0]
        chain1=row_array[1][0:1]
        chain2=row_array[7][0:1]
        chain1_pdb_tm_start=row_array[4]
        chain1_pdb_tm_end=row_array[5]
        chain2_pdb_tm_start=row_array[10]
        chain2_pdb_tm_end=row_array[11]
        hash1arrayx={}
        hash1arrayy={}
        hash1arrayz={}
        hash2arrayx={}
        hash2arrayy={}
        hash2arrayz={}
        hashclosedist={}
        hash1CAarrayx={}
        hash1CAarrayy={}
        hash1CAarrayz={}
        hash2CAarrayx={}
        hash2CAarrayy={}
        hash2CAarrayz={}
        hashCAclosedist={}
        hash1CBarrayx={}
        hash1CBarrayy={}
        hash1CBarrayz={}
        hash2CBarrayx={}
        hash2CBarrayy={}
        hash2CBarrayz={}
        hashCBclosedist={}
        #################residue distance calculation###################################
        #if pdb == "1f88":
        pdbfile=r"/scratch2/zeng/homotypic_data/data/Bind/pdb/%s.pdb" %pdb
        closedist_output_file = r"/scratch2/zeng/homotypic_data/data/Bind/pdb/%s.clodis" %pdb
        closedist_output_file_handle = open(closedist_output_file, 'w')
        pdbfile_handle=open(pdbfile,"r")
        for row1 in pdbfile_handle:
            Atom_line=re.search("^ATOM",row1)
            if Atom_line:
                index=row1[6:11]
                x=row1[30:38]
                y=row1[38:46]
                z=row1[46:54]
                chain=row1[21:22]
                residue_num=row1[22:26]
                residue_name=row1[17:20]
                atom=row1[13:15]
                if (chain == chain1 and int(chain1_pdb_tm_start) <= int(residue_num) and int(residue_num) <= int(chain1_pdb_tm_end)):
                    kk='_'.join([pdb,chain,str(int(residue_num)),residue_name])
                    if atom=="CA":
                        hash1CAarrayx[kk]=x
                        hash1CAarrayy[kk]=y
                        hash1CAarrayz[kk]=z
                    if atom=="CB":
                        hash1CBarrayx[kk]=x
                        hash1CBarrayy[kk]=y
                        hash1CBarrayz[kk]=z
                    if kk not in hash1arrayx.keys():
                        hash1arrayx[kk]=x
                        hash1arrayy[kk]=y
                        hash1arrayz[kk]=z
                    else:
                        hash1arrayx[kk]='_'.join([hash1arrayx[kk],str(float(x))])
                        hash1arrayy[kk]='_'.join([hash1arrayy[kk],str(float(y))])
                        hash1arrayz[kk]='_'.join([hash1arrayz[kk],str(float(z))])
                if (chain == chain2 and int(chain2_pdb_tm_start) <= int(residue_num) and int(residue_num) <= int(chain2_pdb_tm_end)):
                    kk='_'.join([pdb,chain,str(int(residue_num)),residue_name])
                    if atom=="CA":
                        hash2CAarrayx[kk]=x
                        hash2CAarrayy[kk]=y
                        hash2CAarrayz[kk]=z
                    if atom=="CB":
                        hash2CBarrayx[kk]=x
                        hash2CBarrayy[kk]=y
                        hash2CBarrayz[kk]=z
                    if kk not in hash2arrayx.keys():
                        hash2arrayx[kk]=x
                        hash2arrayy[kk]=y
                        hash2arrayz[kk]=z
                    else:
                        hash2arrayx[kk]='_'.join([hash2arrayx[kk],str(float(x))])
                        hash2arrayy[kk]='_'.join([hash2arrayy[kk],str(float(y))])
                        hash2arrayz[kk]='_'.join([hash2arrayz[kk],str(float(z))])
        for key in hash1arrayx.keys():
            arr_xvalue1=hash1arrayx[key].strip().split('_')
            arr_yvalue1=hash1arrayy[key].strip().split('_')
            arr_zvalue1=hash1arrayz[key].strip().split('_')
            for key1 in hash2arrayx.keys():
                for i in range(len(arr_xvalue1)):
                    arr_xvalue2=hash2arrayx[key1].strip().split('_')
                    #print(arr_xvalue2)
                    arr_yvalue2=hash2arrayy[key1].strip().split('_')
                    arr_zvalue2=hash2arrayz[key1].strip().split('_')
                    for j in range(len(arr_xvalue2)):
                        dist=math.sqrt((float(arr_xvalue1[i])-float(arr_xvalue2[j]))**2+(float(arr_yvalue1[i])-float(arr_yvalue2[j]))**2+(float(arr_zvalue1[i])-float(arr_zvalue2[j]))**2)
                        if key not in hashclosedist.keys():
                            hashclosedist[key]=dist
                        else:
                            if dist<hashclosedist[key]:
                                hashclosedist[key]=dist

        for key in sorted(hash1CAarrayx):
            for key1 in sorted(hash2CAarrayx):
                cadist=math.sqrt((float(hash1CAarrayx[key])-float(hash2CAarrayx[key1]))**2+(float(hash1CAarrayy[key])-float(hash2CAarrayy[key1]))**2+(float(hash1CAarrayz[key])-float(hash2CAarrayz[key1]))**2)
                if key not in hashCAclosedist:
                    hashCAclosedist[key]=cadist
                else:
                    if cadist<hashCAclosedist[key]:
                        hashCAclosedist[key]=cadist

        for key in sorted(hash1CBarrayx):
            for key1 in sorted(hash2CBarrayx):
                cbdist=math.sqrt((float(hash1CBarrayx[key])-float(hash2CBarrayx[key1]))**2+(float(hash1CBarrayy[key])-float(hash2CBarrayy[key1]))**2+(float(hash1CBarrayz[key])-float(hash2CBarrayz[key1]))**2)
                if key not in hashCBclosedist:
                    hashCBclosedist[key]=cbdist
                else:
                    if cbdist<hashCBclosedist[key]:
                        hashCBclosedist[key]=cbdist

        for key in sorted(hashclosedist):
            #closedist_output_file_handle.write(str(key)+"\t"+str(hashclosedist[key])+"\n")
            try:
                if key in hashCBclosedist:
                    closedist_output_file_handle.write(str(key.strip().split('_')[0])+"\t"+str(key.strip().split('_')[1])+"\t"+str(key.strip().split('_')[2])+"\t"+str(key.strip().split('_')[3])+"\t"+str(hashclosedist[key])+"\t"+str(hashCAclosedist[key])+"\t"+str(hashCBclosedist[key])+"\n")
                else:
                    closedist_output_file_handle.write(str(key.strip().split('_')[0])+"\t"+str(key.strip().split('_')[1])+"\t"+str(key.strip().split('_')[2])+"\t"+str(key.strip().split('_')[3])+"\t"+str(hashclosedist[key])+"\t"+str(hashCAclosedist[key])+"\t"+"0"+"\n")
            except:
                print("residue" + "%s" %key + "without CB" )
        closedist_output_file_handle.close()