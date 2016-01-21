#!/usr/bin/env python


# TEST

import sys

residues = []
temp = False
T = ""
filename = ""
opdracht = ""

class Model:
    def __init__(self, number):
        self.number = number
        self.res = []

class Residue:
    def __init__(self, res_number):
#        self.res_type = res_type
        self.res_number = int(res_number)
        self.atoms = []
        self.CA = Atom(0, 0)
#        self.AvgB = 0
        
class Atom:
#    def __init__(self, name, b_factor): # use *args, **kwargs to choose between b_factor or x,y,z?
#    def __init__(self, number, name, x, y, z, occupancy, b_factor):
#    def __init__(self, name, *args, **kwargs):
        self.number = int(number)
#        self.name = name
#        self.chain = chain
#        self.x = float(x)
#        self.y = float(y)
#        self.z = float(z)
#        self.occupancy = float(occupancy)
#        self.b_factor = float(b_factor)

#  This method writes all B-factors from the C-alpha's to a file
def Return_CA_Bfactors():
    global filename
    bestand = open(filename[0:-4] +'_CA.txt', 'w')
    if temp:
        bestand.write("# Experimental conditions - Temperature : %s"%T)
    for r in residues:
        bestand.write(str(r.res_number) + " " + str(r.CA.b_factor) + "\n")
    bestand.close()

# This method first calculates, and then writes the average B-factor per residue to a file
def ReturnAvgBfactor():
    global filename
    bestand = open(filename[0:-4] + '_AvgB.txt', 'w')
    for r in residues:
        total = float(0)
        t = int(0)
        for a in r.atoms:
            total += a.b_factor
            t += 1
        AvgB = total/t
        bestand.write(str(r.res_number) + " " + str(AvgB) + "\n")
    bestand.close()
    
# This method calculates the total average B-factor and writes it to the terminal
def ReturnTotalAvgB():
    for r in residues:
        total = float(0)
        nres = 0
        for a in r.atoms:
            if "H" not in a.name and "O" not in a.name:
                total += a.b_factor
                nres += 1
    tavg = float(total/nres)
    print("The average b-factor (without considering hydrogens) is: %s" %tavg)
    
def main():
    global filename
    global opdracht
    print("What PDB-file do you want to be read?")
    filename = sys.stdin.readline()[:-1]
    print("Do you want to extract b-factors? y/n")
    if sys.stdin.readline() == "y"
        print("Do you want the average b-factor per residue (avg), the c-alpha b-factor per residue (ca) or the total average b-factor (tavg)?")
        opdracht = sys.stdin.readline()[:-1]
#    temperature = sys.argv[3]
#    if temperature == "temp":
#        temp = True
    LeesFile()
    if opdracht == "ca":
        Return_CA_Bfactors()
    if opdracht == "avg":
        ReturnAvgBfactor()
    if opdracht == "tavg":
        ReturnTotalAvgB()
    else:
        # doe iets met symmetrie
        
def LeesFile():
    global filename
    f = open(filename, 'r')
    l = f.readline()
    if temp:
        while "REMARK 200  TEMPERATURE           (KELVIN) :" not in l:
            l = f.readline()
        words = l.split()
        T = words[5]
    current_mdl = Model(0)
    while l[0:4] != "ATOM":
        if l[0:6] == "MODEL ":
            mdl_number = l[10:14]
            current_mdl = Model(int(mdl_number))
        l = f.readline()
    while l[0:3] != "TER":
        if l[0:6] == "ENDMDL":
            while l[0:6] != "MODEL ":
                l = f.readline()
            if l[0:3] != "TER":
                break
        while l[0:6] != "ENDMDL":
            if l[0:4] == "ATOM":
                res_number = l[22:26]
                current_residue = Residue(res_number)
    #            current_residue = Residue(l[17:20], res_number)
                current_mdl.res.append(current_residue)
                while l[0:3] != "TER" and l[22:26] == res_number:
                    if l[0:4] == "ATOM":
                        current_atom = Atom(l[12:16], l[60:66]) # change this to the right values for number, x, y and z
    #                    current_atom = Atom(l[6:11], l[12:16], l[21], l[30:38], l[38:46], l[46:54], l[54:60], l[60:66])
                        current_residue.atoms.append(current_atom)
                        if "CA" in current_atom.name:
                            current_residue.CA = current_atom
                    l = f.readline()
            else:
                l = f.readline()
    f.close()
    
if __name__ == "__main__":
    main()