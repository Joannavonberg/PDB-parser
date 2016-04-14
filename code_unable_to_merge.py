#!/usr/bin/env python

import sys

models = []
temp = False
T = ""
filename = ""
opdracht = ""

class Model:
    def __init__(self, number):
        self.number = number
        self.chains = []

class Chain:
    def __init__(self, letter):
        self.letter = letter
        self.res = []

class Residue:
    def __init__(self, number):
#        self.res_type = res_type
        self.number = int(number)
        self.atoms = []
        self.CA = Atom(0, 0)
        
class Atom:
    def __init__(self, number = 0, name=None, x=0, y=0, z=0, b_factor=0):
        self.number = int(number)
        self.name = name
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.b_factor = float(b_factor)
    
def main():
    global filename
    global opdracht
    print("What PDB-file do you want to be read?")
    filename = sys.stdin.readline()[:-1]
    print("Do you want to extract b-factors? y/n If you answer n, the RMSF of different molecules in the unit cell will be calculated.")
    if sys.stdin.readline()[:-1] == "y":
        print("Do you want the average b-factor per residue (avg), the c-alpha b-factor per residue (ca) or the total average b-factor (tavg)?")
        opdracht = sys.stdin.readline()[:-1]
        #print "yes"
    else:
        opdracht = "RMSF"
#    temperature = sys.argv[3]
#    if temperature == "temp":
#        temp = True
#    print opdracht
    ReadFile()
    if opdracht == "ca":
        print "opdracht == ca"
        Return_CA_Bfactors()
    if opdracht == "avg":
        ReturnAvgBfactor()
    if opdracht == "tavg":
        ReturnTotalAvgB()
    if opdracht == "RMSF":
        SymRMSF()
        #cAlphaRMSF()
        
# reading for RMSF now sees MODEL entries and chains as the same thing, while for MD-pdb's each model is a different time step and the chains are the different molecules
def ReadFile():
    print " in ReadFile"
    global filename
    f = open(filename, 'r')
    l = f.readline()
    mdlnum = 0
    current_mdl = Model(mdlnum)
    models.append(current_mdl)
    while l[0:6] != "MODEL " and l[0:4] != "ATOM":
        l = f.readline()
    while l: # until the eof is reached
        if l[0:3] == "TER" or l[0:6] == "CONECT":
            while l and l[0:6] != "MODEL " and l[0:4] != "ATOM":
                l = f.readline()
            if not l:
                break
        if l[0:6] == "MODEL ":
            # if there are no real models yet
            if current_mdl.number == 0:
                models.pop()
            mdlnum += 1
            current_mdl = Model(mdlnum)
            models.append(current_mdl)
            l = f.readline()
        """if not current_mdl:
            current_mdl = Model(0)
            models.append(current_mdl)"""
        chain_letter = l[21]
        current_chain = Chain(chain_letter)
        current_mdl.chains.append(current_chain)
        #print(l)
        while l[0:3] != "TER" and l[0:6] != "CONECT" and l[21] == chain_letter:
            if l[0:4] == "ATOM":
                res_number = l[22:26]
                current_residue = Residue(res_number)
                current_chain.res.append(current_residue)
                while l[0:3] != "TER" and l[22:26] == res_number:
                    if l[0:4] == "ATOM":
                        if opdracht == "RMSF":
                            current_atom = Atom(number=l[6:11], name=l[12:16], x=l[30:38], y=l[38:45], z=l[46:54])
                        else:
                            current_atom = Atom(name=l[12:16], b_factor=l[60:66])
                        if "CA" in current_atom.name:
                            current_residue.CA = current_atom
                        current_residue.atoms.append(current_atom)
                    l = f.readline()
            else:
                l = f.readline()
    f.close()

def SymRMSF():
    # write to three seperate files the x, y and z-coordinates
    #f = open("x.txt", 'w')
    #x = []
    #y = []
    #z = []
    #print("In SymRMSF")
    #print(models[0].number)
    #for c in models[0].chains:
    #print(c.letter)
    for m in models:
        fx = open("x%s.txt" %m.number, 'w')
        fy = open("y%s.txt" %m.number, 'w')
        fz = open("z%s.txt" %m.number, 'w')
        for c in m.chains:
            #print("model")
            for r in c.res:
                for a in r.atoms:
                    fx.write(str(a.x) + "\n")
                    fy.write(str(a.y) + "\n")
                    fz.write(str(a.z) + "\n")
                #x.append(str(a.x))
                #y.append(str(a.y))
                #z.append(str(a.z))
    fx.close()
    fy.close()
    fz.close()

def cAlphaRMSF():
    for m in models:
        fx = open("x%s.txt" %m.number, 'w')
        fy = open("y%s.txt" %m.number, 'w')
        fz = open("z%s.txt" %m.number, 'w')
        for c in m.chains:
            for r in c.res:
                fx.write(str("%8.3f" % r.CA.x) + "\n")
                fy.write(str("%8.3f" % r.CA.y) + "\n")
                fz.write(str("%8.3f" % r.CA.z) + "\n")
        fx.close()
        fy.close()
        fz.close()
        

#  This method writes all B-factors from the C-alpha's to a file
def Return_CA_Bfactors():
    print "In return ca"
    global filename
    bestand = open(filename[0:-4] +'_CA.txt', 'w')
    if temp:
        bestand.write("# Experimental conditions - Temperature : %s"%T)
    for c in models[0].chains:
        for r in c.res:
            bestand.write(str(r.number) + " " + str(r.CA.b_factor) + "\n")
    bestand.close()

# This method first calculates, and then writes the average B-factor per residue to a file
def ReturnAvgBfactor():
    global filename
    bestand = open(filename[0:-4] + '_AvgB.txt', 'w')
    for c in models[0].chains:
        for r in c.res:
            total = float(0)
            t = int(0)
            for a in r.atoms:
                total += a.b_factor
                t += 1
                AvgB = total/t
            bestand.write(str(r.number) + " " + str(AvgB) + "\n")
    bestand.close()
    
# This method calculates the total average B-factor and writes it to the terminal
def ReturnTotalAvgB():
    for c in models[0].chains:
        for r in c.res:
            total = float(0)
            nres = 0
            for a in r.atoms:
                if "H" not in a.name and "O" not in a.name:
                    total += a.b_factor
                    nres += 1
    tavg = float(total/nres)
    print("The average b-factor (without considering hydrogens) is: %s" %tavg)

if __name__ == "__main__":
    main()
