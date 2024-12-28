"""
RSA 4/2/23
https://pynative.com/make-python-class-json-serializable/#:~:text=Use%20toJSON()%20Method%20to%20make%20class%20JSON%20serializable&text=So%20we%20don't%20need,Python%20Object%20to%20JSON%20string.

"""

from maptial.xyz import vectorthree as v3
import pandas as pd

amino_acids = ["ala","arg","asn","asp","cys","gln","glu","gly","his","ile","leu","lys","met","phe","pro","ser","thr","trp","tyr","val"]
class PdbObject(object):
    def __init__(self, pdb_code):
        # PUBLIC INTERFACE        
        self.pdb_code = pdb_code
        self.resolution = -1
        self.exp_method = ""
        self.chains = {}
        self.exc_hetatm = False
        self.lines = []
              
    def __str__(self):
        return f"{self.pdb_code}\t{self.resolution}\t{self.exp_method}"
    
    def lines(self,atoms=False):
        lines = []
        for chain,resdic in self.chains.items():
            for no,res in resdic.items():                
                lines.append(chain + "\t" + str(res))
                for attype,atm in res.atoms.items():
                    if atoms:
                        lines.append(atm)
        return lines
                                
    def add_atoms(self,bio_struc):
        self.bio_struc = bio_struc
        self.resolution = self.bio_struc.header['resolution']        
        self.exp_method = self.bio_struc.header["structure_method"]                
        atomNo = 0
        ridx = 0
        last_bad = ""
        for model in self.bio_struc:
            for chain in model:
                self.chains[chain.id] = {} ####  a chain is a dictionary of residue number to GeoResidue ############################################################
                for residue in chain:
                    aa = residue.get_resname()
                    rid = residue.get_full_id()[3][1]                    
                    resd = PdbResidue(aa,rid,ridx)
                    chain = residue.get_full_id()[2]                    
                    hetatm = residue.get_full_id()[3][0]
                    aah = residue.get_full_id()
                    mutation = aah[3][2] #these are alternative conformations where mutations could occur do not add
                    # in fact remove any that area lready there
                    if mutation.strip() != "":
                        #self.chains[chain][rid] = resd
                        if chain in self.chains:
                            if rid in self.chains[chain]:
                                self.chains[chain].pop(rid)
                    elif mutation.strip() == "":
                        this_bad = chain + "_" + str(rid)
                        if this_bad != last_bad:                                    
                            if rid in self.chains[chain] and self.exc_hetatm:                          
                                # if it is already there something is wrong so lets remove it
                                del self.chains[chain][rid]
                            elif mutation[0] != " " and self.exc_hetatm:                          
                                # if there are mutation insertions they cause us geoemtry problems and we will remove them                        
                                last_bad = chain + "_" + str(rid)
                                if rid in self.chains[chain]: #we know this has alreadybeen checked, it is just a line of code to stop it adding
                                    del self.chains[chain][rid]
                            else:                        
                                if str(hetatm[0:2]) == "H_" and self.exc_hetatm:
                                    if False:
                                        print("HETATM",rid)
                                else:
                                    # only proceed if we explicitly want hetatms                        
                                    ridx = ridx + 1

                                    for atom in residue:
                                        disordered = 'N'
                                        if atom.is_disordered():
                                            disordered = 'Y'
                                            if atom.disordered_has_id("A"):
                                                atom.disordered_select("A")
                                        if atom.get_occupancy() == None:
                                            disordered = 'Y'
                                        elif atom.get_occupancy() < 1:
                                            disordered = 'Y'
                                        atomNo += 1
                                        atom_name = atom.get_name()
                                        occupant = atom.get_full_id()[4][1]
                                        if occupant == ' ':
                                            occupant = 'A'
                                        x = atom.get_vector()[0]
                                        y = atom.get_vector()[1]
                                        z = atom.get_vector()[2]
                                        bfactor = atom.get_bfactor()
                                        occupancy = atom.get_occupancy()
                                        one_atom = PdbAtom(chain,resd,atom_name[0],atom_name,atomNo,disordered,occupancy,bfactor,x,y,z)
                                        self.add_line_string(atomNo, rid, aa, atom_name, chain, occupant, x, y, z, occupancy, bfactor, atom_name[0])
                                        resd.atoms[atom_name] = one_atom
                                    self.chains[chain][rid] = resd
                                    last_bad = ""

    def elementInList(self,element, atomlist):
        for atm in atomlist:
            if element in atm:
                return True, atm
        return False, None
    
    # MAP interface
    def get_coords_key(self,key):
        atm = self.get_atm_key(key)
        if atm == {} or atm == None:
            return "(0,0,0)"
        return f"({atm['x']},{atm['y']},{atm['z']})"

    def get_coords(self,atm):
        if atm == {}:
            return "(0,0,0)"
        return f"({atm['x']},{atm['y']},{atm['z']})"
    
    def get_atm_key(self,key):
        #A:709@C.A
        try:
            sx = key.split("@")
            s01 = sx[0].split(":")
            s23 = sx[1].split(".")
            ch,rid,am,ver = s01[0],s01[1],s23[0],s23[1]            
            for atm in self.lines:                                
                if atm["chain"] == ch:
                    if int(atm["rid"]) == int(rid):
                        if atm["atm"] == am:
                            if atm["version"] == ver: 
                                return atm
        except Exception as e:
            print("Error finding key",key,e)
            return {}


    def dataFrame(self):        
        dicdfs = []        
        for chain,resdic in self.chains.items():
            for no,res in resdic.items():
                for attype,atm in res.atoms.items():     
                    dic={'pdbCode':self.pdb_code,'resolution':self.resolution,
                    'chain':atm.chain,'aa':res.amino_acid,'rid':res.rid,'ridx':res.ridx,                    
                    'atom':atm.atom_name, 'atomNo':atm.atom_no,'element':atm.atom_type,
                    'bfactor':atm.bfactor, 'occupancy':atm.occupancy,
                    'x':atm.x, 'y':atm.y, 'z':atm.z}
                    dicdfs.append(dic)
        return pd.DataFrame.from_dict(dicdfs)
    
    # if we add lines from a cif file I will do something different    
    def add_line_string(self,aid, rid, aa, am, ch, ver, x, y, z, occ, bf, ele):
        #ATOM    435  OE1AGLU A  23      12.418  13.330   0.541  0.58  7.57           O  
        #if "ATOM" in line[:6] or "HETATM" in line[:8]:                        
            #aid = line[6:12].strip()
            #am = line[12:16].strip()
            #ver = line[16:17].strip()
            #aa = line[17:21].strip()
            #ch = line[21:23].strip()
            #rid = line[23:29].strip()
            #x = line[29:38].strip()
            #y = line[38:46].strip()
            #z = line[46:54].strip()
            #occ = line[54:60].strip()
            #bf = line[60:66].strip()
            #ele = line[66:].strip()
        atm = {}
        atm["rid"] = rid
        atm["aid"] = aid
        atm["aa"] = aa
        atm["atm"] = am
        atm["chain"] = ch
        if ver == "":
            ver = "A"
        atm["version"] = ver            
        atm["x"] = x
        atm["y"] = y
        atm["z"] = z
        atm["occupancy"] = occ
        atm["bfactor"] = bf
        atm["element"] = ele        
        self.lines.append(atm)
        #if "HETATM" in line[:8]:
        #    print(line)
        #    print(atm)

    
###################################################################################################################
class PdbResidue:
    """Class for a single residue

    :data amino acid:
    :data atoms:
    """

    def __init__(self, amino_acid,rid,ridx):
        """Initialises a GeoPdb with a biopython structure

        :param biopython_structure: A list of structures that has been created from biopython
        """

        self.atoms = {}
        self.elements = {}
        self.amino_acid = amino_acid
        self.rid = rid
        self.ridx = ridx

    def __str__(self):
       delim = "\t"
       return f"{self.amino_acid}{delim}{self.rid}{delim}{self.ridx}"
    
    def infoResd(self):
       delim = "|"
       return f"{self.amino_acid}{delim}{self.rid}"

###################################################################################################################
class PdbAtom:
    """Class for a single atom

    :data atom_type:
    :data atom_name:
    :data x:
    :data y:
    :data z:
    :data disordered:
    :data occupancy:
    :data bfactor:
    """

    def __init__(self, chain,res,atom_type,atom_name,atom_no,disordered,occupancy,bfactor,x,y,z):
        """Initialises a pdb with a biopython structure

        :param biopython_structure: A list of structures that has been created from biopython
        """

        self.info = {}
        self.chain = chain
        self.res = res
        self.atom_type = atom_type
        self.atom_name = atom_name
        self.atom_no = atom_no
        self.disordered = disordered
        self.occupancy = occupancy
        self.bfactor = bfactor
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        delim = "\t"
        return f"{self.atom_type}{delim}{self.atom_name}{delim}{self.atom_no}{delim}{self.disordered}{delim}{self.occupancy}{delim}{self.bfactor}{delim}({self.x}{self.y}{self.z})"

    def infoAtom(self):
        delim = "|"
        return f"({self.chain}{delim}{self.res.infoResd()}{delim}{self.atom_name}{delim}{self.atom_no})"

    def matchesCriteria(self, criteria,rids=[],dis=-1):        
        if self.disordered == "Y":
            return False
        if criteria == "":
            return True
        crits = criteria.split(",")
        for crit in crits:
            cri = crit.split("|")
            if cri[0].lower() == "aa":
                if self.res.amino_acid.upper() == cri[1].upper() or (cri[1].upper()=="20" and self.res.amino_acid.lower() in amino_acids) :
                    pass
                else:
                    return False
            if cri[0].lower() == "~aa":
                if self.res.amino_acid.upper() == cri[1].upper():
                    return False
                else:
                    pass
            
            elif cri[0].lower() == "rid" and len(rids) > 0:
                if "><" in cri[1]: # this is between the values
                    crits = cri[1].split('>')
                    crit_dis1 = float(crits[0])
                    crit_dis2 = float(crits[1][1:])
                    for rid in rids:
                        if abs(int(rid)-int(self.res.rid)) >= int(crit_dis1) and abs(int(rid)-int(self.res.rid)) <= int(crit_dis2):
                            pass
                        else:
                            return False
                elif "<>" in cri[1]: # these are at the extremes
                    crits = cri[1].split('<')
                    crit_dis1 = float(crits[0])
                    crit_dis2 = float(crits[1][1:])
                    for rid in rids:
                        if abs(int(rid)-int(self.res.rid)) <= int(crit_dis1) or abs(int(rid)-int(self.res.rid)) >= int(crit_dis2):
                            pass
                        else:
                            return False
                elif ">" in cri[1]:
                    crit_dis = cri[1][1:]
                    for rid in rids:
                        if abs(int(rid)-int(self.res.rid)) >= int(crit_dis):
                            pass
                        else:
                            return False
                elif "<" in cri[1]:
                    crit_dis = cri[1][1:]
                    for rid in rids:
                        if abs(int(rid)-int(self.res.rid)) <= int(crit_dis):
                            pass
                        else:
                            return False

            elif cri[0].lower() == "dis" and dis != -1:
                if "><" in cri[1]: #betweem
                    crits = cri[1].split('>')
                    crit_dis1 = float(crits[0])
                    crit_dis2 = float(crits[1][1:])
                    if dis >= crit_dis1 and dis <= crit_dis2:
                        pass
                    else:
                        return False
                elif "<>" in cri[1]: #extremes
                    crits = cri[1].split('<')
                    crit_dis1 = float(crits[0])
                    crit_dis2 = float(crits[1][1:])
                    if dis <= crit_dis1 or dis >= crit_dis2:
                        pass
                    else:
                        return False
                elif "<" in cri[1]:
                    crit_dis = cri[1][1:]
                    if float(dis) <= float(crit_dis):
                        pass
                    else:
                        return False            
                elif ">" in cri[1]:
                    crit_dis = cri[1][1:]
                    if float(dis) >= float(crit_dis):
                        pass
                    else:
                        return False
            
            elif cri[0].lower() == "occ":
                if "=" in cri[1]:
                    crit_occ = cri[1][1:]
                    if float(self.occupancy) == float(crit_occ):
                        pass
                    else:
                        return False            
                elif "<" in cri[1]:
                    crit_occ = cri[1][1:]
                    if float(self.occupancy) <= float(crit_occ):
                        pass
                    else:
                        return False            
                elif ">" in cri[1]:
                    crit_occ = cri[1][1:]
                    if float(self.occupancy) >= float(crit_occ):
                        pass
                    else:
                        return False
                
        return True
