from Axis3 import Axis3
import Constants as const

#This program works only on Autodock PDBQT file format to read structure information

class Atom:
    ''' check pdbqt file format @ http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file
        and PDB file format @ http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
        the format is same as PDB, but with additional info on the partial charge and the atom type
        it is the column position matters, not the spacing!!
        column is numbered starting from 1, not 0
         1~ 6   record name "ATOM  "/"HETATM"
         7~11   serial number, needed for torsion tree branching
        31~38   coordinate X in %8.3f
        39~46   coordinate Y in %8.3f
        47~54   coordinate Z in %8.3f
        71~76   partial charge in %6.3f
        78~79   atom type

        the column in PDBQT is numbered starting from 1 instead of 0, not the same as python
    '''
    def __init__(self,line):
        self.line=line  #the original line is stored, in cases need to print it out in the future
        self.no=int(line[6:11])
        self.coords=Axis3(float(line[31:38]),float(line[38:46]),float(line[46:54]))
        self.charge=float(line[70:76])
        self.atype=line[77:79].strip()  #the atom type may be just using one char followed by a blank space, strip the blank
        #some additional field is need to store the energy info

    def __repr__(self):
        return self.line[:30]+("%8.3f%8.3f%8.3f" %(self.coords.x, self.coords.y, self.coords.z)) + self.line[54:]
        # this output may need to be updated in order to show the energies for dlg

class Branch:
    ''' this branching information is used to define torsions.
        basically, a branch is identified by a rotatable bond (represented by the serial number of the linking atoms)
        a branch have a list of atoms that this rotatable bond will affect
        branching can be nested, i.e. there might be a rotatable bond linked to a flexible part
        for simpler coordinates calculations, work on most nested branch first and then less nested branches.
        rank the branches by their nest level
    '''
    def __init__(self, anchor, link,rank):
        self.torList=[]
        self.anchor=anchor
        self.link=link
        self.rank=rank

    def push(self,serial):
        self.torList.append(serial)

    def __repr__(self):
        return "Branch{%d} (%3d, %3d) " %(self.rank,self.anchor, self.link)+str(self.torList)

class Molecule:
    def __init__(self,filename):
        self.name=filename
        self.torsions=[]
        self.atoms=[]

        stack=[]
        with open(filename,'r') as pdbqt:
            for line in pdbqt:
                if line.startswith('BRANCH'):
                    anchor,link=[int(x) for x in line.split()[1:]]
                    rank=len(stack)
                    stack.append(Branch(anchor,link,rank))
                elif line.startswith('ENDBRANCH'):
                    anchor,link=[int(x) for x in line.split()[1:]]
                    if stack[-1].anchor==anchor and stack[-1].link==link:
                        self.torsions.append(stack.pop())
                    else:
                        # although not specified in anywhere, personally believe this will hold true
                        raise Exception('bad branching nesting!')
                elif line.startswith('ATOM') or line.startswith('HETATM'):
                    self.atoms.append(Atom(line[:-1]))
                    for branch in stack:
                        branch.torList.append(self.atoms[-1].no)

        # reorder the torsion list, put the most nested branch at first
        # this is important in order to make torsions easier
        self.torsions.sort(key=lambda tor:tor.rank,reverse=True)



    def __repr__(self):
        for tor in self.torsions:
            print tor
        return '\n'.join([str(x) for x in self.atoms])

if __name__=='__main__':
    mol=Molecule('Inputs/ind.pdbqt')
    print mol
