import copy
from Axis3 import Axis3
from Quaternion import Quaternion
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
        self.axis=None
        self.base=Axis3(0,0,0)

    def push(self,serial):
        self.torList.append(serial)

    def __repr__(self):
        return "Branch{%d} (%3d, %3d) " %(self.rank,self.anchor, self.link)+str(self.torList)


class Molecule:
    def __init__(self,filename):
        self.name=filename
        self.torsions=[]
        self.atoms=[]
        self.coords=[]
        self.about=None
        self.trans=None
        self.quate=None
        self.rbond=None

        stack=[]
        with open(filename,'r') as pdbqt:
            for line in pdbqt:
                if line.startswith('BRANCH'):
                    anchor,link=[int(x) for x in line.split()[1:]]
                    rank=len(stack)
                    stack.append(Branch(anchor,link,rank))
                    self.torsions.append(stack[-1])
                elif line.startswith('ENDBRANCH'):
                    anchor,link=[int(x) for x in line.split()[1:]]
                    if stack[-1].anchor==anchor and stack[-1].link==link:
                        if self.atoms[anchor-1].no==anchor and self.atoms[link-1].no==link:
                            stack[-1].axis=self.atoms[anchor-1].coords - self.atoms[link-1].coords
                            stack[-1].axis.normalize()
                            stack[-1].base=self.atoms[link-1].coords
                        else:
                            raise Exception('bad atom serial numbering')
                        stack.pop()
                    else:
                        # although not specified in anywhere, personally believe this will hold true
                        raise Exception('bad branching nesting!')
                elif line.startswith('ATOM') or line.startswith('HETATM'):
                    self.atoms.append(Atom(line[:-1]))
                    for branch in stack:
                        if self.atoms[-1].no!=branch.anchor and self.atoms[-1].no!=branch.link:
                            branch.torList.append(self.atoms[-1].no)
        self.coords=[x.coords for x in self.atoms]

        # reorder the torsion list, put the branches with less movable atoms at first
        # this is important in order to make torsions easier
        self.torsions.sort(key=lambda tor:len(tor.torList))

    def setAbout(self,about):
        self.about=about
        for atom in self.atoms:
            atom.coords -= about
        for branch in self.torsions:
            branch.base -= about
        self.coords=[x.coords for x in self.atoms]

    def __repr__(self):
        for tor in self.torsions:
            print tor
        return '\n'.join([str(x) for x in self.atoms])

    def transform(self,move,rot):
        ''' Move the molecule as a whole according to the Axis3 vector for translation and Quaternion rot for rotation'''
        rotParam=rot.getRot()
        for atom in self.atoms:
            atom.coords.transform(rotParam,move)

    def twist(self, angles):
        ''' apply the torsion angles.
            alway do twist before transform anything
        '''
        if len(angles)!=len(self.torsions):
            raise Exception('bad torsion parameter size')
        for angle,torsion in zip(angles,self.torsions):
            quat=Quaternion()
            quat.set_angle_axis(angle,torsion.axis)
            torParam=quat.getRot()
            for atom in [x for x in self.atoms if x.no in torsion.torList]:
                #print torsion.torList
                atom.coords = atom.coords - torsion.base
                atom.coords.transform(torParam,torsion.base)


    def resetCoords(self):
        for atom, coord in zip(self.atoms, self.coords):
            atom.coords=copy.deepcopy(coord)

if __name__=='__main__':
    #mol=Molecule('Inputs/ind.pdbqt')
    mol=Molecule('test/1pgp_lig.pdbqt')
    #print mol
    mol.setAbout(Axis3(22.894,28.598,40.259))
    
    mol.resetCoords()
    mol.twist([const.DEG2RAD*x for x in [92.30, -123.03, -29.07, 23.56, -166.38, 136.34, 37.65, 44.53, 25.79, 118.89, -87.74]])
    mol.transform(Axis3(29.905504,27.014664,43.309807),Quaternion(0.241375, 0.374112,-0.361016,-0.819418))
    print mol
    mol.resetCoords()
    mol.twist([x*const.DEG2RAD for x in [37.65,44.53,25.79,118.89,-87.74,106.41,82.69,-46.02,84.62,-58.47,153.66]])
    mol.transform(Axis3(22.894,28.598,40.259),Quaternion(1, 0,0,0))
    #mol.setAbout(Axis3(0.368900,-0.214800,-4.986500))
    #mol.transform(Axis3(2.056477,5.846611,-7.245407),Quaternion(0.53221, 0.379383,0.612442,0.444674))
    #mol.transform(Axis3(2.742728,5.886342,-7.713194),Quaternion(0.636998, 0.470398,-0.503061,-0.346249))
    print mol
