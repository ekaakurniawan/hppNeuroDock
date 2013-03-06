#This is an alternative implementation of Autodock.

The goal is to make the use of FPGAs to accelerate the docking.

We are also doing a python partial implementation of Autodock, specifically from the designing variables to the energy score part, that is:  
*Input:*  
a structure of the ligand S,  
a vector of designin variables V,  
and a scoring function F,  
*Output:*  
S'=transform(S,V)  
and F(S')  
