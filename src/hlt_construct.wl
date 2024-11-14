(* ::Package:: *)

ClearAll[\[Sigma]0, \[Sigma]1, \[Sigma]2, \[Sigma]3, Jsymp];  
\[Sigma]0 = {{1,0},{0,1}}; 
\[Sigma]1 = {{0,1},{1,0}}; 
\[Sigma]2 = {{0,-I},{I,0}}; 
\[Sigma]3 = {{1,0},{0,-1}}; 

(* Function to expand a 2 x 2 matrix in terms of the Pauli matrices and identity with complex coefficients *) 
MatPauli[M_]:=(1/2){Tr[M], Tr[M.\[Sigma]1], Tr[M.\[Sigma]2], Tr[M.\[Sigma]3]}; 

(* Commutator of two matrices *) 
Com[A_, B_] := A.B - B.A; 


(* Construct the real space Hamiltonian for a system specified by the on-site matrix M and the hopping matrix J 
for Nsites unit cells. The optional boolean argument pbc sets the periodic/open boundary conditions. *)
ConstructRealSpaceHlt[J_, M_, Nsites_,pbc_:False] := Module[ {Htemp},
	Htemp = Table[ J\[ConjugateTranspose]KroneckerDelta[i+1,j] + J  KroneckerDelta[i,j+1]+M KroneckerDelta[i,j],{i,Nsites},{j,Nsites}];
	If[pbc, 
		Htemp[[Nsites, 1]]+=J\[ConjugateTranspose]; 
		Htemp[[1,Nsites]]+=J; 
	]; 
	Return[ArrayFlatten[Htemp]];
];




(* Construct the real space Hamiltonian for a system consisting of a set of distinct phases specified by 
hopping and on-site matrices Jlist = {J1,...,Jn} and Mlist = {M1,...,Mn} with identical dimensions. 
The "interface" can consist of arbitrary on-site matrices listed in IMlist connected to the system to the
left/right by the matrices listed in IJList and IJRlist. The number of unit cells for each phase is listed
in Nlist, while dNlist 

If the overall system is periodic, then the number of interfaces is equal to the number of individual 
systems. Thus, if IMlist and Mlist are of equal length, the system is assumed to have periodic boundary
conditions, and open boundary conditions otherwise.  *)
ConstructGenInterfaceHlt[Jlist_, Mlist_, IJLlist_, IJRlist_, IMlist_, Nlist_] := Module[
{Nsys, pbc, Htemp, Hblock, JSize, MSize, IJLSize, IJRSize, IMSize, Sz, JLpad, JRpad, PadArrL, PadArrR, Ran, i},

	Nsys = Dimensions[Mlist][[1]];   (* Number of systems *) 
	pbc = Dimensions[IMlist][[1]] - Nsys; 
	(* pbc = 0 for periodic boundary conditions, -1 otherwise. *)

	(* Compute sizes of all set of matrices. Need to implement consistency checks here at some point *) 
	JSize = Table[Dimensions[Jlist[[i]]][[1]],{i,Nsys}]; 
	MSize = Table[Dimensions[Mlist[[i]]][[1]],{i,Nsys}]; 
	IJLSize = Table[Dimensions[IJLlist[[i]]],{i,Nsys+pbc}]; 
	IJRSize = Table[Dimensions[IJRlist[[i]]],{i,Nsys+pbc}]; 
	IMSize = Table[Dimensions[IMlist[[i]]][[1]],{i,Nsys+pbc}]; 

	(* Sizes of individial blocks, with one block for each single (continuous) phase as well as one for each interface *) 
	Sz = Riffle[MSize*Nlist, IMSize]; 

	(* Construct block Hamiltonians for regions, indexed by m *) 
	Ran[m_]:= 1;;Nlist[[m]]MSize[[m]]; 
	Hblock[m_] := ArrayFlatten[Table[ Jlist[[m]]\[ConjugateTranspose]KroneckerDelta[i+1,j] + Jlist[[m]]  KroneckerDelta[i,j+1]+Mlist[[m]] KroneckerDelta[i,j],{i,Nlist[[m]]},{j,Nlist[[m]]}]][[Ran[m], Ran[m]]];


	(* Start with a template consisting of block matrices with zeros. Will later replace the required blocks. *) 
	Htemp = Table[ConstantArray[0, {Sz[[i]], Sz[[j]]}], {i,1,2Nsys+pbc}, {j,1,2Nsys+pbc}];

	(* Put in the blocks for each phase *) 
	For[i=1, i<= Nsys, i++, 
		Htemp[[2i-1,2i-1]] = Hblock[i]; 
	]; 

	(* Add in the interface matrices. Note that we need to suitably pad the hopping matrices to make their size commensurate with Hblock[]. Generically, Subscript[M, left], Subscript[M, right] and Subscript[M, interface] might all have different sizes.  *) 
	For[i=1, i<= Nsys-1, i++, 
		(* On site matrix at the interface. *) 
		Htemp[[2i,2i]] = IMlist[[i]]; 

		(* Hopping matrix to the left.  *) 
		PadArrL = Join[ConstantArray[0,Nlist[[i]]-1],{1}]; 
		JLpad = ArrayFlatten[{Table[PadArrL[[n]] IJLlist[[i]], {n, 1,Nlist[[i]]}]}]; 		
		Htemp[[2i-1,2i]] = JLpad\[ConjugateTranspose];
		Htemp[[2i,2i-1]] = JLpad; 

		(* Hopping matrix to the right. *) 
		PadArrR = Join[{1}, ConstantArray[0,Nlist[[i+1]]-1]]; 
		JRpad = ArrayFlatten[Table[{PadArrR[[n]] IJRlist[[i]]}, {n, 1,Nlist[[i+1]]}]]; 
		Htemp[[2i,2i+1]] =JRpad\[ConjugateTranspose];
		Htemp[[2i+1,2i]] = JRpad; 
	]; 

	If[pbc == 0, 
		(* Wrapping around the last site. *) 
		Htemp[[2Nsys, 2Nsys]] = IMlist[[Nsys]]; 

		(* Hopping matrix to the left.  *) 
		PadArrL = Join[ConstantArray[0,Nlist[[Nsys]]-1],{1}]; 
		JLpad = ArrayFlatten[{Table[PadArrL[[n]] IJLlist[[Nsys]], {n, 1,Nlist[[Nsys]]}]}]; 		
		Htemp[[2Nsys-1,2Nsys]] = JLpad\[ConjugateTranspose];
		Htemp[[2Nsys,2Nsys-1]] = JLpad; 

		(* Hopping matrix to the right. *) 
		PadArrR = Join[{1}, ConstantArray[0,Nlist[[1]]-1]]; 
		JRpad = ArrayFlatten[Table[{PadArrR[[n]] IJRlist[[Nsys]]}, {n, 1,Nlist[[1]]}]]; 
		Htemp[[2Nsys,1]] = JRpad\[ConjugateTranspose];
		Htemp[[1,2Nsys]] = JRpad; 
	];
	
	(* Finally flatten out the Hamiltonian. *) 
	Htemp = ArrayFlatten[Htemp];
	Return[Htemp]; 
];

(*  Construct the real space Hamiltonian for the "standard" interface, i.e, the interface between 
two phases on a torus, with the only tunable parameter being the coupling between the two surfaces.   *) 
ConstructStdInterfaceHlt[J_, Mlist_, \[Kappa]_, Nlist_] := Module[{Jlist, IJLlist, IJRlist, IMlist},
	Jlist = {J,J}; 
	IMlist = Mlist; 
	IJLlist = {J, J}; 
	IJRlist = {\[Kappa] J, \[Kappa] J};
	Return[ConstructGenInterfaceHlt[Jlist, Mlist, IJLlist, IJRlist, IMlist,Nlist]];
];
