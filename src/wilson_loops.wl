(* ::Package:: *)

(* Some matrix functions required for the computation of Wilson loop. *) 

(* Returns the "num" lowest energy (filled) states for the numerical matrix mat. *) 
OccStates[mat_, num_] := Module[{esys, ord},
	esys = Eigensystem[mat];
	ord = Ordering[esys[[1]], num]; 
	Return[esys[[2, ord]]\[Transpose]];
];

(* Returns the "unitary" part of the SVD for a given matrix mat. *) 
ExtractUMat[mat_] := Module[{U,X,V}, 
	{U,X,V} = SingularValueDecomposition[mat]; 
	Return[U.V\[ConjugateTranspose]]; 
];

(* Stable matrix multiplication routine using the QR decomposition at each step. *) 
MultiplyMat[matlist_] := Module[{len, prod, Q, R, i, temp},
	len = Length[matlist]; 
	prod = IdentityMatrix[Length[matlist[[1]]]]; 
	temp = prod; 
	For[i=len, i >= 1, i--, 
		{Q,R} = QRDecomposition[matlist[[i]].temp]; 
		prod = R . prod; 
		temp = Q\[ConjugateTranspose]; 
	]; 
	Return[temp.prod]; 
];


(* Computation of the Wilson loop. *) 

(* Computes the Wilson loop for a given effectively 1d system. The inputs are a Bloch Hamiltonian H[k] (with k the momentum along the Wilson loop), 
number of occupied bands Noc and number of pts Npts to sample the momentum at. For multidimensional Hamiltonians, pass an anonymous function 
HBloch[#, ky0, kz0, ...]& as the first argument to compute Wilson loop along ky=ky0, kz=kz0 etc.  *) 
ComputeWilson[H_, Noc_, Npts_] := Module[{ulist,Flist,Wilson},
	ulist = Table[OccStates[1.0H[2 Pi*n/Npts], Noc], {n, 0, Npts-1}];     (* Wavefunctions sampled at Npts points in [0,2Pi]. *) 
	Flist = Table[ExtractUMat[ulist[[Mod[n+1, Npts,1]]]\[ConjugateTranspose].ulist[[n]]], {n, Npts}];
	Wilson = MultiplyMat[Reverse[Flist]]; 
	Return[Wilson]; 
];

(* Computes the Wilson loop for the occupied bands along a path specified by a list of k-points with u=u(0) at the end points. *) 
ComputeWilsonPath[H_, Noc_, klist_, u0_] := Module[{ulist,Flist,  W},
	(* If klist does not have any points *) 
	If[Length[klist] <= 1, Return[{{1}}]];
 
	ulist = Table[OccStates[1.0H[k], Noc], {k,klist}]; 
	ulist[[1]] = ulist[[-1]] = u0;
	Flist = Table[ExtractUMat[ulist[[n+1]]\[ConjugateTranspose].ulist[[n]]], {n, Length[klist]-1}];
	W = MultiplyMat[Reverse[Flist]]; 	
	Return[W]
];


(* Computes the (nonabelian) Chern number of a 2D Hamiltonian H for a set of Noc bands (with Noc counted from bottom/top for Noc positive/negative.) *) 
CalculateChern[H_, Noc_, Npts_]:= Module[{Ugrid, Fgrid, Wline}, 
	(* Compute the occupied states for points on a grid covering the Brillouin zone. *) 
	Ugrid = Table[OccStates[1.0H[{(2\[Pi] ind1)/Npts, (2\[Pi] ind2)/Npts}], Noc], {ind1,  Npts}, {ind2,  Npts}]; 
	
	(* Function to compute the Wilson line between two points on the momentum grid.  *) 
	Wline[{m1_, n1_},{m2_, n2_}]:= #/Abs[#]& @  Det[Ugrid[[Mod[m1,Npts,1], Mod[n1, Npts,1]]]\[ConjugateTranspose].Ugrid[[Mod[m2, Npts,1], Mod[n2, Npts,1]]]]; 
	
	(* Compute the flux across each plaquette (indexed by their bottom-left grid point).  *) 
	Fgrid =Table[Log[Wline[{ind1,ind2},{ind1+1,ind2}] Wline[{ind1+1,ind2},{ind1+1,ind2+1}] Wline[{ind1+1,ind2+1},{ind1,ind2+1}] Wline[{ind1,ind2+1},{ind1,ind2}]], {ind1,  Npts}, {ind2,  Npts}];
	
	Return[Chop[1/(2 \[Pi] I) Total[Flatten[Fgrid]]]];
];
