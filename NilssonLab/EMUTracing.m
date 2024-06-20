
(* 
 * Algorithms for Elementary Metabolite Unit analysis, developed by Maciek Antoniewicz
 *
 * This implementation allows multiple chemical elements (N, C, ... ) which yields multidimensional MIDs.
 * EMU networks are parameterized by a single flux vector, consisting of uptake, forward,
 * and reverse reactions (concatenated).
 *)


BeginPackage["NilssonLab`EMUTracing`",
	{"NilssonLab`MetabolicNetwork`Common`", "NilssonLab`AtomMap`Common`", "NilssonLab`Utilities`",
	 "NilssonLab`SparseArrays`", "NilssonLab`IsotopeDistributions`"}];


EMUFluxMap;
EMUStoichiometry;
EMUFreeFluxBasis;
EMU;
EMUSize;
EMUOrderedQ;
AtomNumbers;
AtomsToEMU;
MetaboliteEMUs;
EMUReaction;
SubstrateEMU;
ProductEMU;
CondensationQ;
EquivalentEMUs;
FindEMUReactions;
EMUAlgorithm;
EMUSystem;
EMUEquation;
EMULists;
EMUTensors;
InternalEMUs;
InternalTensor;
SubstrateEMUs;
SubstrateTensor;
CondensationEMUs;
CondensationTensor;
Condense;
EMUMarginalMID;
EMUSimulate;
EMUSolution;
InternalMIDs;
SubstrateMIDs;
CondensationMIDs;
EMUDifferentiation;
EMURatioTransform;
EMURatioEquation;
EMURatioSystem;
SimplexConstraints;
EMURatioSimulate;
EMURatioDifferentiation;
EMUGradient;
EMUJacobianMatrix;
OutflowTensor;
ReduceSingleInputs;
MIDResidual;
MIDJacobian;
EMUDiffEquation;


Begin["NilssonLab`EMUTracing`Private`"];



(* Map a flux state to/from an internal vector v used in the EMU system.
   Here n and nu are the number of reactions and uptake fluxes in the original network,
   and revi is the index of reversible reactions *)

EMUFluxMap[M_MetabolicNetwork] :=
	EMUFluxMap[
		Length[Reactions[M]],
		Length[InputMetabolites[M]],
		Flatten[Position[ReversibleQ[Reactions[M]],True]]]	


(* Vector to FluxState. This needs the underlying MetabolicNetwork to recover the release fluxes
   which are not stored in the vector. *)

EMUFluxMap[M_MetabolicNetwork, v_?VectorQ] :=
	EMUFluxMap[M][M,v]

EMUFluxMap[n_, nu_, revi_][M_MetabolicNetwork, v_?VectorQ] :=
	Block[{vrev},
		vrev = Table[0, {n}];
		vrev[[revi]] = Take[v, -Length[revi]];
		FluxState[M, Take[v, nu + {1, n}], vrev]]


(* FluxState to Vector. This defines the order of fluxes in the EMU representation. *)

EMUFluxMap[n_, nu_, revi_][f_FluxState] :=
	Join[BoundaryUptake[f], ForwardFlux[f], ReverseFlux[f][[revi]]]


(*Return a list of flux names*)

EMUFluxMap[M_MetabolicNetwork, "reactions"] :=
	EMUFluxMap[M][M, "reactions"]

EMUFluxMap[n_, nu_, revi_][M_MetabolicNetwork, "reactions"] :=
	Block[{B, nr, S, Sb},
		B = BoundaryNetwork[M];
		(* list of uptake/fwd/rev/release fluxes *)
		Join[
			Take[BoundaryFluxNames[B], nu],
			Table[StringJoin[ReactionID[r], "_f"], {r, Reactions[M]}],
			Table[StringJoin[ReactionID[r], "_r"], {r, Select[Reactions[M], ReversibleQ]}]]]

(*Return forward and reverse index lists*)

EMUFluxMap[n_, nu_, revi_]["index"] := List @@ EMUFluxMap[n, nu, revi][Range[n]]

(*Total number of fluxes*)

NumberOfFluxes[EMUFluxMap[n_, nu_, revi_]] := nu + n + Length[revi]


(* Stoichiometry matrix over the EMU flux vector *)

EMUStoichiometry[M_MetabolicNetwork] :=
	EMUStoichiometry[M, EMUFluxMap[M]]

EMUStoichiometry[M_MetabolicNetwork, EMUFluxMap[n_, nu_, revi_]] :=
	Block[{S, B, nr},
		S = Stoichiometry[M];
		B = Take[Stoichiometry[BoundaryNetwork[M]], Length[Metabolites[M]]];
		ArrayFlatten[{{Take[B,All, nu], S, -S[[All, revi]]}}]]


(* Construct a free flux basis for the EMU flux vector as defined by EMUFluxMap. *)

EMUFreeFluxBasis[M_MetabolicNetwork] :=
	EMUFreeFluxBasis[M, EMUFluxMap[M]]

EMUFreeFluxBasis[M_MetabolicNetwork, EMUFluxMap[n_, nu_, revi_]] :=
	Block[{S, B, nr, free, K},
		S = Stoichiometry[M];
		B = Take[Stoichiometry[BoundaryNetwork[M]], Length[Metabolites[M]]];
		nr = Length[First[B]] - nu;
		{free, K} = FreeFluxBasis[
			SparseArray[ArrayFlatten[{{Take[B, All, -nr], Take[B,All, nu], S, -S[[All, revi]]}}]]];
		{free - nr, Drop[K, nr]}]


(* An EMU is defined by a metabolite m and a list of fragments of atoms,
   denoting an equivalence class fo subsets of atoms induced by molecular symmetry.
   If there is not symmetries, a single fragment is specified. Each fragments is a list
   of integers for each element (atom type) considered. All equivalent fragments must have
   the same dimensions. *)

emuIndexString[al:{{{_Integer...}..}..}] :=
	ListToString[emuIndexString /@ al, "|"]

emuIndexString[al:{{_Integer...}..}] :=
	ListToString[
		Table[
			"(" <> ListToString[
				a /. {i_Integer /; i >= 10 :> FromCharacterCode[97 + i - 10]}, ""] <> ")",
			{a, al}],
		 ","]

emuString[EMU[m_Metabolite, al:{{{_Integer...}..}..}]] :=
	MetaboliteShortName[m] <> "{" <> emuIndexString[al] <> "}"

Format[e:EMU[Metabolite[_String, _String, _], {{{_Integer...}..}..}]] :=
	"<" <> emuString[e] <> ">"

Metabolite[EMU[m_Metabolite, al_]] := m

AtomNumbers[EMU[m_Metabolite, al_]] := al

(* The EMU dimension is the vector of # atoms for each element.*)

Dimensions[EMU[m_Metabolite, {al:{{_Integer...}..}, ___}]] ^:= Length /@ al

Dimensions[Condense[e1_EMU, e2_EMU]] ^:= Dimensions[e1] + Dimensions[e2]

(* EMU size is defined as the number of atoms total, for all elements*)

EMUSize[e_EMU] := Total[Dimensions[e]]

EMUSize[Condense[e1_EMU, e2_EMU]] := EMUSize[e1] + EMUSize[e2]



(* Converts a list of atoms belong to a single metabolite to the corresponding EMU *)

AtomsToEMU[al:{{_Atom...}..}] :=
	EMU[
		First[Cases[al, Atom[m_, _, _] -> m, {2}]],
		{al /. {Atom[m_,e_,x_] -> x}}]


(* Creates a list of whole metabolite EMUs from a atom list. Useful for choosing targets.
   If no elements are specified, all elements present in the atom map are used. *)

MetaboliteEMUs[atomMap_AtomMap] := 
	MetaboliteEMUs[AtomList[atomMap]]

MetaboliteEMUs[atomMap_AtomMap, elements_List] :=
	MetaboliteEMUs[AtomList[atomMap], elements]

MetaboliteEMUs[atomMaps:{_AtomMap..}] := 
	MetaboliteEMUs[AtomList[atomMaps]]

MetaboliteEMUs[atomMaps:{_AtomMap..}, elements_List] :=
	MetaboliteEMUs[AtomList[atomMaps], elements]

MetaboliteEMUs[atoms:{_Atom..}] := 
	MetaboliteEMUs[atoms, Union[AtomElement /@ atoms]]

MetaboliteEMUs[atoms:{_Atom..}, elements_List] := 
	Block[{al, ae},
		al = GatherBy[atoms, Metabolite];
		Flatten[Table[
			ae = Table[Cases[a, Atom[_, e, _]], {e, elements}];
			(* skip empty atom lists *)
			If[Max[Length /@ ae] > 0, AtomsToEMU[ae], {}],
			{a, al}]]]


(* An EMUReaction is a subset of a carbon map for a given reactions, describing formation
   of particular EMUs. The substrate s may be an EMU or a Condense[]. The last argument is
   a multiplier in case of equivalent EMUs *)

SubstrateEMU[EMUReaction[_Integer, s_, _EMU, _]] := s

ProductEMU[EMUReaction[_Integer, _, p_EMU, _]] := p

CondensationQ[e_EMUReaction] := 
	TrueQ[Head[SubstrateEMU[e]] == Condense]


(* The size of any EMU in a reaction must be equal; hence this is well defined *)

Dimensions[r_EMUReaction] ^:= Dimensions[ProductEMU[r]]

EMUSize[r_EMUReaction] := EMUSize[ProductEMU[r]]


(* For multidimensional EMUs, there is no complete ordering of EMU subsystems,
   but we may defined a partial order as the component order of the dimensions.
   This is sufficient to guarantee that the cascading system is solvable
   (all convolutions will exist when required) *)

EMUOrderedQ[e1_EMU, e2_EMU] := EMUOrderedQ[Dimensions[e1], Dimensions[e2]]

EMUOrderedQ[dim1_?VectorQ, dim2_?VectorQ] :=
	And @@ Thread[LessEqual[dim1, dim2]]


(* Given a symmetry list as produced by GetSymmetry, return the equivalence classes of an EMU *)

EquivalentEMUs[EMU[m:Metabolite[metId_, c_,x_], {{al_List}}], symmetry:{Rule[_String,_List]..}] :=
	EMU[m, Sort[List /@ RotationalSymmetry[al, metId /. symmetry /. {_String->{}}]]]


(* MIDs are currently represented by a multidimensional, full array whose depth equals
   the number of elements, e.g. carbon-nitrogen networks have depth 2. For each EMU,
   the associated MID has the dimensions given by EMUSize + 1 (for the unlabeled state). *)


(* A condensation of two EMUs *)

Format[Condense[e1_EMU, e2_EMU]] := 
	"<" <> emuString[e1] <> " x " <> emuString[e2] <> ">"


(* MIDs for condensation EMUs is obtained by convolution, e.g. {a,b} x {c,d,e} = {ac,ad+bc,ae+bd,be}.
   Currently, only binary condensations are supported: higher-order condensation reactions like
   A + B + C --> D
   must be split into intermediate steps, e.g.
   A + B -- >AB   and   AB + C --> D
 *)


(* Note: this is not very efficient; we might try an approach with precomputed sparse tensors
   to handle the summation *)

Condense[A_?ArrayQ, B_?ArrayQ] :=
	Block[{nA, nB, X},
		nA = Partition[Dimensions[A], 1];
		nB = Partition[Dimensions[B], 1];
		X = Table[0, 
			Evaluate[Sequence @@ Partition[Dimensions[A] + Dimensions[B] - 1, 1]]];
		Do[
			Part[X, Sequence @@ (ai+bi-1)] += 
				Part[A, Sequence @@ ai]*Part[B, Sequence @@ bi],
			{ai, Tuples[Table[Range[k], {k, Dimensions[A]}]]},
			{bi, Tuples[Table[Range[k], {k, Dimensions[B]}]]}];
		X]



(* Determine the "marginal" MID for an EMU (a fragment) given an isotopomer distribution of the
  underlying metabolite. Since all equivalent fragments of an EMU must have the same
  isotopomer distribution, we consider only the first fragment. *)

EMUMarginalMID[e:EMU[_Metabolite, fl_List], dist_] :=
	MarginalizeToMID[FragmentIsotopomerDistribution[dist, First[fl]]]


(* EMU decomposition of carbon maps.
   Find all EMU reactions producing a given EMU in a given metabolic reaction
   (assuming that the EMU is indeed produced by the reaction!) 
   We always keep carbon indices sorted, since e,g,
   EMU[x, {1,2}] and EMU[x, {2,1}] are equivalent.
   In most cases this returns a single EMU reaction; however, cleavages like A --> 2B
   produce multiple EMU reactions.
   This only considers the left-to-right direction of the supplied atom map.
   The id argument is returned in the EMUReaction(s) found the reaction but is not
   actually used by this function. *)

FindEMUReactions::noprod = "No such product `1` for reaction `2`";


(* Case of single-fragment target EMU (no symmetry in target).
   This accepts a function equiv[e_EMU] which must return an EMUs of all fragments equivalent to e,
   including e itself. Using equiv = Identity means no symmetry. *)

FindEMUReactions[id_Integer, rmap_AtomMap, elements_List,
				EMU[m_Metabolite, {al:{{_Integer...}..}}], equiv_]:=
	Block[{atomPattern, prodMult, eductList, atomLists, emuList},
		(* Check multiplicity on the product side = # EMU reactions to generate *)
		prodMult = Max[Table[
			Cases[ProductList[rmap],
				{Atom[m, elements[[i]], Alternatives @@ al[[i]]], k_Integer} -> k],
			{i, Length[elements]}]];
		If[!prodMult > 0, Return[{}]];
		(* find list of all educts matching the product emu *)
		eductList = Table[
			Extract[EductList[rmap], 
				Position[ProductList[rmap],
					{Alternatives @@ Table[
						Atom[m, elements[[i]], Alternatives @@ al[[i]]],
						{i, Length[elements]}], k}]],
			{k, prodMult}];
		Table[
			(* for each EMU reaction (product multiplicity) *)
			(* group atoms by metabolite and educt multiplicity *)
			atomLists = GatherBy[educts, {Metabolite[#[[1]]], #[[2]]}&][[All,All,1]];
			(* create reactant EMUs from atoms *)
			emuList = Table[
				EMU[Metabolite[eal[[1]]],
					{Table[
						Sort[Cases[eal, Atom[_, e, x_] -> x]],
						{e, elements}]}],
				{eal, atomLists}];
			(* create EMU reaction from EMUs, applying symmetry function *)
			(* in case of symmetric cleavages all reactions have coefficient 1 *)
			EMUReaction[id,
				If[Length[emuList] == 1,
					(* single EMU \[Equal]> simple EMU reaction *)
					equiv[emuList[[1]]],
					(* > 1 emu ==> condensation reaction *)
					Condense @@ (equiv /@ emuList)],
				EMU[m, {al}], 1],
			{educts, eductList}]]


(* Symmetry case (> 1 equivalent fragments).  Here we first find EMU reactions
   for each equivalent fragment, and then merge duplicates, adjusting coefficients accordingly *)

FindEMUReactions[id_Integer, rmap_AtomMap, elements_List,
				EMU[m_Metabolite, al:{{{_Integer...}..}, __}], equiv_] :=
	Block[{erl, erg, c},
		(* find EMU reactions for each equivalent fragmetnt of the product EMU *)
		erl = Join @@ Table[
			FindEMUReactions[id, rmap, elements, EMU[m, {a}], equiv],
			{a, al}];
		If[Length[erl] == 0, Return[{}]];
		(* replace all product EMUs with the given EMU *)
		erl = Table[
			EMUReaction[First[er], SubstrateEMU[er], EMU[m, al], Last[er]],
			{er, erl}];
		(* group reactions with the same substrate and add up coefficients *)
		{erg, c} = Transpose[Tally[erl]];
		Table[
			EMUReaction[Sequence @@ Most[erg[[k]]], Last[erg[[k]]] * c[[k]] / Length[erl]],
			{k, Length[erg]}]]

(*Find all EMU reactions in a network producing a given EMU.*)

FindEMUReactions[amap:{_AtomMap..}, elements_List, e_EMU, equiv_] :=
	Join @@ Table[
		FindEMUReactions[i, amap[[i]], elements, e, equiv],
		{i, Length[amap]}]

(*Find a list of EMU reactions for each target in el*)

FindEMUReactions[amap:{_AtomMap..}, elements_List, el:{_EMU..}, equiv_] :=
	Table[FindEMUReactions[amap, elements, e, equiv], {e, el}]


(* The EMU decomposition algorithm. Yields a list of EMU reactions sufficient to simulate target EMU t.
   Input EMUs are not included here; they can be found by calling FindEMUReactions on the set
   of EMUs produced by this algorithm. 
   NOTE: reaction indices produced by this algorithm refer to the reaction order used by EMUFluxMap
   *)

EMUAlgorithm::badsize = "Encountered EMU of wrong size! Atom map corrupt?"

(* With rotational symmetry lists, as produced by GetSymmetry
   NOTE: this has not been tested with multiple elements. *)

EMUAlgorithm[M_MetabolicNetwork, amap:{_AtomMap..}, bmap:{_AtomMap..}, elements_List,
			targets:{_EMU..}, symmetry:{Rule[_String,_List]..}] :=
	Block[{equiv},
		(* define equivalece function from symmetry *)
		equiv[EMU[Metabolite[metId_, c_,x_], {{al_List}}]] :=
			EMU[Metabolite[metId,c,x],
				Sort[List /@ RotationalSymmetry[al, metId /. symmetry /. {_String->{}}]]];
		(* applying equiv on an EMU equivalence class returns itself *)
		equiv[e:EMU[_Metabolite, {_List, _List..}]] := e;
		(* run EMU algorithm *)
		EMUAlgorithm[Metabolites[M], Reactions[M], amap, bmap, elements, targets, equiv]]

(*With no symmetry given, the equivalance function is assumed to be Identity*)

EMUAlgorithm[M_MetabolicNetwork, amap:{_AtomMap..}, bmap:{_AtomMap..}, elements_List,
			targets:{_EMU..}]:=
	EMUAlgorithm[Metabolites[M], Reactions[M], amap, bmap, elements, targets, Identity]

(*Here equiv is a function determining EMU equivalance; see FindEMUReactions. *)

EMUAlgorithm[M_MetabolicNetwork, amap:{_AtomMap..}, bmap:{_AtomMap..}, elements_List,
			targets:{_EMU..}, equiv_]:=
	EMUAlgorithm[Metabolites[M], Reactions[M], amap, bmap, elements, targets, equiv]

EMUAlgorithm[ml:{_Metabolite..}, rl:{_Reaction..},
			amap:{_AtomMap..}, bmap:{_AtomMap..}, elements_List,
			targets:{_EMU..}, equiv_]:=
	Block[{erl, maxSize, cp, kp, np, er, allMap},
		maxSize = Max[EMUSize /@ targets];
		erl = {};
		kp = {}; (* known products *)
		cp = targets; (* current products *)
		(* make a single atom map for input, forward, and reverse reactions *)
		allMap = Join[
			Take[bmap, Length[InputMetabolites[ml]]],
			amap,
			Reverse /@ Pick[amap, ReversibleQ /@ rl]];
		(* find all EMU reactions producing the current products *)
		While[cp != {},
			cp = Union[equiv /@ cp]; (* process equivalent EMUs only once *)
			If[Max[EMUSize /@ cp] > maxSize,
				Message[EMUAlgorithm::badsize]; Return[erl]];
			er = Join @@ FindEMUReactions[allMap, elements, cp, equiv];
			erl = Join[erl, er];
			kp = Join[kp, cp];
			(* find next round of products, splitting up condensations *)
			cp = Complement[
				Flatten[(SubstrateEMU /@ er) /. {Condense[e__] -> {e}}],
				kp]];
		erl]


(* Generate systems of EMU distribution equations from a list of EMU reactions.
   This represents a tensor equation A v X + B b Y + C v Z = D v X, where v is the vector
   of absolute fluxes and b is the input boundary exchange.
   With EMUs containing > 1 element, the MIDs are multidimensional.
   We obtain a separatesystem for each distinct dimensionality, since EMU dimensions
   cannot change within a system (without condensing) *)

EMUSystem[emuReactions:{_EMUReaction...}, emap_EMUFluxMap] :=
	EMUSystem[emuReactions, NumberOfFluxes[emap]]

EMUSystem[emuReactions:{_EMUReaction...}, nr_Integer]:=
	Block[{erg, size, br, il, A, sl, B, cl, C},
		(* arrange by increasing fragment size *)
		erg = GatherBy[emuReactions, Dimensions[ProductEMU[#]]&];
		erg = SortBy[erg, Dimensions[ProductEMU[First[#]]]&];
		(* loop over subnetworks *)
		Table[
			(* find substrate reactions for this level *)
			size = EMUSize[ProductEMU[First[er]]];
			(* A tensor, internals *)
			il = Union[ProductEMU /@ er]; (* internal emus *)
			A = makeInternalTensor[il, er, nr];
			(* B tensor, substrates *)
			sl = Cases[Union[SubstrateEMU /@ er], EMU[m_, _] /;
				 Compartment[m] == "Input"];
			B = makeSubstrateTensor[il, sl, er, nr];
			(* C array, condensations *)
			cl = Cases[Union[SubstrateEMU /@ er], _Condense];
			C = makeCondensationTensor[il, cl, er, nr];
			EMUEquation[{il, A}, {sl, B}, {cl, C}],
			{er, erg}]]

makeInternalTensor[internalEMUs:{_EMU..}, er:{_EMUReaction..}, nr_Integer] :=
	Block[{ni, d, rules},
		ni = Length[internalEMUs];
		d = Dispatch[Thread[Rule[internalEMUs, Range[ni]]]];
		(* forward reactions with internals as substrates *)
		rules = Cases[er, EMUReaction[j_, s:Alternatives @@ internalEMUs, p_, c_]
						:> Rule[{p, s, j}, c]] /. d;
		SparseArray[rules, {ni, ni, nr}]]

makeSubstrateTensor[internalEMUs:{_EMU..}, substrateEMUs:{_EMU..},
				 er:{_EMUReaction..}, nr_Integer] :=
	Block[{ni, ns, d, rules},
		ni = Length[internalEMUs];
		ns = Length[substrateEMUs];
		d = Dispatch[Join[
			Thread[Rule[substrateEMUs, Range[ns]]], 
			Thread[Rule[internalEMUs, Range[ni]]]]];
		(* boundary uptake EMU reactions *)
		rules = Cases[er, EMUReaction[j_, s:Alternatives @@ substrateEMUs, p_, c_]
							 :> Rule[{p, s, j}, c]] /. d;
		SparseArray[rules, {ni, ns, nr}]]

makeSubstrateTensor[_, {}, _, _] := {}

makeCondensationTensor[internalEMUs:{_EMU..}, cond:{_Condense..},
					er:{_EMUReaction..}, nr_Integer] :=
	Block[{ni, nc, d, rules},
		ni = Length[internalEMUs];
		nc = Length[cond];
		d = Dispatch[Join[
			Thread[Rule[cond, Range[nc]]], Thread[Rule[internalEMUs, Range[ni]]]]];
		(* reactions with condensations as substrates *)
		rules = Cases[er, EMUReaction[j_, s:Alternatives @@ cond, p_, c_]
						:> Rule[{p, s, j}, c]] /. d;
		SparseArray[rules, {ni, nc, nr}]]

makeCondensationTensor[_, {}, _, _] := {}


(* An EMUEquation consists of the tensors that specify the flux-dependent subsystem equations,
   plus the corresponding EMU lists. *)

Format[EMUEquation[{il_List, A_}, {sl_List, B_}, {cl_List, C_}]] :=
	"<EMUEquation (" <> ListToString[Dimensions[First[il]],","] <> "), " <>
		ToString[Length[il]]<>" internals, " <>
		ToString[Length[sl]] <> " substrates, " <> ToString[Length[cl]] <>
		" condensations >"

EMULists[EMUEquation[{il_List, A_}, {sl_List, B_}, {cl_List, C_}]] := {il, sl, cl}

EMUTensors[EMUEquation[{il_List, A_}, {sl_List, B_}, {cl_List, C_}]] := {A, B, C}

InternalEMUs[eq_EMUEquation] := EMULists[eq][[1]]

InternalTensor[eq_EMUEquation] := EMUTensors[eq][[1]]

SubstrateEMUs[eq_EMUEquation] := EMULists[eq][[2]]

SubstrateTensor[eq_EMUEquation] := EMUTensors[eq][[2]]

CondensationEMUs[eq_EMUEquation] := EMULists[eq][[3]]

CondensationTensor[eq_EMUEquation] := EMUTensors[eq][[3]]

EMUSize[eq_EMUEquation] := EMUSize[First[InternalEMUs[eq]]]

Dimensions[eq_EMUEquation] ^:= Dimensions[First[InternalEMUs[eq]]]

NumberOfFluxes[eq_EMUEquation] := Last[Dimensions[InternalTensor[eq]]]


(* Calculate the diagonal tensors dDdv. These depend on forward, reverse and boundary fluxes,
  returns a triple {Dfwd, Drev, Dsub} *)

OutflowTensor[eq_EMUEquation] := OutflowTensor @@ EMUTensors[eq]

(*Calculate row sums of tensors*)

OutflowTensor[A_SparseArray, B_, C_] :=
	 diagonalTensor[
		Plus @@ Transpose[A] + 
		If[Length[C] > 0, Plus @@ Transpose[C], 0] +
		If[Length[B] > 0, Plus @@ Transpose[B], 0]]


(* This takes an (m x n) sparse matrix and creates a (m x m x n) tensor with the matrix
   on the "diagonal" in the m-m direction.*)

diagonalTensor[d_SparseArray] :=
	Block[{m, n, ar, i, x},
		{m,n} = Dimensions[d];
		ar = Most[ArrayRules[d]];
		{i,x} = {ar[[All,1]],ar[[All,2]]};
		SparseArray[
			Thread[Rule[Transpose[{i[[All,1]],i[[All,1]],i[[All,2]]}],x]],
			{m,m,n}]]


(* This encapsulates a solved EMU system, providing MID values for all EMUs.
   This structure is used for EMU equations in both fluxes and ratios.
   It provides access to MID matrices X, Y and Z. *)

Format[EMUSolution[emus_List, rules_]] :=
	"<EMUSolution, " <> ToString[Length[emus]] <> " subsystems>"

EMULists[EMUSolution[emus_List, rules_]] := emus

InternalMIDs[e:EMUSolution[emus_List, rules_], k_Integer] :=
	emus[[k,1]] /. rules

InternalMIDs[e:EMUSolution[emus_List, rules_]] :=
	emus[[All,1]] /. rules

SubstrateMIDs[e:EMUSolution[emus_List, rules_], k_Integer] :=
	emus[[k,2]] /. rules

CondensationMIDs[e:EMUSolution[emus_List, rules_], k_Integer] :=
	emus[[k,3]] /. rules


(* Applying solution to EMUs *)

EMUSolution[emus_List, rules_][e:(_EMU|_Condense)] := e /. rules

EMUSolution[emus_List, rules_][el:{(_EMU|_Condense)..}] := el /. rules


(* Simulate steady-state MID distributions for an EMU network by solving
   the cascade of linear equations, starting from size 1. At size 1, no convolutions exist;
   moving upward in EMU size, convolutions are always over solutions from previous systems,
   so they evaluate to constants. Substrates are always constants.
   New right-hand side for next EMU size is obtained from convolution, e.g.
    {a,b} x {c,d,e} = {ac,ad+bc,ae+bd,be}
   Here inputs is a list of rules EMU[...] -> {m0,m1, ...} providing the mass isotopomer distribution
   of network substrates, and vr is a list of rules.
   This assumes that the outflow matrix D is nonsingular, i.e. the net flow through every EMU is nonzero.*)

EMUSimulate::nfluxes = "Wrong flux dimensions";

EMUSimulate::neg = "Encountered negative MIDs. Mass balance errors?";

EMUSimulate[eqns:{_EMUEquation..}, inputRules_, v_?VectorQ]:=
	Block[{e, A, B, C, D, X, Y, Z, RHS, LHS},
		If[Length[v] != NumberOfFluxes[First[eqns]],
			Message[EMUSimulate::nfluxes]; $Failed];
		e = inputRules;
		Do[
			(* For each subnetwork: *)
			A = InternalTensor[eq];
			B = SubstrateTensor[eq];
			C = CondensationTensor[eq];
			D = OutflowTensor[eq];
			(* get Y and Z matrices by pattern replace *)
			Y = Flatten /@ (SubstrateEMUs[eq] /. inputRules);
			(* this calculates convolutions *)
			Z = Flatten /@ (CondensationEMUs[eq] /. e);
			LHS = (A-D).v;
			RHS = -If[Length[B] > 0, (B.v).Y, 0]
				- If[Length[C] > 0, (C.v).Z, 0];
			X = LinearSolve[LHS, RHS];
			X = Map[ArrayReshape[#, Dimensions[eq]+1]&, X];
			(* Mass balance errors can cause negative MIDs *)
			If[Chop[Min[X]] < 0,
				Print["Min [X] = ", Chop[Min[X]]];
				Message[EMUSimulate::neg]];
			e = Join[e, Thread[Rule[InternalEMUs[eq], X]]],
			{eq, eqns}];
		EMUSolution[Table[EMULists[eq], {eq,eqns}], Dispatch[e]]]


(* Gradients of MIDs with respect to fluxes.
   These are needed for sensitivity analysi and local optimization (search) methods. *)

Format[EMUGradient[emus_List, rules_]] :=
	"<EMUGradient, " <> ToString[Length[emus]] <> " subsystems>"

EMUJacobianMatrix[EMUGradient[emus_List, rules_]] := Join @@ (Flatten[emus] /. rules)

EMUJacobianMatrix[EMUGradient[emus_List, rules_], targets:{_EMU..}] :=
	Join @@ (targets /. rules)

EMULists[EMUGradient[emus_List, rules_]] := emus

EMUGradient[emus_List, rules_][e_EMU] := e /. rules


(* Analytical differentiation of the EMU equations to find the matrix dx/dv,
   evaluated at the current simulated x(v). This again proceeds sequentially,
   from size 1 EMUs and upward. At the first level, derivatives for EMU substrates are zero. *)

(* This calculate the derivate of all EMU isotopomer distributions w.r.t. flux j,
   at the current point (v, x(v)) given by the rule lists vr, xr. *)

EMUDifferentiation[eqns:{_EMUEquation..}, sol_EMUSolution, v_?VectorQ] :=
	Block[{dr, drk},
		(* all network inputs have derivative zero *)
		dr = {};
		Do[
			drk = EMUDifferentiation[k, eqns[[k]], v, dr, sol];
			dr = Join[dr, drk],
			{k, Length[eqns]}];
		EMUGradient[EMULists[sol], Dispatch[dr]]]

(*For one subnetwork*)

EMUDifferentiation[k_Integer, eq:EMUEquation[{il_List, A_}, {sl_List, B_}, {cl_List, C_}],
				 v_?VectorQ, dr_List, sol_EMUSolution] :=
	Block[{n, midSize, D, Xk, Yk, Zk, Cv, dXdv, dZdv, sr, L, RHS},
		D = OutflowTensor[eq];
		n = Length[v];
		Xk = Flatten /@ InternalMIDs[sol, k];
		Yk = Flatten /@ SubstrateMIDs[sol, k];
		Zk = Flatten /@ CondensationMIDs[sol, k];
		midSize = Length[First[Xk]];
		(* all network inputs have derivative zero *)
		sr = Table[
			Rule[es, SparseArray[{}, Join[Dimensions[eq]+1, {n}]]],
			{es, sl}];
		(* compute condensation derivatives dZ/dv by the chain rule *)
		dZdv = condensationDiff[
					C, cl, sr, dr, Length[Dimensions[eq]], n];
		(* compute internal derivatives dX/dv*)
		L = LinearSolve[(A-D).v]; (* Factorize for repeated solve *)
		Cv = If[Length[C] > 0, C.v, {}];
		dXdv = SparseArray[Table[
			(* solve for derivative w.r.t flux j k *)
			(* dB/dv_j * Y + dC/dv_j * Z + C * dZ/dv_j - dA/dv * X *)
			RHS = SparseArray[-(A[[All,All,j]]-D[[All,All,j]]).Xk +
				If[Length[B] > 0, -B[[All,All,j]].Yk, 0] +
				If[Length[C] > 0, -C[[All,All,j]].Zk - Cv.dZdv[[j]], 0]];
			Chop[SparseArray[L[RHS]]],
			{j, n}]];
		dXdv = Map[
			ArrayReshape[#, Join[Dimensions[eq]+1, {n}]]&,
			Transpose[dXdv, {3,1,2}]];
		If[Length[C] > 0, 
			dZdv = Map[
				ArrayReshape[#, Join[Dimensions[eq]+1, {n}]]&,
				Transpose[dZdv, {3,1,2}]]];
		Join[sr, Thread[Rule[il, dXdv]], Thread[Rule[cl, dZdv]]]]


(* Differeniation for condensations.
   This returns an array of dimensions  fluxes EMUs x Flatten[element1 x element2 ... ]  *)

condensationDiff[C_, cl_, sr_List, dr_List, d_Integer, n_Integer] :=
	Block[{m, x1, x2},
		If[Length[C] > 0,
			m = Length[cl];
			x1 = Map[TransposeToFirst[#, d+1]&, (First /@ cl) /. sr /. dr];
			x2 = Map[TransposeToFirst[#, d+1]&, (Last /@ cl) /. sr /. dr];
			SparseArray[Table[
				Flatten[Condense[x1[[i,j]], x2[[i,j]]]],
				{j, n}, {i, m}]],
			{}]]


(* This uses a list of EMUEquations to generate a list of differential equations in variables
   x[emu, mi][t] where x and t are supplied symbol.
   This can be used by NDSolve to simulate MI time courses. *)

EMUDiffEquation[eqns:{_EMUEquation..}, x_Symbol, t_Symbol, inputRules_, v_?VectorQ, c:{_?VectorQ..}]:=
	Block[{r, m, Am, Bm, Cm, Dm, X, Y, Z, RHS, LHS},
		r = inputRules;
		Table[
			(* For each subnetwork: *)
			Am = InternalTensor[eqns[[k]]];
			Bm = SubstrateTensor[eqns[[k]]];
			Cm = CondensationTensor[eqns[[k]]];
			Dm = OutflowTensor[eqns[[k]]];
			(* variables for internal emus *)
			m = EMUSize[First[InternalEMUs[eqns[[k]]]]];
			X = Table[x[e,j][t], {e, InternalEMUs[eqns[[k]]]}, {j, 0, m}];
			(* get Y and Z matrices by pattern replace *)
			Y = Flatten /@ (SubstrateEMUs[eqns[[k]]] /. inputRules);
			(* this calculates convolutions *)
			Z = Flatten /@ (CondensationEMUs[eqns[[k]]] /. r);
			(* add new variables x[i,j] to replacement list *)
			r = Join[r, Thread[Rule[InternalEMUs[eqns[[k]]], X]]];
			(* Create diff equation *)
			MapThread[Equal, {
				c[[k]]*D[X,t],
				((Am-Dm).v).X + If[Length[Bm]>0, (Bm.v).Y, 0] + If[Length[Cm]>0, (Cm.v).Z, 0]},
				2],
			{k, Length[eqns]}]]


End[];

EndPackage[];
