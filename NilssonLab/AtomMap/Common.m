
(* This package provides atom-to-atom maps for metabolic networks. *)

BeginPackage["NilssonLab`AtomMap`Common`",
	{"NilssonLab`MetabolicNetwork`Common`","NilssonLab`Utilities`"}];

Atom;
AtomElement;
AtomNumber;
CheckAtomMap;
AtomList;
AtomMap;
EductList;
EductAtoms;
EductIndices;
EductArity;
ProductList;
ProductAtoms;
ProductArity;
ProductIndices;
GetStochiometry;
SubsetAtomMap;
SubsetByEducts;
SubsetByProducts;
AtomToAtomAdjacencyMatrix;
RemoveUnreachable;
ShortestAtomPath;
ReduceToBinaryReactions;
BoundaryAtomMaps;
AtomHashMap;
RenumberAtoms;
FindFreeAtomFluxes;
RotationalSymmetry;
AtomRotations;


Begin["NilssonLab`AtomMap`Common`Private`"];

(* An Atom is represented by a Metabolite, an element and an atom number, which can be any positive integer. *)

Format[Atom[m_Metabolite, e_String, nr_Integer]] :=
	"<" <> MetaboliteShortName[m] <> "(" <> e <> ToString[nr] <> ")>"

Metabolite[Atom[m_Metabolite, e_String, nr_Integer]] ^:= m

AtomElement[Atom[m_Metabolite, e_String, nr_Integer]] ^:= e

AtomNumber[Atom[m_Metabolite, e_String, nr_Integer]] := nr


(* For a given Reaction, an AtomMap is an undirected one-to-one mapping
  between atoms in the reaction's substrates and its products. If a metabolites occurs
  with higher stoichiometry than 1, multiple copies are present in the AtomMap, and therefore
  each educt and product has an index that ranges from 1 to the stoichiometry coefficient
  of the metabolite. All stoichiometry must be integers. 
  Note: while the mapping is undirected, by convention the educt-product order always agrees
  with the reference direction for the Reaction (left-to-right). This simplifies lookup. *)


Format[AtomMap[educts:{{_Atom, _Integer}...}, products:{{_Atom, _Integer}...}]] :=
	"<AtomMap, " <> ToString[Length[educts]] <> " transitions>"

EductList[AtomMap[educts:{{_Atom, _Integer}...}, products:{{_Atom, _Integer}...}]] := educts

EductAtoms[r_AtomMap] := EductList[r][[All,1]]

EductIndices[r_AtomMap] := EductList[r][[All,2]]

ProductList[AtomMap[educts:{{_Atom, _Integer}...}, products:{{_Atom, _Integer}...}]] := products

ProductAtoms[r_AtomMap] := ProductList[r][[All,1]]

ProductIndices[r_AtomMap] :=ProductList[r][[All,2]]


(* Reverses the order of a reaction map, exchanging educts and products *)

Reverse[AtomMap[educts:{{_Atom, _Integer}...}, products:{{_Atom, _Integer}...}]] ^:= 
	Block[{ix},
		(* sort by product, then by atom *)
		ix = Ordering[products /. {{Atom[m_Metabolite,e_, a_],i_} ->{m,i,e,a}}];
		AtomMap[products[[ix]], educts[[ix]]]]


(* Generates an ordered list of all Atoms in a map *)

AtomList[rmap_AtomMap] :=
	Union[Flatten[{EductAtoms[rmap], ProductAtoms[rmap]}]]

AtomList[amap:{_AtomMap..}] :=
	Union[Flatten[Table[
		{EductAtoms[rmap], ProductAtoms[rmap]},
		{rmap, amap}]]]


(* To find arity (number of distinct metabolites) we must count multiplicities,
   e.g. 2A + B --> C is ternary on the educt side*)

EductArity[rmap_AtomMap] :=
	Total[Max /@ GatherBy[
		Thread[{Metabolite/@EductAtoms[rmap], EductIndices[rmap]}],First][[All,All,2]]]

ProductArity[rmap_AtomMap] :=
	Total[Max /@ GatherBy[
		Thread[{Metabolite/@ProductAtoms[rmap], ProductIndices[rmap]}],First][[All,All,2]]]


(* Return a pair of lists {metabolites, coefficients} describing the stoichiometry
   of the entries in this map*)

GetStochiometry[rmap_AtomMap] :=
	Reverse[Transpose[Join[
		First /@ Sort /@ GatherBy[
			Thread[{-EductIndices[rmap], Metabolite/@EductAtoms[rmap]}], Last],
		Last /@ Sort /@ GatherBy[
			Thread[{ProductIndices[rmap], Metabolite/@ProductAtoms[rmap]}], Last]]]]


(* Retain only the transitions in which both educt and product are members of the provided atom list *)

SubsetAtomMap[rmap:AtomMap[educts_, products_], atoms:{_Atom..}] :=
	Block[{t},
		t = Cases[Thread[{educts, products}],
			{{Alternatives @@ atoms, _}, {Alternatives @@ atoms, _}}];
		AtomMap[t[[All,1]], t[[All,2]]]]


(* Retain only the transitions in which educts are members of the provided educt list *)

SubsetByEducts[AtomMap[educts_List, products_List], el:{{_Atom, _Integer}..}] :=
	Block[{i},
		i = Flatten[Position[educts, Alternatives @@ el]];
		AtomMap[educts[[i]], products[[i]]]]


(* Retain only the transitions in which educts are members of the provided educt list *)

SubsetByProducts[AtomMap[educts_List, products_List], pl:{{_Atom, _Integer}..}] :=
	Block[{i},
		i = Flatten[Position[products, Alternatives @@ pl]];
		AtomMap[educts[[i]], products[[i]]]]


(* Verify consistency of an atom map: each reaction that has an atom map entry must provide
  a transition for every atom of every metabolite it interacts with. For reactions with 
  metabolite stoichiometries > 1 all "instances" are checked. *)

CheckAtomMap[amap:{_AtomMap..}] :=
	Block[{ar, em, pm, diff},
		ar = GatherBy[AtomList[amap], Metabolite];
		ar = Dispatch[Thread[Rule[Metabolite /@ First /@ ar, ar]]];
		diff = Table[
			{k, 
			 If[Length[EductList[amap[[k]]]] > 0,
				em = Union[
					EductList[amap[[k]]] /. {Atom[m_Metabolite, __], i_Integer} -> {m, i}];
				Join @@ Table[
					DeleteCases[
						Thread[{e[[1]] /. ar, e[[2]]}],
						Alternatives @@ EductList[amap[[k]]]],
					{e, em}],
				{}],
			 If[Length[ProductList[amap[[k]]]] > 0,
				pm = Union[
					ProductList[amap[[k]]] /. {Atom[m_Metabolite, __], i_Integer} -> {m, i}];
				Join @@ Table[
					DeleteCases[
						Thread[{p[[1]] /. ar, p[[2]]}],
						Alternatives @@ ProductList[amap[[k]]]],
					{p, pm}],
				{}]},
			{k, Length[amap]}];
		Cases[diff, Except[{_, {}, {}}], {1}]]



(* Generate an atom map for boundary exchange corresponding to BoundaryNetwork[M],
   resulting in one entry for each boundary flux. Non-mapped boundary reactions will have
   empty atom lists. *)

BoundaryAtomMaps[M_MetabolicNetwork, atomMap:{_AtomMap..}] :=
	BoundaryAtomMaps[M, AtomList[atomMap]]

BoundaryAtomMaps[M_MetabolicNetwork, atoms:{_Atom..}] :=
	Flatten[Join[
		Table[
			If[MemberQ[{"Uptake", "Free"}, ExchangeStatus[m]],
			AtomMap[
					Cases[atoms, Atom[m, e_, a_] -> 
						{Atom[Metabolite[MetaboliteID[m], "Input", "Boundary"], e, a], 1}],
					Cases[atoms, Atom[m, e_, a_] -> {Atom[m, e, a], 1}]],
				{}],
			{m, ExchangeMetabolites[M]}],
		Table[
			If[MemberQ[{"Release", "Free"}, ExchangeStatus[m]],
				AtomMap[
					Cases[atoms, Atom[m, e_, a_] -> {Atom[m, e, a], 1}],
					Cases[atoms, Atom[m, e_, a_] -> 
						{Atom[Metabolite[MetaboliteID[m], "Output", "Boundary"], e, a], 1}]],
				{}],
			{m, ExchangeMetabolites[M]}]]]



(* Generate a sparse adjacency matrix A where Subscript[A, ij]=1 iff atom i can be converted
   to atom j by some reaction(s), representing a directed graph. Reversible reactions yield
   bidirectional edges. *)

AtomToAtomAdjacencyMatrix[M_MetabolicNetwork, amap:{_AtomMap..}] :=
	AtomToAtomAdjacencyMatrix[M, amap, AtomList[amap]]

(* Over a given list of atoms (also saves some calculation time if atoms is known) *)

AtomToAtomAdjacencyMatrix[M_MetabolicNetwork, amap:{_AtomMap..}, atoms:{_Atom..}] :=
	AtomToAtomAdjacencyMatrix[amap, atoms, ReversibleQ /@ Reactions[M]]

(* Weighting edges by fluxes *)

AtomToAtomAdjacencyMatrix[M_MetabolicNetwork, amap:{_AtomMap..}, flux_FluxState] :=
	AtomToAtomAdjacencyMatrix[amap, AtomList[amap], ReversibleQ /@ Reactions[M],
		 {ForwardFlux[flux], ReverseFlux[flux]}]

AtomToAtomAdjacencyMatrix[amap:{_AtomMap..}, atoms:{_Atom..}, rev:{(True|False)..}] :=
	AtomToAtomAdjacencyMatrix[amap, atoms, rev, {Table[1, {Length[rev]}], Table[1, {Length[rev]}]}]

AtomToAtomAdjacencyMatrix[amap:{_AtomMap..}, atoms:{_Atom..}, rev:{(True|False)..}, {f_, r_}] :=
	Block[{na, rules},
		na = Length[atoms];
		rules = Flatten[Join[Table[
			{Thread[Rule[Thread[{EductAtoms[amap[[i]]], ProductAtoms[amap[[i]]]}], f[[i]]]],
			 If[rev[[i]],
				Thread[Rule[Thread[{ProductAtoms[amap[[i]]], EductAtoms[amap[[i]]]}], r[[i]]]],
				{}]},
			{i, Length[rev]}]]];
		SparseArray[
			rules /. Dispatch[Thread[Rule[atoms, Range[na]]]],
			{na, na}]]


(* Find the shortest path between two atoms, if one exists; otherwise {} is returned. *)

ShortestAtomPath[amap:{_AtomMap..}, a1_Atom, a2_Atom] :=
	ShortestAtomPath[
		{AtomList[amap], AtomToAtomAdjacencyMatrix[amap]},
		a1, a2]

ShortestAtomPath[{atomList:{_Atom..}, adj_SparseArray}, a1_Atom, a2_Atom] :=
	FindShortestPath[AdjacencyGraph[atomList, adj], a1, a2]


(* Remove all atoms and corresponding atom-atom transitions that are not reachable
   from an Uptake metabolite. *)

RemoveUnreachable[M_MetabolicNetwork, amap:{_AtomMap..}] :=
	RemoveUnreachable[M, amap, 
		Cases[Metabolites[M], Metabolite[_, _, "Uptake"|"Free"]]]

RemoveUnreachable[M_MetabolicNetwork, amap:{_AtomMap..}, inmet:{_Metabolite..}] :=
	Block[{adj, atoms, reachable, rv, newmap},
		atoms = AtomList[amap];
		adj = AtomToAtomAdjacencyMatrix[M, amap, atoms];
		(* make binary vector over all atoms, 1 = reached *)
		rv = atoms /. {Atom[Alternatives @@ inmet, _, _] -> 1, Atom[_,_,_] -> 0};
		rv = FixedPoint[UnitStep[(# + #.adj) - 1]&, rv];
		reachable = Pick[atoms, rv, 1];
		(* subset map (this is not very efficient!)  *)
		Table[SubsetAtomMap[rmap, reachable], {rmap, amap}]]


(* Reduce arbitrary metabolic network so that there are at most two reactants and two products
   in the associated atom map. The new metabolites introduced have compartment = "Virtual".
   This conversion is required for export to the FTBL format. *)

ReduceToBinaryReactions::missing = "No reaction `1` in network'";

(*This reduces all reactions with arity > 3 at input or output to binary reactions*)

ReduceToBinaryReactions[M_MetabolicNetwork, amap:{_AtomMap..}] :=
	Module[{ne, np, red, Mred, amapRed},
		(* Find reactions with > 2 educts or > 2 products in the atom map *)
		ne = EductArity /@ amap;
		np = ProductArity /@ amap;
		(* these reactions will be removed and replaced with new *)
		red = Union[Pick[Reactions[M], ne, x_/; x > 2],
					Pick[Reactions[M], np, x_/; x > 2]];
		{Mred, amapRed} = {M, amap};
		Do[
			{Mred, amapRed} = ReduceToBinaryReactions[Mred, amapRed, r],
			{r, red}];
		{Mred, amapRed}]


(*This reduces one given reaction to binary reactions*)
(*NOTE: there may be a special case unaccounted for when # products = 1. Check this.*)

ReduceToBinaryReactions[
		M:MetabolicNetwork[rl_List, ml_List, S_SparseArray],
		amap:{_AtomMap..}, r_Reaction] :=
	Module[{ri, rmap, educts, nonAtomSubstrates, products, nonAtomProducts,
			met, coeff, im, Mred, amapRed, newrmap},
		ri = Position[rl, r];
		If[Length[ri] != 1,
			Message[ReduceToBinaryReaction::missing, ReactionID[r]];
			Return[$Failed]];
		ri = ri[[1,1]];
		(* find metabolites present in the reaction atom map *)
		rmap = amap[[ri]];
		(* educts, with indexes for multiplicities *)
		educts = Union[Thread[{Metabolite /@ EductAtoms[rmap], EductIndices[rmap]}]];
		nonAtomSubstrates = DeleteCases[
			ReactionSubstrates[M, r], {Alternatives @@ educts[[All,1]], _}];
		products = Union[Thread[{Metabolite /@ ProductAtoms[rmap], ProductIndices[rmap]}]];
		nonAtomProducts = DeleteCases[
			ReactionProducts[M, r], {Alternatives @@ products[[All,1]],_}];
		(* Remove original reaction *)
		Mred = DeleteReactions[M, {r}];
		amapRed = Delete[amap, ri];
		(* first condensing step must handle nonAtomSubstrates and stochiometry > 1 *)
		If[Length[educts] >= 2,
			{met, coeff} = getStoichiometry[{{educts[[1,1]],1}, {educts[[2,1]],1}}];
			im = mergedMetabolite[educts[[1;;2,1]]];
			Mred = AddReactions[Mred,
				condReaction[r, 1],
				Join[met, nonAtomSubstrates[[All,1]], {im}],
				Join[-coeff, -nonAtomSubstrates[[All,2]], {1}]];
			amapRed = Append[amapRed, 
				condReactionMap[rmap, educts[[1]], educts[[2]], True]];
			newrmap = replaceEducts[rmap, educts[[1]], educts[[2]], im],
			(* else single educt: create identity map but remove nonAtomSubstrates *)
			im = Metabolite[MetaboliteID[educts[[1,1]]], "Virtual", "Internal"];
			Mred = AddReactions[Mred,
				condReaction[r, 1],
				Join[{educts[[1,1]]}, nonAtomSubstrates[[All,1]], {im}],
				Join[{-1}, -nonAtomSubstrates[[All,2]], {1}]];
			amapRed = Append[amapRed, 
				AtomMap[EductList[rmap], EductList[rmap] /. {educts[[1,1]] -> im}]];
			newrmap = rmap];
		(* do intermediate condensing and cleaving steps *)
		{Mred, amapRed, newrmap} = reduceSubstrates[Mred, amapRed, 
			newrmap,
			Join[{{im, 1}}, Drop[educts, Min[2, Length[educts]]]], products, r, 2];
		(* final cleaving step, return result *)
		{met, coeff} = getStoichiometry[{{products[[-2,1]],1}, {products[[-1,1]],1}}];
		(* substrate for last cleavage is the last two products merged,
		   or the last condensed metabolite is # products <= 2 *)
		im = If[Length[products] > 2,
			mergedMetabolite[Take[products[[All,1]],-2]],
			Metabolite[First[ProductAtoms[Last[amapRed]]]]];
		{AddReactions[Mred, cleavReaction[r, Length[products]-1],
			Join[{im}, met, nonAtomProducts[[All,1]]],
			Join[{-1}, coeff, nonAtomProducts[[All,2]]]],
		 Append[amapRed, AtomMap[EductList[newrmap],
			 If[Length[met] == 1,
				ProductList[newrmap] /. {{a_,n_/; n>2} -> {a,2}},
				ProductList[newrmap] /. {{a_,n_/; n>1} -> {a,1}}]]]}]


(*Recursive step to reduce substrates*)

reduceSubstrates[M_MetabolicNetwork, amap:{_AtomMap..}, rmap_AtomMap,
				educts:{{_Metabolite, _Integer}..}, products:{{_Metabolite, _Integer}..},
				r_Reaction, counter_Integer] :=
	Module[{im},
		If[Length[educts] >= 2,
			im = mergedMetabolite[educts[[1;;2,1]]];
			(* recurse *)
			reduceSubstrates[
				(* generate new condensation reaction *)
				AddReactions[M, condReaction[r, counter],
					Append[educts[[1;;2,1]], im],
					{-1,-1, 1}],
				(* generated map for current step *)
				Append[amap, 
					condReactionMap[rmap, educts[[1]], educts[[2]], False]],
				(* AtomMap remaining for next step *)
				replaceEducts[rmap, educts[[1]], educts[[2]], im],
				(* remaining educts, same products *)
				Join[{{im, 1}}, Drop[educts, 2]], products,
				r, counter+1],
		(* else condensation is complete, next do cleaving steps *)
		reduceProducts[M, amap, rmap, educts[[1,1]], products, r, 1]]]

(* Make a reaction map for a condensation that merges the specified two educts {e1, c1} and {e2, c2}
   NOTE: the coefficients c1, c2 are always 1 except possibly at the first recursion step where c2 = 2
   may occur.*)

condReactionMap[rmap_AtomMap, {e1_Metabolite, c1_Integer}, {e2_Metabolite, c2_Integer},
				step1:True|False] :=
	Block[{eductSub, ne, im},
		(* new carbon map entry is a subset on the educt side *)
		eductSub = Cases[EductList[rmap], {Atom[e1, _, _], c1} | {Atom[e2, _, _], c2}];
		If[!step1, eductSub[[All,2]] = 1];
		ne = Length[eductSub];
		im = mergedMetabolite[{e1, e2}];
		AtomMap[
			eductSub,
			Thread[Thread[{Atom[im, AtomElement /@ eductSub[[All,1]], Range[ne]], 1}]]]]

(* Recursive reduction on the product side. *)

reduceProducts[M_MetabolicNetwork, amap:{_AtomMap..}, rmap_AtomMap,
				educt_Metabolite, products_, r_Reaction, counter_Integer] :=
	Module[{im,newmap},
		newr = cleavReaction[r, counter];
		If[Length[products] > 2,
			(* generate new cleavage reaction *)
			im = mergedMetabolite[Rest[products[[All,1]]]];
			newmap = cleaveReactionMap[rmap, products, r, counter];
			(* recurse *)
			reduceProducts[
				AddReactions[M, cleavReaction[r, counter],
					{educt, products[[1,1]], im},
					{-1, 1, 1}], 
				(* generated map for current step *)
				Append[amap, newmap],
				(* for next step, replace educts with intermediate and remove first product *)
				AtomMap[
					Cases[ProductList[newmap], {Atom[im, _, _], 1}],
					Extract[ProductList[rmap], 
						Position[ProductList[newmap], {Atom[im, _, _], 1}]]],
				im, Rest[products],
				r, counter+1],
		(* else recursion done *)
		(* we return rmap here as it is needed in the final step *)
		{M, amap, rmap}]]


(*Generate a reaction map for a cleavage reaction.*)

cleaveReactionMap[rmap_AtomMap, products_, r_Reaction, counter_Integer] :=
	replaceOtherProducts[rmap,
		First[products], mergedMetabolite[Rest[products[[All,1]]]]]



(* Get stoichiometry for a pair of educts, reducing the case xA + yA to (x+y)A *)

getStoichiometry[{{a_Metabolite, 1}, {b_Metabolite, 1}}] := {{a,b}, {1,1}}

getStoichiometry[{{a_Metabolite, 1}, {a_Metabolite, 1}}] := {{a},{2}}

condReaction[r_Reaction, counter_Integer] :=
	Reaction[ReactionID[r] <> "_cond_" <> ToString[counter], ReversibleQ[r]]

cleavReaction[r_Reaction, counter_Integer] :=
	Reaction[ReactionID[r] <> "_cleav_" <> ToString[counter], ReversibleQ[r]]


(*Produce an A_B_C ..  string for an intermediate "merged" metabolite*)

mergedMetabolite[products:{_Metabolite..}] :=
	Metabolite[mergedMetaboliteID[MetaboliteID /@ products], "Virtual", "Internal"]

mergedMetaboliteID[productIDs:{_String..}] :=
	If[Length[productIDs] > 1,
		productIDs[[1]] <> "_" <> mergedMetaboliteID[Rest[productIDs]],
		productIDs[[1]]]


(* Replace components of atom maps (several special cases, should have a more general function for this) *)

replaceEducts[AtomMap[educts_List, products_List], 
		{e1_Metabolite, c1_Integer}, {e2_Metabolite, c2_Integer}, im_Metabolite] :=
	AtomMap[ 
		educts /. Module[{k = 1}, 
			{Atom[e1, elem_, _], c1} | {Atom[e2, elem_, _], c2} :> {Atom[im, elem, k++], 1}],
		 products]


(*TODO: this does not replace elements properly!*)

replaceOtherProducts[
		AtomMap[educts_List, products_List], product_, im_Metabolite] :=
	AtomMap[educts, 
		Replace[products,
			Module[{k = 1}, 
				{Except[{Atom[product[[1]], _, _], product[[2]]}] :> {Atom[im, "?", k++], 1},
				 {Atom[product[[1]], e_, a_], product[[2]]} -> {Atom[product[[1]], e, a], 1}}],
			1]]


(* Renumber all atoms in an atom map so that every metabolite uses consecutive number 1 ... n
   Some software requires contiguous atom numbers, and InChI-based atom numbering of molecules
   may include atom numbers that are not traced. *)

RenumberAtoms[amap:{_AtomMap..}] :=
	Block[{atoms, rules},
		atoms = GatherBy[AtomList[amap], Metabolite];
		rules = Flatten[Table[
			Thread[Rule[al, Thread[Atom[Metabolite /@ al, AtomElement /@ al, Range[Length[al]]]]]],
			{al, atoms}]];
		amap /. Dispatch[rules]]


(* Finds free fluxes in a metabolic network, considering only the reactions that are present
   in an atom map. See also FindFreeFluxes in package MetabolicNetwork *)

FindFreeAtomFluxes[M_MetabolicNetwork, amap:{_AtomMap..},
			B_BoundaryNetwork, bmap:{_AtomMap..}] :=
	Block[{ncmet,Mdel,Bdel,ri,bi},
		ncmet = Complement[Metabolites[M], Union[Metabolite /@ AtomList[amap]]];
		Mdel = DeleteMetabolites[M, ncmet];
		Bdel = BoundaryNetwork[Mdel];
		{ri,bi} = FindFreeFluxes[Mdel, Bdel];
		{ReindexList[Reactions[Mdel][[ri]], Reactions[M]],
		 ReindexList[BoundaryFluxNames[Bdel][[bi]], BoundaryFluxNames[B]]}]


(* Generate n-degree rotational symmetry given a list of rotations
   {a_1 -> a_2 -> ... -> a_n, b_1 -> b_2 -> ... -> b_n,  ...}
   and the list of all atoms in the molecule. The equivalent EMUs are, for each subset
   of the rotating atoms, collection of atoms after rotating 1,2,3 ... n times.
   For example, with a rotation list {1,2,3} and the atom set {1,2}, the equivalence classes are {2,3}, {1,3}, {1,2}
   *)

RotationalSymmetry[atomList:{_Integer..}, rl:{{_Rule..}..}] :=
	DeleteDuplicates[Table[Sort[Replace[atomList,r,1]], {r,rl}]]

RotationalSymmetry[atomList:{_Integer..}, {}] := {atomList}


(*From a list of rotation lists, as generated by GetSymmetryMap*)

RotationalSymmetry[atomList:{_Integer..}, rotList:{_List..}] :=
	RotationalSymmetry[atomList, symmetryRules[rotList]]

(* Generate replacement rules from a list of atoms lists describing rotations.
   For example, for fumarate, {{1,2},{3,4}} indicates the rotations 1 -> 2 and 3 -> 4, and gives the rules
   {{1->1,2->2,3->3,4->4}, {1->2,2->1,3->4,4->3}}   *)

symmetryRules[rotList:{{_Integer..}..}]:=
	Table[
		Thread[Rule[
			Flatten[rotList],
			Join @@ Table[RotateLeft[r,k], {r,rotList}]]],
		{k, 0, Length[First[rotList]] - 1}]


End[];

EndPackage[];
