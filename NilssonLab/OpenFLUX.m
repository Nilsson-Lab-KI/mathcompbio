
(* This package provides functions for importing/exporting metabolic network models from/to OpenFLUX *)


BeginPackage["NilssonLab`OpenFLUX`",
	{"NilssonLab`Utilities`",
	 "NilssonLab`MetabolicNetwork`Common`", "NilssonLab`AtomMap`Common`","NilssonLab`EMUTracing`"}];


OpenFLUXExport::usage = "OpenFLUXExport[fileName, net, atomMaps, targets] exports a metabolite network model to OpenFLUX format (tab-separated text)";
OpenFLUXImport::usage = "OpenFLUXImport[fileName] ... ";


Begin["NilssonLab`OpenFLUX`Private`"];


(*
 * Export functions
 *)

(* In Openflux, metabolites are either balanced (internal) or not ("free"), 
   there are no "uptake" (Sv > 0) or "release" (Sv < 0) metabolites. To model these,
   we must add "input" reactions (suffix "_IN" ) for each "uptake" metabolite and
   "output" reactions (suffix "_OUT") for each "release" metabolite.
   Free metabolites have both an "_IN" and "_OUT" reaction. *)

OpenFLUXExport::stoich = "Stoichiometry > 1 is not supported";


(* Sort a list of atoms (with indices) according to a metabolite order.
   The order may include metabolites not present in the atom list. *)

atomOrdering[atomList:{{_Atom,_Integer}..}, order:{{_Metabolite,_Integer}..}] :=
	Block[{productOrder, nProducts, rule},
		productOrder = Join @@ Table[{p[[1]],i}, {p,order}, {i,p[[2]]}];
		nProducts = Length[productOrder];
		rule = Thread[Rule[productOrder, Range[nProducts]]];
		SortBy[Range[Length[atomList]],
			Replace[{Metabolite[atomList[[#,1]]], atomList[[#,2]]}, rule]&]]

atomOrdering[{}, order:{{_Metabolite,_Integer}..}] := {}


(* Generates atom strings in abcdef... format from metabolites listed in the atom map,
  in the order given by the substrates and products lists. *)

atomMapToStrings[amap_AtomMap,
			substrates:{{_Metabolite,_Integer}..}, products:{{_Metabolite,_Integer}..},
			symmetry:{Rule[_String,_List]...}]:=
	Block[{na, eductOrder, el, pl, productOrder, rotation},
		na = Length[EductList[amap]];
		(* ensure educt list is in substrate order *)
		el = EductList[amap];
		eductOrder = atomOrdering[el, substrates];
		(* educt atom list, grouped by substrate *)
		el = GatherBy[
			Transpose[{el[[eductOrder]], Take[CharacterRange["a","z"], na]}],
			{Metabolite[#[[1,1]]], #[[1,2]]} &];
		el = Map[{#[[1,1]], #[[2]]}&, el, {2}];
		(* products atom list, sorted and grouped by product and atom nr  *)
		pl = ProductList[amap][[eductOrder]];
		productOrder = atomOrdering[pl, products];
		pl = GatherBy[
			Transpose[{pl, Take[CharacterRange["a","z"], na]}][[productOrder]],
			{Metabolite[#[[1,1]]], #[[1,2]]} &];
		pl = Map[{#[[1,1]], #[[2]]}&, pl, {2}];
		pl = Table[SortBy[p, AtomNumber[#[[1]]]&], {p, pl}];
		(* generate lists of {metabolite, atom string} pairs, for substrates and products *)
		{Table[
			{MetaboliteShortName[Metabolite[e[[1,1]]]], StringJoin @@ e[[All,2]]},
			{e, el}],
		 Table[
			rotation = Replace[
				Replace[MetaboliteID[Metabolite[p[[1,1]]]], symmetry],
				 _String -> {}];
			productToStrings[p, rotation],
			{p, pl}]}];


(*Generate product atom strings, accounting for symmetry*)

productToStrings[product:{{_Atom, _String}..}, rotations:{_List..}] :=
	Block[{metName, rl, coeff},
		metName = MetaboliteShortName[Metabolite[product[[1,1]]]];
		rl = AtomRotations[AtomNumber /@ First /@ product, rotations];
		coeff = 1/Length[rl];
		Map[ListToString[#, " + "]&,
			Transpose[Table[
				{ToString[N[coeff]] <> " " <> metName,
				 ToString[N[coeff]] <> " " <> (StringJoin @@ product[[r,2]])},
				{r, rl}]]]]

productToStrings[product:{{_Atom, _String}..}, {}] :=
	{MetaboliteShortName[Metabolite[product[[1,1]]]], StringJoin @@ product[[All,2]]}


(* Convert atom map for one reaction to OpenFLUX format *)

(*forward direction*)

openFluxAtomMap[amap_AtomMap, net_MetabolicNetwork, reaction_Reaction, 1, symmetry:{Rule[_String,_List]...}]:=
	openFluxAtomMap[ReactionID[reaction], amap,
		ReactionSubstrates[net,reaction], ReactionProducts[net,reaction], symmetry]

(*reverse direction*)

openFluxAtomMap[amap_AtomMap, net_MetabolicNetwork, reaction_Reaction, -1, symmetry:{Rule[_String,_List]...}]:=
	openFluxAtomMap[ReactionID[reaction]<>"_R", Reverse[amap],
		ReactionProducts[net,reaction], ReactionSubstrates[net,reaction], symmetry]

(*From a list of reaction substrates and products (to add metabolites with no atoms mapped)*)

openFluxAtomMap[reactionId_String, amap_AtomMap, 
			substrates:{{_Metabolite,_Integer}..}, products:{{_Metabolite,_Integer}..},
			symmetry:{Rule[_String,_List]...}] :=
	Block[{ml,el,pl,rsl,rpl,order},
		(* generate string format from atom map *)
		{el, pl} = atomMapToStrings[amap, substrates, products, symmetry];
		(* find reactants and products not in the atom map *)
		{rsl, rpl} = {
			DeleteCases[substrates,
				{Alternatives@@Union[Metabolite/@First/@EductList[amap]],_}],
			DeleteCases[products,
				{Alternatives@@Union[Metabolite/@First/@ProductList[amap]],_}]};
		(* merge reactants, empty atom lists marked by "X" *)
		el =Join[el, Join@@Table[
			 Table[{MetaboliteShortName[First[m]],"X"},{Last[m]}],
			{m,rsl}]];
		pl=Join[pl,Join@@Table[
			Table[{MetaboliteShortName[First[m]],"X"},{Last[m]}],
			{m,rpl}]];
		(* generate table *)
		{reactionId,
		(* reaction stoichiometry *)
		ListToString[First/@el," + "]<>" = "<>ListToString[First/@pl," + "],
		(* atom map *)
		ListToString[Last/@el," + "]<>" = "<>ListToString[Last/@pl," + "]}]


(* Generate a table for all internal reactions.
   For now we ignore the rates, basis, and deviation fields *)

openFluxAtomInternalMapTable[atomMaps:{_AtomMap..}, net_MetabolicNetwork, symmetry:{Rule[_String,_List]...}]:=
	Join @@ Table[
		If[ReversibleQ[Reactions[net][[k]]],
			{Join[openFluxAtomMap[atomMaps[[k]],net,Reactions[net][[k]],1,symmetry],{"","FR","",""}],
			 Join[openFluxAtomMap[atomMaps[[k]],net,Reactions[net][[k]],-1,symmetry],{"","R","",""}]},
			{Join[openFluxAtomMap[atomMaps[[k]],net,Reactions[net][[k]],1,symmetry],{"","F","",""}]}],
		{k,Length[Reactions[net]]}]

(*For boundary reactions*)

openFluxBoundaryAtomMapTable[boundaryMaps:{_AtomMap..}, bnet_BoundaryNetwork]:=
	Block[{ml, S, fluxNames, j, substrate, product},
		ml = Cases[Metabolites[bnet], Metabolite[_,_,"Uptake"|"Release"|"Free"]];
		S =  Stoichiometry[bnet];
		fluxNames = BoundaryFluxNames[bnet];
		Join @@ Table[
			(* iterate over reactions for metabolite m (two reactions for Free metabolites) *)
			Table[
				substrate = First[Pick[Metabolites[bnet], S[[All,j]], -1]];
				product = First[Pick[Metabolites[bnet], S[[All,j]], 1]];
				Join[
					openFluxAtomMap[fluxNames[[j]], boundaryMaps[[j]], {{substrate, 1}}, {{product, 1}}, {}],
					{"", If[TrueQ[substrate == m], "B", "F"], "", ""}],
				{j, BoundaryFluxIndex[bnet, m]}],
		{m, ml}]]


(*Full table with header*)

openFluxAtomMapTable[
		atomMaps:{_AtomMap..}, net_MetabolicNetwork, boundaryMaps:{_AtomMap..}, bnet_BoundaryNetwork,
		symmetry:{Rule[_String,_List]...}]:=
	Join[
		{{"RxnID","rxnEq","rxnCTrans","rates","rxnType","basis","deviation"}},
		openFluxBoundaryAtomMapTable[boundaryMaps, bnet],
		openFluxAtomInternalMapTable[atomMaps, net, symmetry]]



(* Generate list of all metabolites that are not included in the S matrix.
   In our case, the boundary input and output metabolites are excluded. *)

openFluxExcludedMetabolites[bnet_BoundaryNetwork]:=
	Join[
		{{"##","excludedMetabolites"}},
		Thread[{"#",MetaboliteShortName/@BoundaryMetabolites[bnet]}]]


(* This corresponds to the EMU targets (typically the measured metabolites, in any compartment) *)

openFluxSimulatedMDVs[emus:{_EMU..}]:=
	Join[
		{{"##","simulatedMDVs"}},
		Table[
			{"#",MetaboliteShortName[Metabolite[e]]<>"#"<>(StringJoin@@Table["1",{EMUSize[e]}])},
			{e,emus}]]



(* Lists all network substrates, based on the boundary atom map *)

openFluxInputSubstrates[bmap:{_AtomMap..}]:=
	Block[{ml},
		ml=Cases[
			Union[Metabolite/@MetaboliteEMUs[bmap]],
			Metabolite[_,"Input",_]];
		Join[
			{{"##","inputSubstrates"}},
			Thread[{"#",MetaboliteShortName/@ml}]]]



(* Measurements table, not used*)

openFluxMeasurementsTable[]:={{"##","measurements"}}


(* Error table, not used *)

openFluxErrorTable[]:={{"##","error"}}


(*
 * Import functions
 *)

(* Split the openflux text format into sections *)

splitBySection[openflux_List] :=
	Block[{ix},
		ix=Join[{1},Flatten[Position[openflux[[All,1]],"##"]],{Length[openflux]+1}];
		Table[Take[openflux,{ix[[i]],ix[[i+1]]-1}],{i,Length[ix]-1}]]


(* Create a MetabolicNetwork from the openflux reaction table *)

importMetabolicNetwork[reactionTable_List] :=
	Block[{reactionIds, reversible, reactions, stoichiometry, atoms, net, balIndex, atomMaps},
		reactionIds = reactionTable[[All,1]];
		reversible = Replace[reactionTable[[All,5]], {"FR" -> True, _String -> False}, {1}];
		(* split reaction strings on = into substrate and product lists *)
		reactions = StringSplit[reactionTable[[All,2]], "="];
		(* Split reactant lists by +, then by space for explicit coefficients *)
		reactions = Map[StringSplit[#,"+"]&, reactions,{2}];
		reactions = Map[StringSplit[#," "]&, reactions,{3}];
		(* extract coefficients *)
		reactions = Replace[reactions, 
			{{m_String} -> {m, 1}, {c_String, m_String} :> {m, ToExpression[c]}}, {3}];
		(* get atom lists; similar to the above, but ignore coefficients  *)
		atoms = StringSplit[reactionTable[[All,3]], "="];
		atoms = Map[StringSplit[#,"+"]&, atoms,{2}];
		atoms = Map[Last[StringSplit[#," "]]&, atoms,{3}];
		(* gather by metabolite and sum coefficients *)
		stoichiometry = Map[GatherBy[#, First]&, reactions, {2}];
		stoichiometry = Map[{#[[1,1]], Round[Total[#[[All,2]]]]}&, stoichiometry, {3}];
		(* create MetabolicNetwork and AtomMaps *)
		{balIndex, net} = importMetabolicNetwork[reactionIds, reversible, stoichiometry];
		atomMaps = importAtomMaps[net, reactions[[balIndex]], atoms[[balIndex]]];
		{net, atomMaps}]

importMetabolicNetwork[reactionIds_List, reversible_List, stoichiometry_List] :=
	Block[{inputIndex, outputIndex, balIndex, S, balMetIds, inputMetIds, outputMetIds, exchange},
		(* input, output and mass balanced (internal) reactions *)
		inputIndex = Flatten[Position[reactionIds, s_String/;StringMatchQ[s,"*_IN"]]];
		outputIndex = Flatten[Position[reactionIds, s_String/;StringMatchQ[s,"*_OUT"]]];
		balIndex = Complement[Range[Length[reactionIds]],
			inputIndex, outputIndex,
			Flatten[Position[reactionIds, s_String/;StringMatchQ[s,"*_R"]]]];
		(* get balanced metabolite list *)
		balMetIds = Union[Flatten[Map[First, stoichiometry[[balIndex]], {3}]]];
		(* get metabolite exchange status from IN / OUT reactions *)
		inputMetIds = Union[Flatten[Map[First, stoichiometry[[inputIndex]], {3}]]];
		outputMetIds = Union[Flatten[Map[First, stoichiometry[[outputIndex]], {3}]]];
		exchange = Replace[balMetIds, Dispatch[Join[
			Thread[Rule[Intersection[inputMetIds, outputMetIds], "Free"]],
			Thread[Rule[inputMetIds, "Uptake"]],
			Thread[Rule[outputMetIds, "Release"]],
			{_String -> "Internal"}]], {1}];
		(* make stoichiometry matrix *)
		S = stoichiometry[[balIndex]] /. Dispatch[Thread[Rule[balMetIds, Range[Length[balMetIds]]]]];
		S = SparseArray[ 
			Flatten[Table[
				{S[[i,1]]/.{mi_Integer,c_Integer}:>Rule[{mi,i},-c], S[[i,2]]/.{mi_Integer,c_Integer}:>Rule[{mi,i},c]},
				{i, Length[balIndex]}]],
			{Length[balMetIds], Length[balIndex]}];
		{balIndex, MetabolicNetwork[
			Thread[Reaction[reactionIds[[balIndex]], reversible[[balIndex]]]],
			Thread[FromMetaboliteShortName[balMetIds, exchange]],
			S]}]

importTargets[net_MetabolicNetwork, targetsTable_List] :=
	Block[{metNames, atoms},
		{metNames, atoms} = Transpose[Map[StringTrim[StringSplit[#, "#"]]&, targetsTable[[All, 2]]]];
		Thread[EMU[
			Table[FindMetabolite[net, m], {m, metNames}],
			Table[{{First/@StringPosition[a,"1"]}}, {a, atoms}]]]]

importAtomMaps[net_MetabolicNetwork, reactions_List, atoms_List] :=
	Block[{metList, metNameMap},
		metList = Join[Metabolites[net], BoundaryMetabolites[net]];
		metNameMap = Thread[Rule[MetaboliteShortName /@ metList, metList]];
		Table[
			importAtomMap @@ filterAtoms[Map[First, reactions[[i]], {2}] /. metNameMap, atoms[[i]]],
			{i, Length[reactions]}]]


(* Create an AtomMap from a OpenFLUX-style list of educt and product atoms.
   Here we assume no stoichiometry is > 1 !  Any repeated metabolites in the reaction list
    will be skipped, including those with non-integer stoichiometry encoding symmetry. *)

atomLetterRules = Table[Rule[FromCharacterCode[96 + i], i], {i, 26}];	

importAtomMap[reaction_List, atoms:{{_String..},{_String..}}] :=
	Block[{eductLetters, eductNumbers, eductAtoms, productIndex, productNumbers, productAtoms, metAtoms},
		(* educt atoms *)
		eductLetters = Map[Characters, atoms[[1]]];
		eductNumbers = eductLetters /. atomLetterRules /. {x:{_Integer..} :> x - Min[x] + 1};
		eductAtoms = Flatten @ Map[
			Thread[Atom[First[#], "C", Last[#]]]&,
			Transpose[{reaction[[1]], eductNumbers}]];
		(* product atoms *)
		productNumbers = Map[Range[StringLength[#]]&, atoms[[2]]]; 
		productAtoms = Flatten @ Map[
			Thread[Atom[First[#], "C", Last[#]]]&,
			Transpose[{reaction[[2]], productNumbers}]];
		(* reorder product atoms according to atomsDel list *)
		productIndex = Characters[StringJoin @@ Flatten[atoms[[2]]]] /. atomLetterRules;
		productAtoms[[productIndex]] = productAtoms;
		(* generate AtomMap *)
		AtomMap[Thread[{eductAtoms, 1}], Thread[{productAtoms, 1}]]]

importAtomMap[{{},{}}, {{},{}}] := AtomMap[{}, {}]


(* Drop duplicates and non-mapped metabolites ("X") *)

filterAtoms[reaction_List, atoms:{{_String..},{_String..}}] :=
	Block[{keepIndex},
		Transpose[Table[
			keepIndex = Intersection[
				ReindexList[DeleteDuplicates[reaction[[i]]], reaction[[i]]],
				Flatten[Position[atoms[[i]], Except["X"], Heads->False]]];
			{reaction[[i, keepIndex]], atoms[[i, keepIndex]]},
			{i, 2}]]]


(* Export with symmetry map *)

OpenFLUXExport[fileName_String,
			net_MetabolicNetwork, atomMaps:{_AtomMap..}, symmetry:{Rule[_String,_List]...}, targets:{_EMU..}]:=
	OpenFLUXExport[fileName, net, atomMaps, BoundaryNetwork[net], BoundaryAtomMaps[net, atomMaps], symmetry, targets]

OpenFLUXExport[fileName_String,
		net_MetabolicNetwork, atomMaps:{_AtomMap..}, bnet_BoundaryNetwork, boundaryMaps:{_AtomMap..},
		symmetry:{Rule[_String,_List]...}, targets:{_EMU..}]:=
	Export[fileName,
		Join[
			openFluxAtomMapTable[atomMaps, net, boundaryMaps, bnet, symmetry],{{""}},
			openFluxExcludedMetabolites[bnet],{{""}},
			openFluxSimulatedMDVs[targets],{{""}},
			openFluxInputSubstrates[boundaryMaps],{{""}},
			openFluxMeasurementsTable[],{{""}},
			openFluxErrorTable[]],
		"TSV"]


(* Import .*)

OpenFLUXImport[fileName_String] :=
	Block[{openflux,reactions,excluded,targets,inputs, net, atomMaps},
		openflux = Import[fileName, "TSV"];
		(* remove any empty rows *)
		openflux = DeleteCases[openflux,{""...}];
		(* split into sections, take the first four, drop headers *)
		{reactions,excluded,targets,inputs} = Rest /@ Take[splitBySection[openflux],4];
		{net, atomMaps} = importMetabolicNetwork[reactions];
		{net, atomMaps, importTargets[net, targets]}]


End[];

EndPackage[];
