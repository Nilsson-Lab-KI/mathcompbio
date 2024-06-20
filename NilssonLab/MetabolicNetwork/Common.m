(* This package contains data structures for representing metabolic networks. *)

BeginPackage["NilssonLab`MetabolicNetwork`Common`",	{"NilssonLab`Utilities`"}];

MetabolicNetwork::usage = "This symbol is used for data structures describing a metabolic network model.";
Reactions::usage = "Reaction[M] returns the list of reactions in a MetabolicNetwork M";
Metabolites::usage = "Metabolites[M] returns the list of reactions in a MetabolicNetwork M";
Reaction::usage = "This symbol is used to describe a reaction.";
ReactionID::usage = "ReactionID[r] where r is a Reaction returns the ID string of the reaction.";
ReversibleQ::usage = "ReversibleQ[r] where r is a Reaction returns True if r is reversible.";
ReactionSubstrates::usage = "ReactionSubstrates[M, r] where M is a MetabolicNetwork and r is a Reaction returns a list of reaction substrates and their stoichometry coefficients.";
ReactionProducts::usage = "ReactionProducts[M, r] where M is a MetabolicNetwork and r is a Reaction returns a list of reaction products and their stoichometry coefficients.";
Metabolite::usage = "This symbol is used to describe a metabolite.";
MetaboliteID::usage = "MetaboliteID[m] where m is a Metabolite returns the ID string of the metabolite. Note that this is not a unique identifier in case of multiple compartments; see MetaboliteShortName.";
MetaboliteProducers::usage = "MetaboliteProducers[M, m] where M is a MetabolicNetwork and m is a Metabolite returns a list of reactions producing m and their stoichometry coefficients.";
MetaboliteConsumers::usage = "MetaboliteConsumers[M, m] where M is a MetabolicNetwork and m is a Metabolite returns a list of reactions consuming m and their stoichometry coefficients.";
ExchangeStatus::usage = "ExchangeStatus[m] returns the exchange status of Metabolite m, which is either \"Free\", \"Uptake\", \"Release\", or \"Internal\".";
Compartment::usage = "Compartment[m] returns the compartment of Metabolite m.";
Stoichiometry::usage = "Stoichiometry[M] returns the stoichiometry matrix of a MetabolicNetwork M.";
ReverseStoichiometry::usage = "ReverseStoichiometry[M] returns the stoichiometry matrix describing the reversible reactions in M, in the reverse direction.";
MetaboliteShortName::usage = "MetaboliteShortName[m] returns an string providing a unique identifier of the Metabolite m, considering compartments.";
FromMetaboliteShortName;
FindMetabolite::usage = "FindMetabolite[M, s] returns the metabolite m in the MetabolicNetwork M with MetaboliteShortName[m] == s, if it exists.";
DeleteReactions::usage = "DeleteReactions[M, r] where M is a MetabolicNetwork returns a new network where the reaction(s) r have been deleted.";
DeleteMetabolites::usage = "DeleteReactions[M, m] where M is a MetabolicNetwork returns a new network where the metabolites(s) m have been deleted.";
AddReactions::usage = "AddReactions[M, r] where M is a MetabolicNetwork returns a new network where the reaction(s) r have been added.";
MergeReactions::usage = "MergeReactions[M, rl] where M is a MetabolicNetwork and rl is a List of reactions a new network where the reactions have been merged into the first reaction.";
SetExchange::usage = "SetExchange[M, m] returns a new MetaboliteNetwork where the exchange status of Metabolite m is as specified in m.";
ReactionFormula::usage = "ReactionFormula[M, r] returns a reactions formula in string form for a Reaction r.";
ExchangeMetabolites::usage = "ExchangeMetabolites[M] returns the list of metabolites in a MetabolicNetwork M with exchange status \"Free\", \"Uptake\", or \"Release\".";
BoundaryReactions;
InputMetabolites;
OutputMetabolites;
BoundaryMetabolites;
BoundaryNetwork;
NumberOfBoundaryFluxes;
BoundaryFluxNames;
UptakeFluxIndex;
BoundaryFluxIndex;
CompleteStoichiometry;
MetaboliteReactionGraph;
SetReversible;
SetIrreversible;
TransporterQ;
NetFluxState;
NetFlux;
MassBalance;
MassBalanceQ;
MassBalanceNorm;
FluxState;
ForwardFlux;
ReverseFlux;
BoundaryUptake;
BoundaryRelease;
BoundaryFlux;
ExchangeFlux;
ExchangeCoefficient;
FindFreeFluxes;
FreeFluxBasis;
FindViableReactions;
MetaboliteTotalFlux;
NumberOfFluxes;


(* ::Subsection:: *)
(*Begin private context*)


Begin["NilssonLab`MetabolicNetwork`Common`Private`"];


(*
 * Metabolic networks
 *)


(*
 * A MetabolicNetwork is a list {reactions, metabolites, stochiometry}.
 * We currently require integer stoichiometry, since  atom maps cannot be constructed otherwise.
 *)

Format[MetabolicNetwork[rl_List, ml_List, S_SparseArray]] :=
	"<MetabolicNetwork, "<>ToString[Length[rl]]<>" reactions>"

Reactions[MetabolicNetwork[rl_List, ml_List, S_SparseArray]] ^:= rl

Metabolites[MetabolicNetwork[rl_List, ml_List, S_SparseArray]] := ml


(* This gives the stoichiometry matrix of a reaction *)

Stoichiometry[MetabolicNetwork[rl_List, ml_List, S_SparseArray]] := S


(* This gives a row of the stoichiometry matrix for a given metabolite *)

Stoichiometry::nomet = "No such metabolite `1`";

Stoichiometry[MetabolicNetwork[rl_List, ml_List, S_SparseArray], m_Metabolite] :=
	Block[{i},
		i = Position[ml, m];
		If[Length[i] != 1, 
			Message[Stoichiometry::nomet, m]; $Failed,
			Part[S, i[[1,1]]]]]


(* This gives a stoichiometry submatrix for multiple metabolites *)

Stoichiometry::metmatch = "Some metabolites do not match.";

Stoichiometry[MetabolicNetwork[rl_List, ml_List, S_SparseArray], m:{_Metabolite..}] :=
	Block[{i},
		i = ReindexList[m, ml];
		If[MemberQ[i, Missing], Message[Stoichiometry::metmatch]];
		Part[S, Cases[i, _Integer]]]


(*This gives a column of the stoichiometry matrix for a given reaction*)

Stoichiometry::noreact = "No such reaction `1`";


Stoichiometry[MetabolicNetwork[rl_List, ml_List, S_SparseArray], r_Reaction] :=
	Block[{i},
		i = Position[rl, r];
		If[Length[i] != 1, 
			Message[Stoichiometry::noreact, r]; $Failed,
			S[[All, i[[1,1]]]]]]


(*This gives a submatrix for multiple reactions*)

Stoichiometry::reactmatch = "Some reactions do not match.";

Stoichiometry[MetabolicNetwork[rl_List, ml_List, S_SparseArray], r:{_Reaction..}] :=
	Block[{i},
		i = Cases[ReindexList[r, rl], _Integer];
		If[MemberQ[Length[i], Missing], 
			Message[Stoichiometry::reactmatch]];
		Part[S, All, Cases[i, _Integer]]]


(*This is the stoichiometry matrix corresponding to the reverse of the reversible reactions.*)

ReverseStoichiometry[M_MetabolicNetwork] :=
	Transpose[Transpose[Stoichiometry[M]]
		* SparseArray[ReversibleQ[Reactions[M]] /. {True-> -1, False-> 0}]]



(*
 * A Reaction has an identifier and a reversibility flag.
 * We pretty-print reaction with a "(rev)" decoration if they are reversible.
  *)

Format[Reaction[id_String, r:(True|False)]] :=
	"<Reaction " <> id <> If[r, "(rev)", ""] <> ">"

ReactionID[Reaction[id_String, r:(True|False)]] := id

ReactionID[rl:{_Reaction..}] := rl[[All,1]]

ReversibleQ[Reaction[id_String, r:(True|False)]] := r

ReversibleQ[rl:{_Reaction..}] := rl[[All,2]]


(*
 * A Metabolite has an ID, a compartment, and an exchange status (Internal, Uptake, Release, or Free)
 *)

MetaboliteID[Metabolite[id_,c_,e_]] := id

MetaboliteID[ml:{_Metabolite..}] := ml[[All,1]]

Compartment[Metabolite[id_,c_,e_]] := c

Compartment[ml:{_Metabolite..}] := ml[[All,2]]

ExchangeStatus[Metabolite[id_,c_,e_]] := e

ExchangeStatus[ml:{_Metabolite..}] := ml[[All,3]]

(* This produces a String which can be used to unique identify a metabolite in a given compartment. *)

MetaboliteShortName[Metabolite[id_,c_,e_]] := MetaboliteShortName[{id, c}]


(* List of known compartments. The Input and Output compartments represent boundary metabolites
 * which are strictly not part of the metabolic network but may be required for isotopic simulations.
 * The Virtual compartment is used for dummy metabolites which arise in some network representations,
 * notably decomposition into binary networks.
 *)

MetaboliteShortName[{met_String, "Cytosol"}] := met <> "_c"
MetaboliteShortName[{met_String, "EndoplasmicReticulum"}] := met <> "_r"
MetaboliteShortName[{met_String, "Extracellular"}] := met <> "_e"
MetaboliteShortName[{met_String, "Golgi"}] := met <> "_g"
MetaboliteShortName[{met_String, "Lysosome"}] := met <> "_l"
MetaboliteShortName[{met_String, "Mitochondria"}] := met <> "_m"
MetaboliteShortName[{met_String, "Nucleus"}] := met <> "_n"
MetaboliteShortName[{met_String, "Peroxisome"}] := met <> "_x"
MetaboliteShortName[{met_String, "Input"}] := met <> "_i"
MetaboliteShortName[{met_String, "Output"}] := met <> "_o"
MetaboliteShortName[{met_String, "Virtual"}] := met <> "_v"

(*Parse a metabolite short name to a Metabolite*)

FromMetaboliteShortName[metName_String] :=
	Block[{metId, suffix},
		metId = StringDrop[metName, -2];
		suffix = StringTake[metName, -1];
		Metabolite[metId, suffixToCompartment[suffix], suffixToExchange[suffix]]]


FromMetaboliteShortName[metName_String, exhange_String] :=
	Block[{metId, suffix},
		metId = StringDrop[metName, -2];
		suffix = StringTake[metName, -1];
		Metabolite[metId, suffixToCompartment[suffix], exhange]]

suffixToCompartment[suffix_String] :=
	Replace[suffix, {
		"c" -> "Cytosol", "r" -> "EndoplasmicReticulum", "e" -> "Extracellular", "g" -> "Golgi",
		"l" -> "Lysosome", "m" -> "Mitochondria", "x" -> "Peroxisome",
		"i" -> "Input",  "o" -> "Output", "v" -> "Virtual"}]


(* Set default exchange status to "Free" or "Boundary" *)

suffixToExchange[suffix_String] :=
	Replace[suffix, {
		"c" | "r" | "e" | "g" | "l" | "m" | "x" -> "Free",
		"i" | "o" | "v" -> "Boundary"}]


(* Locate a Metabolite in a MetabolicNetwork from a MetaboliteShortName *)

FindMetabolite::notfound = "Metabolite `1` not found.";

FindMetabolite[MetabolicNetwork[rl_List, ml_List, S_SparseArray], mshort_String] :=
	Block[{mi},
		mi = Position[MetaboliteShortName /@ ml, mshort];
		If[Length[mi] == 0,
			Message[FindMetabolite::notfound, mshort]; $Failed,
			ml[[mi[[1,1]]]]]]

FindMetabolite::missing = "There were missing metabolites.";

FindMetabolite[MetabolicNetwork[rl_List, ml_List, S_SparseArray], mshort:{_String..}] :=
	Block[{mi},
		mi = ReindexList[mshort, MetaboliteShortName /@ ml];
		If[MemberQ[mi, Missing],
			Message[FindMetabolite::missing];
			mi = Cases[mi, _Integer]];
		ml[[mi]]]


(* Alter a MetabolicNetwork so that a given reaction ID is reversible. Returns a modified network. *)

SetReversible[M:MetabolicNetwork[rl_List, ml_List, S_SparseArray], rids_List] :=
	 MetabolicNetwork[
		rl /. Thread[Rule[
			Thread[Reaction[rids, _]],
			Thread[Reaction[rids ,True]]]],
		ml, S]

SetReversible[M_MetabolicNetwork, rid_String] := SetReversible[M, {rid}]


(*Same, indexing into the reaction list. *)

SetReversible[MetabolicNetwork[rl_List, ml_List, S_SparseArray], i_Integer] :=
	 MetabolicNetwork[
		Join[Take[rl,i-1], 
			{Reaction[ReactionID[rl[[i]]], True]},
			Take[rl,{i+1, Length[rl]}]],
		 ml, S]


(* Alter a MetabolicNetwork so that a given reaction ID is irreversible. Retuns a modified network. *)

SetIrreversible[M:MetabolicNetwork[rl_List, ml_List, S_SparseArray], rids_List] :=
	 MetabolicNetwork[
		rl /. Thread[Rule[
			Thread[Reaction[rids, _]],
			Thread[Reaction[rids, False]]]],
		ml, S]

SetIrreversible[M_MetabolicNetwork, rid_String] := SetIrreversible[M, {rid}]

SetIrreversible[MetabolicNetwork[rl_List, ml_List, S_SparseArray], i_Integer] :=
	 MetabolicNetwork[
		Join[Take[rl,i-1], 
			{Reaction[ReactionID[rl[[i]]], False]},
			Take[rl,{i+1, Length[rl]}]],
		 ml, S]

(* Delete reactions and any metabolites that become unconnected as a result *)

DeleteReactions[
		M:MetabolicNetwork[rl_List, ml_List, S_SparseArray], rldel:{_Reaction..}] :=
	DeleteReactions[M, ReindexList[rldel, rl]]

DeleteReactions[M:MetabolicNetwork[rl_List, ml_List, S_SparseArray], rids:{_String..}] :=
	DeleteReactions[M, ReindexList[rids, ReactionID[rl]]]

DeleteReactions[M:MetabolicNetwork[rl_List, ml_List, S_SparseArray], ri:{_Integer..}] :=
	Block[{ric, mi},
		ric = Complement[Range[Length[rl]], ri];
		mi = Flatten[Position[Total /@ Abs[S[[All, ric]]], _?Positive]];
		MetabolicNetwork[rl[[ric]], ml[[mi]], S[[mi, ric]]]]


(* Delete metabolites and any reactions that become unconnected as a result *)

DeleteMetabolites[M_MetabolicNetwork, {}] := M

DeleteMetabolites[
		M:MetabolicNetwork[rl_List, ml_List, S_SparseArray], mldel:{_Metabolite..}] :=
	DeleteMetabolites[M, ReindexList[mldel, ml]]

DeleteMetabolites[M:MetabolicNetwork[rl_List, ml_List, S_SparseArray], mi:{_Integer..}] :=
	Block[{mic, ri},
		mic = Complement[Range[Length[ml]], mi];
		ri = Flatten[Position[Total[Abs[S[[mic]]]], _?Positive]];
		MetabolicNetwork[rl[[ri]], ml[[mic]], S[[mic, ri]]]]


(* Create a new metabolic network by adding a reaction.
 * Here met is a list of the metabolites involved, and c is hte corresponding stoichiometry coefficients.
 *)

AddReactions[M:MetabolicNetwork[rl_List, ml_List, S_SparseArray],
			r_Reaction, met:{_Metabolite..}, s_?VectorQ] :=
	Block[{nr, nm, newMet, mi},
		nr = Length[rl];
		nm = Length[ml];
		newMet = Complement[met, ml];
		mi = ReindexList[met, Join[ml, newMet]];
		MetabolicNetwork[
			Append[rl, r], Join[ml, newMet],
			SparseArray[
				Join[Most[ArrayRules[S]],
					Thread[Rule[Thread[{mi, nr+1}],s]]],
				{nm + Length[newMet], nr + 1}]]]


(* Merge or "lump" 2 or more reactions, and delete any metabolites that become disconnected as a result. *)

MergeReactions[
		M:MetabolicNetwork[rl_List, ml_List, S_SparseArray], rl:{_Reaction..}] :=
	MergeReactions[M, ReindexList[rl, rl]]

MergeReactions[
		M:MetabolicNetwork[rl_List, ml_List, S_SparseArray], ri:{_Integer..}] :=
	Block[{Snew},
		If[Length[ri] < 2, Return[M]]; (* nothing to do *)
		Snew = Stoichiometry[M];
		Snew[[All,First[ri]]] = Total /@ Snew[[All,ri]];
		DeleteReactions[MetabolicNetwork[rl, ml, Snew], Rest[ri]]]


(* A reaction is a transporter if its products and substrates cancel when compartments are discarded. *)

TransporterQ[M_MetabolicNetwork, r_Reaction] :=
	Block[{i, Si, ml},
		Si = Stoichiometry[M, r];
		ml = MetaboliteID[Metabolites[M]];
		TrueQ[Si.ml == 0]]


(* Return a list of metabolites that are substrates (reactants, educts) for a given reaction,
 * together with their stoichiometry coefficients
 *)

ReactionSubstrates[M:MetabolicNetwork[rl_List, ml_List, S_SparseArray], r_Reaction] :=
	Cases[Thread[{ml, -Stoichiometry[M, r]}], {_, _?Positive}]



(* Return a list of metabolites that are substrates (reactants, educts) for a given reaction,
 * together with their stoichiometry coefficients
 *)

ReactionProducts[M:MetabolicNetwork[rl_List, ml_List, S_SparseArray], r_Reaction] :=
	Cases[Thread[{ml, Stoichiometry[M, r]}], {_, _?Positive}]

(* Return a list of reactions for which the given metabolite is a product, together with
   their stoichiometry coefficients. NOTE: this does not account for reversibility! *)

MetaboliteProducers[M:MetabolicNetwork[rl_List, ml_List, S_SparseArray], m_Metabolite] :=
	Cases[Thread[{rl, Stoichiometry[M, m]}], {_, _?Positive}]

(* Return a list of reactions for which the given metabolite is a substrate (reactant, educt),
   together with their stoichiometry coefficients *)

MetaboliteConsumers[M:MetabolicNetwork[rl_List, ml_List, S_SparseArray], m_Metabolite] :=
	Cases[Thread[{rl, -Stoichiometry[M, m]}], {_, _?Positive}]


(* Return the list of metabolites that may exchange with the system boundary.
   Note that this list does not correspond BoundaryMetabolites[M] since metabolites in M
   may have both uptake and release.*)

ExchangeMetabolites[M_MetabolicNetwork] :=
	Cases[Metabolites[M], Metabolite[_, _, "Uptake" | "Release" | "Free"]]


(* Format reactions into readable formulas *)

ReactionFormula[M_MetabolicNetwork, r_Reaction]:=
	First[ReactionFormula[M, {r}]]

ReactionFormula[M_MetabolicNetwork] :=
	ReactionFormula[M, Reactions[M]]

ReactionFormula[M_MetabolicNetwork, rl:{_Reaction..}] :=
	Block[{S, mids}, 
		S = Stoichiometry[M, rl];
		mids = MetaboliteShortName /@ Metabolites[M];
		Table[
			StringJoin[
				ToString[NegativePart[S[[All,i]]].mids],
				If[ReversibleQ[rl[[i]]], " <=> ", " ==> "],
				ToString[PositivePart[S[[All,i]]].mids]],
			{i, Length[rl]}]]



Options[FindViableReactions] = {Tolerance -> 0.001};

FindViableReactions[M_MetabolicNetwork, opts:OptionsPattern[]] :=
	Block[{bs, em, Se, nr, xmax = 1000, ul, c},
		nr = Length[Reactions[M]];
		(* Stoichiometry matrix minus the Free metabolites *)
		em = Cases[Metabolites[M], Metabolite[_,_,Except["Free"]]];
		Se = N[Stoichiometry[M, em]];
		(* right hand side from metabolite boundary status *) 
		bs = em /. {
			Metabolite[_,_,"Internal"]->{0,0},
			Metabolite[_,_,"Release"]->{0,1},
			Metabolite[_,_,"Uptake"]->{0,-1}};
		(* reaction bounds from reversibility *) 
		ul = ReversibleQ[Reactions[M]] /. {True -> {-xmax, xmax}, False -> {0, xmax}};
		(* improve scaling
		c = Clip[Mean /@ Se, {1, \[Infinity]}] *)
		Table[
			NetFluxState[
				Check[
					Chop[LinearProgramming[
						s*UnitVector[nr, j],
						Se, bs, ul, opts]],
					(* return zero vector if infeasible *)
					Table[0,{nr}]]],
			{j, nr}, {s, {1, -1}}]]


(* Make a list of metabolites in the "virtual" boundary (input/output) compartments.
   These are not part of the stoichiometry system (they do not satisfy mass balance)
   but they are needed as inputs for isotope tracing algorithms. Note that there may be
   more BoundaryMetabolites than ExchangeMetabolites since each internal metabolite may have
   both uptake and release boundary reactions.*)

BoundaryMetabolites[M_MetabolicNetwork] :=
	Join[InputMetabolites[M], OutputMetabolites[M]]

InputMetabolites[M_MetabolicNetwork] :=
	InputMetabolites[Metabolites[M]]

InputMetabolites[ml:{_Metabolite..}] :=
	Cases[ml, Metabolite[id_, _, "Uptake" | "Free"] -> Metabolite[id, "Input", "Boundary"]]

OutputMetabolites[M_MetabolicNetwork] :=
	OutputMetabolites[Metabolites[M]]

OutputMetabolites[ml:{_Metabolite..}] :=
	Cases[ml, Metabolite[id_, _, "Release" | "Free"] -> Metabolite[id, "Output", "Boundary"]]



(* Determines a linear basis in terms of free fluxes using the row reduction method.
   Return a pair {vi, bi}, where vi are indices into the internal fluxes, and bi are indices
   into the boundary fluxes. *)

FindFreeFluxes[M_MetabolicNetwork] :=
 	FindFreeFluxes[M, BoundaryNetwork[M]]

FindFreeFluxes[M_MetabolicNetwork, B_BoundaryNetwork] :=
 	Block[{nm, nr, nb, Sfull, Sred, ix},
  		{nm, nr} = Dimensions[Stoichiometry[M]];
  		nb = NumberOfBoundaryFluxes[B];
  		Sfull = ArrayFlatten[{{
      			Stoichiometry[M], Take[Stoichiometry[B], nm]}}];
  		ix = FindFreeFluxes[Sfull];
  		{Cases[ix, x_ /; x <= nr], Cases[ix, x_ /; x > nr -> x - nr]}]

FindFreeFluxes[S_SparseArray] :=
 	Block[{nr, Sred, ix},
  		nr = Length[First[S]];
  		Sred = SparseArray[RowReduce[S]];
  		ix = NonZeroIndices /@ Sred;
  		ix = Flatten[First /@ Cases[ix, x_List /; Length[x] > 0]];
  		Complement[Range[nr], ix]]


(* Find a null space basis in terms of free fluxes *)

FreeFluxBasis[M_MetabolicNetwork] :=
	FreeFluxBasis[M, BoundaryNetwork[M]]

FreeFluxBasis[M_MetabolicNetwork, B_BoundaryNetwork] :=
	Block[{nm, nr, nb, Sfull, Sred, ix, K},
		{nm, nr} = Dimensions[Stoichiometry[M]];
		nb = NumberOfBoundaryFluxes[B];
		Sfull = ArrayFlatten[{{
			Stoichiometry[M], Take[Stoichiometry[B], nm]}}];
		{ix, K} = FreeFluxBasis[Sfull];
		{Cases[ix, x_ /; x <= nr], Cases[ix, x_ /; x > nr -> x - nr], K}]

FreeFluxBasis[S_SparseArray] :=
	Block[{m, n, Sred, dep, nd, free, nf, K},
		{m,n} = Dimensions[S];
		(* calculate via row reduction *)
		Sred = SparseArray[RowReduce[S]];
		(* first nonzero position for each row are dependent *)
		dep = Flatten[First /@ Cases[NonZeroIndices /@ Sred, Except[{}]]];
		nd = Length[dep];
		free = Complement[Range[n], dep]; (* free fluxes *)
		nf = Length[free];
		K = SparseArray[{}, {n, nf}];
		K[[dep]] = -Take[Sred, nd][[All, free]];
		K[[free]] = SparseIdentity[nf];
		{free, K}]


(* This gives a bipartite directed graph connecting reactions with metabolites.
   Reversible reactions are split into separate nodes for forward and reverse reactions,
   as otherwise it is not possible to recover the reactions from the graph. *)

MetaboliteReactionGraph[M_MetabolicNetwork] :=
	Block[{S, revi, m, n},
		S = Stoichiometry[M];
		revi = Flatten[Position[ReversibleQ[Reactions[M]],True]];
		If[Length[revi] > 0,
			S = ArrayFlatten[{{S, ReverseStoichiometry[M][[All, revi]]}}]];
		{m,n} = Dimensions[S];
		AdjacencyGraph[
			Join[MetaboliteShortName /@ Metabolites[M],
				ReactionID /@ Reactions[M],
				Table[ReactionID[r] <> "_r", {r, Reactions[M][[revi]]}]],
			Sign[ArrayFlatten[
				{{SparseZero[m,m], NegativePart[S]},
				 {Transpose[PositivePart[S]], SparseZero[n, n]}}]]]]


(* Set the exchange status of a given metabolite; may be either "Internal", "Uptake", "Release" or "Free". *)

SetExchange[M_MetabolicNetwork, m:Metabolite[id_,c_,e_]] :=
	Block[{mi, ml},
		ml = Metabolites[M] /. {Metabolite[id, c, _] -> m};
		MetabolicNetwork[Reactions[M], ml, Stoichiometry[M]]]

SetExchange[M_MetabolicNetwork, mlist:{_Metabolite..}] :=
	Block[{mi, ml},
		ml = Metabolites[M] /. Table[
			Metabolite[MetaboliteID[m], Compartment[m], _] -> m,
			{m, mlist}];
		MetabolicNetwork[Reactions[M], ml, Stoichiometry[M]]]


(* A BoundaryNetwork is a stochiometric network handling boundary fluxes.*)

Format[BoundaryNetwork[ml:{_Metabolite..}, S_SparseArray]] :=
	"<BoundaryNetwork, "<> ToString[Length[First[S]]] <> " boundary fluxes>"

Metabolites[BoundaryNetwork[ml:{_Metabolite..}, S_SparseArray]] := ml

Stoichiometry[BoundaryNetwork[ml:{_Metabolite..}, S_SparseArray]] := S

NumberOfBoundaryFluxes[B_BoundaryNetwork] :=
	Dimensions[Stoichiometry[B]][[2]]

BoundaryMetabolites[B_BoundaryNetwork] :=
	Take[Metabolites[B], -NumberOfBoundaryFluxes[B]]


(*Construct a BoundaryNetwork that contains input and output fluxes.*)

BoundaryNetwork[M_MetabolicNetwork] := 
	Block[{nm, ml, uptakeIndex, releaseIndex, nu, nr},
		ml = Metabolites[M];
		nm = Length[ml];
		uptakeIndex = Flatten[Position[Metabolites[M], Metabolite[_, _, "Uptake" | "Free"]]];
		nu = Length[uptakeIndex];
		releaseIndex = Flatten[Position[Metabolites[M], Metabolite[_, _, "Release" | "Free"]]];
		nr = Length[releaseIndex];
		BoundaryNetwork[
			Join[Metabolites[M], BoundaryMetabolites[M]],
			SparseArray[
				Flatten[
					{Table[
						{Rule[{uptakeIndex[[i]], i}, 1], Rule[{nm+i, i}, -1]},
						{i, nu}],
					Table[
						{Rule[{releaseIndex[[i]], nu+i}, -1], Rule[{nm+nu+i, nu+i}, 1]},
						{i, nr}]}],
				{nm + nu + nr, nu + nr}]]]


(* Returns the uptake boundary flux index of a metabolite m, if it is an Uptake or Free metabolite,
   otherwise Null *)

UptakeFluxIndex::nomet = "No such metabolite `1`";

UptakeFluxIndex[B_BoundaryNetwork, m:Metabolite[id_,_,e_]] :=
	Block[{Sm},
		If[MemberQ[{"Uptake","Free"}, e],
			Sm = Extract[Stoichiometry[B], Position[Metabolites[B], m]];
			If[Length[Sm] == 0,
				Message[UptakeFluxIndex::nomet, m]; Return[$Failed]];
			NonZeroIndices[PositivePart[Sm[[1]]]][[1,1]]]]


(* Return the boundary flux indexes of a metabolite m *)

BoundaryFluxIndex[B_BoundaryNetwork, m_Metabolite] :=
	Block[{Sm},
		Sm = Extract[Stoichiometry[B], Position[Metabolites[B], m]];
		If[Length[Sm] == 0,
			Message[BoundaryFluxIndex::nomet, m]; Return[$Failed]];
		Flatten[NonZeroIndices[Sm[[1]]]]]

(* Returns names for boundary fluxes, like metabolite_IN and metabolite_OUT. *)

BoundaryFluxNames[B_BoundaryNetwork] :=
	BoundaryFluxNames[BoundaryMetabolites[B]]

BoundaryFluxNames[ml:{Metabolite[_,_,"Boundary"]...}] :=
	StringReplace[MetaboliteShortName /@ ml, {"_i" -> "_IN", "_o" -> "_OUT"}]


(* The "complete" stoichiometry matrix, defined over
   (internal + boundary metabolites) x (internal + boundary reactions) *)

CompleteStoichiometry[M_MetabolicNetwork, B_BoundaryNetwork] := 
	Block[{S, Sb, nr, nm, nbm, nbr},
		S = Stoichiometry[M];
		Sb = Stoichiometry[B];
		{nm, nr} = Dimensions[S];
		If[Length[Sb] > 0,
			{nbm, nbr} = Dimensions[Sb];
			ArrayFlatten[{{Join[S, SparseArray[{}, {nbm - nm, nr}]], Sb}}],
			(* else return S itself *)
			S]]


(* A net flux state, completely determined by the net fluxes. *)

Format[NetFluxState[v_?VectorQ]] :=
	"<NetFluxState, " <> ToString[Length[v]] <>" reactions>"

NetFlux[NetFluxState[v_?VectorQ]] := v


(* Conversion from FluxState (see below) *)

NetFluxState[f_FluxState] :=
	NetFluxState[ForwardFlux[f] - ReverseFlux[f]]


(* Determines the net boundary exchange for a given net flux state.
   NOTE: This function does not verify that mass balance constraints are satisfied! See MassBalanceQ *)

MassBalance[M_MetabolicNetwork, f_FluxState] := 
	MassBalance[M, NetFluxState[f]]

MassBalance[M_MetabolicNetwork, f_NetFluxState] := 
	Stoichiometry[M].NetFlux[f]

MassBalance[M_MetabolicNetwork, f_NetFluxState, m_Metabolite] := 
	Stoichiometry[M, m].NetFlux[f]

MassBalance[M_MetabolicNetwork, f_NetFluxState, m:{_Metabolite..}] := 
	Stoichiometry[M, m].NetFlux[f]


(* Check if mass balance is satisfied, up to machine precision *)

MassBalanceQ[M_MetabolicNetwork, f:(_NetFluxState|_FluxState)] := 
	Chop[MassBalanceNorm[M,f]] == 0

(* Norm of mass balance differences for internal metabolites *)

MassBalanceNorm[M_MetabolicNetwork, f:(_NetFluxState|_FluxState)] := 
	Norm[Pick[MassBalance[M,f], Thread[ExchangeStatus[Metabolites[M]]=="Internal"]]]

(* A flux state is completely determined by all forward and reverse reactions
   together with boundary input and output. An alternative representation is
   a net flux state plus exchange fluxes. *)

Format[FluxState[vfwd_?VectorQ, vrev_?VectorQ, bin_?VectorQ, bout_?VectorQ]] :=
	"<FluxState, " <> ToString[Length[vfwd]] <>" reactions, " <> 
		ToString[Length[bin] + Length[bout]] <>" b. exch>"


(* Given forward and reverse fluxes, calculate matching boundary exchange vectors.
   Note that boundary exchange is not uniquely determined in cases where both uptake and release
   exists for a given metabolite; in such cases this function chooses either uptake or release to be zero *)

FluxState[M_MetabolicNetwork, vfwd_?VectorQ, vrev_?VectorQ] :=
	Block[{b},
		(* net boundary exchange *)
		b = MassBalance[M, NetFluxState[vfwd - vrev], ExchangeMetabolites[M]];
		(* Choose positive and negative parts *)
		FluxState[vfwd, vrev, 
			Pick[NegativePart[b], ExchangeMetabolites[M], 
				Metabolite[_, _, "Uptake" | "Free"]],
			Pick[PositivePart[b], ExchangeMetabolites[M],
				Metabolite[_, _, "Release" | "Free"]]]]

(* Create from a net flux state, setting all exchange fluxes to zero *)

FluxState[M_MetabolicNetwork, NetFluxState[v_?VectorQ]] := 
	FluxState[M, PositivePart[v], NegativePart[v]]


(* Arithmetic on flux states *)

FluxState /: Plus[fl:_FluxState..] := FluxState @@ Plus @@ Map[(List@@#&), {fl}]

FluxState /: Times[f_FluxState, x_] := FluxState @@ (x * List @@ f)

FluxState /: Times[x_, f_FluxState] := FluxState @@ (x * List @@ f)


(* This vector corresponds to the InputMetabolites[] list of metabolites *)

BoundaryUptake[FluxState[vfwd_?VectorQ, vrev_?VectorQ, bin_?VectorQ, bout_?VectorQ]] :=
	bin

(* This vector corresponds to the OutputMetabolites[] list of metabolites *)

BoundaryRelease[FluxState[vfwd_?VectorQ, vrev_?VectorQ, bin_?VectorQ, bout_?VectorQ]] :=
	bout


ForwardFlux[FluxState[vfwd_?VectorQ, vrev_?VectorQ, bin_?VectorQ, bout_?VectorQ]] :=
	vfwd

ForwardFlux[M_MetabolicNetwork, f_FluxState, r_Reaction] :=
	Pick[ForwardFlux[f], Reactions[M], r]

ReverseFlux[FluxState[vfwd_?VectorQ, vrev_?VectorQ, bin_?VectorQ, bout_?VectorQ]] :=
	vrev

NetFlux[f_FluxState] := ForwardFlux[f] - ReverseFlux[f]


ExchangeFlux[f_FluxState] := Min /@ Transpose[{ForwardFlux[f], ReverseFlux[f]}]

BoundaryFlux[f_FluxState] := Join[BoundaryUptake[f], BoundaryRelease[f]]


(* This is the scaled or "compactified" exchange flux, as defined by Wiechert.
   In many cases the parameter \[Beta] is set to 1, which is suitable for fluxes of this magnitude *)

ExchangeCoefficient[f_FluxState, \[Beta]_?Positive] := 
	Block[{vex},
		vex = ExchangeFlux[f];
		vex / (\[Beta] + vex)]

ExchangeCoefficient[f_FluxState] := 
	ExchangeCoefficient[f, 1.0]


(* The total flux through a metabolite pool, assuming mass balance.
   This quantity must be positive for a metabolite to evaluate its isotopomer distribution. *)

MetaboliteTotalFlux[M_MetabolicNetwork, f_FluxState] :=
	PositivePart[Stoichiometry[M]].ForwardFlux[f] + 
	NegativePart[Stoichiometry[M]].ReverseFlux[f] +
	Take[
		PositivePart[Stoichiometry[BoundaryNetwork[M]]].BoundaryFlux[f],
		Length[Metabolites[M]]]


End[];

EndPackage[];
