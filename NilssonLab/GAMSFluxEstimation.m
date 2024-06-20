(* This package allows export of flux analysis models to the GAMS optimization format*)

BeginPackage[
	"NilssonLab`GAMSFluxEstimation`",
	{
		"NilssonLab`Utilities`",
		"NilssonLab`MetabolicNetwork`Common`",
		"NilssonLab`AtomMap`Common`",
		"NilssonLab`EMUTracing`"
	}
];


(* Public symbols / usage strings *)

GAMSExport;
GAMSExportRatioCI;
ImportGAMSFluxes;
ImportGAMSIsotopes;
ImportGAMSEmuMixing;
ImportGAMSFluxCI;
ImportGAMSFluxRatioCI;
GAMSExportRatioBalanceProblem;


Begin["NilssonLab`GAMSFluxEstimation`Private`"];


(*Helper functions*)


(* Remove special characters from EMU ID strings and encodes atom numbers as 1.. 9, a, b ...  *)

(*For EMUs with multiple equivalent fragments, we here use only the first fragment as identifier.*)


emuString[EMU[m_Metabolite, {al:{{_Integer...}..}, ___}]] := 
	MetaboliteShortName[m] <> "_" <> ListToString[Table[
		ListToString[
			a /. {i_Integer /; i >= 10 :> FromCharacterCode[97 + i - 10]}, ""],
		{a, al}], "_"]


massIndex[sym_String, midDepth_Integer] := 
	ListToString[Table[sym <> ToString[k], {k,midDepth}], ","]


massIndexVector[ix:{_Integer..}] := 
	ListToString[Table["'" <> ToString[i] <> "'", {i, ix}], ","]


massFunctionIndex[sym_String, k_Integer] := 
		"mass" <> ToString[k] <> "(" <> sym <> ToString[k] <> ")"


massIndexConstraint[emuSym_String, massSym_String, midDepth_Integer] := 
	ListToString[Table[
		massFunctionIndex[massSym, k] <> " <= emudim"<>ToString[k]<>"(" <> emuSym <> ")",
		{k, midDepth}], " and "]


massIndexConstraintCond[midDepth_Integer] := 
	ListToString[Table[
		massFunctionIndex["kpp", k] <> " <= emudim" <> ToString[k] <> "(e) " <>
			"+ emudim" <> ToString[k] <> "(ep)",
		{k, midDepth}], " and "]


(* Replace elements `n` in a list of strings using StringTemplate *)

formatLines[lines:{_String..}, repl_List] :=
	Table[StringTemplate[line]@@repl, {line, lines}]


(* This defines the order of fluxes used in the GAMS representation.
 * It is equal to the order used by EMUFluxMap except that release fluxes are added at the end.
 *)

gamsMetaboliteList[B_BoundaryNetwork] :=
	MetaboliteShortName /@ Metabolites[B];

gamsReactionList[M_MetabolicNetwork, B_BoundaryNetwork] :=
	Block[{nu, nr},
		(* list of uptake/fwd/rev/release fluxes *)
		nu = Length[InputMetabolites[M]];
		nr = Length[OutputMetabolites[M]];
		Join[
			Take[BoundaryFluxNames[B], nu],
			Table[StringJoin[ReactionID[r], "_f"], {r, Reactions[M]}],
			Table[StringJoin[ReactionID[r], "_r"], {r, Select[Reactions[M], ReversibleQ]}],
			Take[BoundaryFluxNames[B], -nr]]]

gamsStoichiometry[M_MetabolicNetwork, B_BoundaryNetwork] :=
	Block[{nu, nr, revi, S, Sb},
		(* list of uptake/fwd/rev/release fluxes *)
		nu = Length[InputMetabolites[M]];
		nr = Length[OutputMetabolites[M]];
		revi = Flatten[Position[ReversibleQ[Reactions[M]],True]];
		S = Stoichiometry[M];
		Sb = Stoichiometry[B];
		SparseArray[ArrayFlatten[{{
			Take[Sb,All,nu],
			Join[S, SparseArray[{}, {Length[Sb] - Length[S], Length[Reactions[M]]}]],
			Join[-S[[All, revi]], SparseArray[{}, {Length[Sb] - Length[S], Length[revi]}]],
			Take[Sb,All,-nr]}}]]]


(*
 * Creating GAMS Files
 *)


(*This adds a comment block and GAMS parameters*)

MakeHeader[] :=
	{"**",
	 "**  GAMS model exported from Mathematica",
	 "**",
	 (* ignore digits beyond GAMS machine precision *)
	 "$offdigit"}

(*Stoichiometry block*)


(*
 * Here we define the metabolite and reaction lists and stoichiometry matrix.
 * The full matrix including uptake and release fluxes as well as input and output metabolites
 * is used to allow placing constraints on uptake and release fluxes.
 * The subset of metabolites excluding input and output is used for the stoichiometry constraint
 * Sv = 0 later on.
 *)

MakeStoichiometryBlock[M_MetabolicNetwork, B_BoundaryNetwork] :=
	Block[{ml, rl, S, Sij, Sx},
		ml = gamsMetaboliteList[B];
		(* list of uptake/fwd/rev/release fluxes *)
		rl = gamsReactionList[M, B];
		(* build stoichiometry matrix over fwd/rev/boundary fluxes *)
		S = gamsStoichiometry[M, B];
		Sij = NonZeroIndices[S];
		Sx = NonZeroElements[S];
		Join[
			{"** MakeStoichiometryBlock[] **",
			 "Set i \"metabolites\" / " <> ListToString[ml,", "] <> " / ;",
			 "Set massbal(i) \"metabolites for mass balance\" / " 
				<> ListToString[MetaboliteShortName /@ Metabolites[M],", "] <> " / ;",
			 "Alias (i,ir) ;",
			 "Set j \"reactions\" / " <> ListToString[rl,", "] <> " / ;",
			 "Alias (j,jp,jr) ;",
			 "Parameter S(i,j) \"stoichiometry matrix\" ;"},
			Table[
				"  S('" <> ml[[Sij[[k,1]]]] <> "','" <> rl[[Sij[[k,2]]]] <> "') = " 
				 <> ToString[Sx[[k]]] <>" ;",
				{k, Length[Sij]}]]]


(*This defines the list of EMUs*)

MakeEMUBlock[emuEq:{_EMUEquation..}] :=
	Block[{ei,es, midDepth, maxSizes},
		ei = Flatten[InternalEMUs /@ emuEq];
		es = Flatten[SubstrateEMUs /@ emuEq];
		midDepth = Length[Dimensions[First[ei]]];
		maxSizes = Last /@ Sort /@ Transpose[Dimensions /@ Join[ei,es]];
		Join[
			{"** MakeEMUBlock[] **",
			 "Set e \"EMUs\" / " <>
				ListToString[emuString /@ Join[es,ei], ", "] <> " / ;",
			 "Set intemu(e) \"internal EMUs\" / " <>
				ListToString[emuString /@ ei, ", "] <> " / ;",
			 "Set subemu(e) \"substrate EMUs\" / " <>
				ListToString[emuString /@ es, ", "] <> " / ;",
			 "Alias(e, ep, epp) ;"},
			Flatten[Table[
				{"Set k" <> ToString[k] <> " MIDs / 0*" <> ToString[maxSizes[[k]]] <> " / ;",
				 "Alias(k" <> ToString[k] <> ", kp" <> ToString[k] <> ", kpp" <> ToString[k] <> ") ",
				 "Parameter "<> massFunctionIndex["k", k] <> " / " <> 
					ListToString[Table[
						ToString[i]<>" "<>ToString[i], {i,0,maxSizes[[k]]}], ", "] <> " / ;",
				"Parameter emudim" <> ToString[k] <> "(e) / "<>
				ListToString[Table[
					emuString[e]<>" "<>ToString[Dimensions[e][[k]]], {e, Join[es,ei]}], ", "]
					<> " / ;"},
				{k, midDepth}]]]]


(* Measurements block (data)
 * This requires a standard deviation (normal noise) for each MID.
 *)

MakeMeasurementsBlock[
		measFluxNames:{_String..}, measReaction:{_String..},
		measFluxMean_?VectorQ, measFluxSD_?VectorQ, fluxMixing_SparseArray,
		expName:{_String..}, measName:{_String..}, measMID_List, midsd_List,
		measEMU:{_EMU..}, emuMixing_SparseArray] :=
	Join[
		{"** MakeMeasurementsBlock[] **"},
		fluxMeasurementsBlock[measFluxNames, measReaction, measFluxMean, measFluxSD, fluxMixing],
		midMeasurementsBlock[expName, measName, measMID, midsd, measEMU, emuMixing]]


(* Using measured fluxes only: *)

MakeMeasurementsBlock[measFluxNames:{_String..}, measReaction:{_String..},
					measFluxMean_?VectorQ, measFluxSD_?VectorQ, fluxMixing_SparseArray] :=
	Join[
		{"** MakeMeasurementsBlock[] **"},
		fluxMeasurementsBlock[measFluxNames, measReaction, measFluxMean, measFluxSD, fluxMixing]]


(* Using MID data only (fluxes must then be constrained by a "sink" equation") *)

MakeMeasurementsBlock[
			expName:{_String..}, measName:{_String..}, measMID_List, midsd_List,
			measEMU:{_EMU..}, emuMixing_SparseArray] :=
	Join[
		{"** MakeMeasurementsBlock[] **"},
		midMeasurementsBlock[expName, measName, measMID, midsd, measEMU, emuMixing]]


(*Part of measurement block for fluxes*)

fluxMeasurementsBlock[measFluxNames:{_String..}, measReaction:{_String..},
					measFluxMean_?VectorQ, measFluxSD_?VectorQ, fluxMixing_SparseArray] :=
	Block[{fluxMixingIndex, fluxMixingCoeff},
		Join[
			(* matrix of measurement-flux coefficients *)
			fluxMixingIndex = NonZeroIndices[fluxMixing];
			fluxMixingCoeff = NonZeroElements[fluxMixing];	
			{"Set msFlux \"Measured Fluxes\" / " <> ListToString[measFluxNames,", "] <> " /"},
			{"Parameter measFluxReact(msFlux, j) \"measurement-to-flux map\" ;"},
			Table[
				"  measFluxReact('" <> measFluxNames[[ fluxMixingIndex[[k,1]] ]] <> "','" <> 
					measReaction[[fluxMixingIndex[[k,2]]]] <> "') = " <>
					ToString[CForm[fluxMixingCoeff[[k]]]] <>" ;",
				{k, Length[fluxMixingIndex]}],
			{"Parameter vobs(msFlux) \"mean of measured flux\" ;"},
			Table[
				"  vobs('" <> measFluxNames[[i]] <> "') = " <> ToString[CForm[measFluxMean[[i]]]] <>" ;",
				{i, Length[measFluxNames]}],
			{"Parameter vsd(msFlux) \"std.dev for measured flux\" ;"},
			 Table[
					"  vsd('" <> measFluxNames[[i]] <> "') = " <> ToString[CForm[measFluxSD[[i]]]] <>" ;",
					{i, Length[measFluxNames]}]]]


(* Part of measurement block for MIDs *)

midMeasurementsBlock[
			expName:{_String..}, measName:{_String..}, measMID_List, midsd_List,
			measEMU:{_EMU..}, emuMixing_SparseArray] :=
	Block[{midDepth, measDim, EMij, EMx},
		midDepth = Length[Dimensions[First[measEMU]]];
		(* measurement array dimensions is given by the emuMixing matrix *)
		measDim = Table[
			Dimensions[First[Pick[measEMU, em, _?Positive]]],
			{em, emuMixing}];
		EMij = NonZeroIndices[emuMixing];
		EMx = NonZeroElements[emuMixing];	
		Join[
			{"Set expmt \"Tracer experiments\" / " <> ListToString[expName,", "] <> " /"},
			{"Set ms \"Measured MIDs\" / " <> ListToString[measName,", "] <> " /"},
			{"Set msEmu(ms,e) \"Measurement-emu pairs\" / " <> ListToString[Table[
				measName[[ EMij[[k,1]] ]] <> "." <> emuString[measEMU[[EMij[[k,2]]]] ],
				{k, Length[EMij]}], ", "] <> " /"},
			{"Positive Variable measEmu(ms,e) \"EMU measurement map\" ;"},
			Table[
				"  measEmu.l('" <> measName[[ EMij[[k,1]] ]] <> "','" <> emuString[measEMU[[EMij[[k,2]]]] ] <> "') = " 
					<> ToString[CForm[EMx[[k]]]] <>" ;",
				{k, Length[EMij]}],
			{"Parameter midsd(expmt, ms) \"std.dev for measured MIDs\" ;"},
			Flatten[Table[
				"  midsd('" <> expName[[k]] <> "', '" <> measName[[i]] <> "') = " <> ToString[CForm[midsd[[k, i]]]] <>" ;",
				{k, Length[expName]}, {i, Length[measName]}]],
			{"Parameter xobs(expmt, ms," <> massIndex["k", midDepth] <> ") \"measured MIDs\" ;"},
			Flatten[Table[
				"  xobs('" <> expName[[k]] <> "', '" <> measName[[i]] <> "', " <>  massIndexVector[t-1] <>  ") = " 
					<> ToString[CForm[measMID[[k, i, Sequence @@ t]]]] <> " ;",
				{k, Length[expName]}, {i, Length[measName]},
				{t, Tuples[Map[Range[#+1]&, measDim[[i]] ]]} ]]]]


(* This creats the EMU reaction network *)

MakeEMUReactionsBlock[er:{_EMUReaction...}, rl:{_String..}] :=
	Join[
		{"** MakeEMUReactionsBlock[] **",
		 "Parameter er(e,ep,j) emu reactions ;"},
		Cases[er, EMUReaction[j_Integer, e1_EMU, e2_EMU, c_] :> 
			"er('" <> emuString[e1] <> " ', '" <> emuString[e2] <> "', '" <>
			rl[[j]] <> "') = " <> ToString[N[c]] <> " ;"]]


MakeEMUCondensationsBlock[er:{_EMUReaction...}, emuEq:{_EMUEquation..}, rl:{_String..}] :=
	Block[{cond},
		cond = Flatten[CondensationEMUs/@emuEq];
		Join[
		{"** MakeEMUCondensationsBlock[] **",
		 (* list of condensations *)
		 "Set cond(e,ep) \"condensations\" / " <>
			ListToString[
				cond /. Condense[e1_,e2_] :> 
				("'" <> emuString[e1] <> "'.'" <> emuString[e2]) <> "'", ", "]
			<> " / ;"},
		(* EMU reactions for condensations *)
		{"Parameter ecr(e,ep,epp,j) \"emu condensation reactions\" ;"},
		Cases[er, EMUReaction[j_Integer, Condense[e1_EMU, e2_EMU], e3_EMU, c_] :> 
			"  ecr('" <> emuString[e1] <> "', '" <> emuString[e2] <> "', '" <>
			emuString[e3] <> "', '" <> rl[[j]] <> "') = " <> ToString[N[c]] <> " ;"]]]


MakeVariablesBlock[emuEq:{_EMUEquation..},
				expName:{_String..}, substrateMIDs:{{_Rule..}..}, vmin_, vmax_] :=
	Block[{es, mid, midDepth},
		es = Flatten[SubstrateEMUs /@ emuEq];
		midDepth = Length[Dimensions[First[es]]];
		Join[formatLines[
			{"** MakeVariablesBlock[] **",
			(* declare variables *)
			 "Variables",
			 "  v(j) \"flux through reaction j\"",
			 "  f \"objective function\" ",
			 "  vopt(j) \"flux through reaction j at optimal solution\"",
			 "  rn(j) \"numerator coefficients for range problem\"",
			 "  rd(j) \"denominator coefficients for range problem\"",
			 "  rc \"denominator constant for range problem\"",
			 "  x(expmt, e, `1`) \"MID for EMU e and mass k\"",
			 "  z(expmt, e,ep, `2`) \"MID for condensation EMU e x ep and mass kpp\"",
			 "  xfit(expmt, ms, `1`) \"Fit to measured MIDs\"",
			 "  xopt(expmt, e, `1`) \"MID for EMUs at optimum\"",
			 "  zopt(expmt, e,ep, `2`) \"MID for Condensation EMUs at optimum\"",
			 "  ci(i, j) \"coefficients for range problem\"",
			 "  fr \"range objective function\"",
			 "  fci \"flux ratio CI objective function\"",
			 "  fbound \"objective function bound\" ;",
			 "Positive Variables v, x, z ;",
			(* bounds on fluxes *)
			 "v.lo(j) = `3` ;",
			 "v.up(j) = `4` ;"},
			 {massIndex["k", midDepth], massIndex["kpp", midDepth], vmin, vmax}],
			(* fix MIDs for substrate EMUs *)
			Flatten[Table[
				mid = es /. substrateMIDs[[k]];
				Table[
					"x.fx('" <> expName[[k]] <> "', '" <> 
					emuString[es[[i]]] <> "'," <>  massIndexVector[t-1] <>  ") = " <>
					ToString[CForm[mid[[i, Sequence @@ t]]]] <> " ;",
					{i, Length[es]},
					{t, Tuples[Map[Range[#+1]&, Dimensions[es[[i]]]]]}],
				{k, Length[expName]}]]]]


(* Variables Without isotopes: *)

MakeVariablesBlock[vmin_, vmax_] :=
	{"** MakeVariablesBlock[] **",
	(* declare variables *)
	 "Variables",
	 "  v(j) \"flux through reaction j\"",
	 "  f \"objective function\" ",
	 "  vopt(j) \"flux through reaction j at optimal solution\"",
	 "  rn(j) \"numerator coefficients for range problem\"",
	 "  rd(j) \"denominator coefficients for range problem\"",
	 "  rc \"denominator constant for range problem\"",
	 "  fr \"range objective function\"",
	 "  fbound \"objective function bound\" ;",
	 "Positive Variables v ;",
	(* bounds on fluxes *)
	 "v.lo(j) = " <> ToString[vmin] <> " ;",
	 "v.up(j) = " <> ToString[vmax] <> " ;"}


(* Equations block for MID data and measured fluxes *)

MakeEquationsBlock[midDepth_Integer] :=
	Join[
		{"** MakeEquationsBlock[] **",
		 "Equations"},
		basicEquationsBlock[midDepth],
		midFluxObjectiveEquationsBlock[midDepth],
		midBalanceEquationsBlock[midDepth]]


(* Equations without MID data, stoichiometry constraints only *)

MakeEquationsBlock[] :=
	Join[
		{"** MakeEquationsBlock[] **",
		 "Equations"},
		basicEquationsBlock[],
		fluxObjectiveEquationsBlock[]]		


(* Equations With isotopes and sink equation *)

MakeEquationsBlock[midDepth_Integer, fixedReactions:{_String...}] :=
	Join[
		{"** MakeEquationsBlock[] **",
		 "Equations"},
		basicEquationsBlock[midDepth],
		midObjectiveEquationsBlock[midDepth],
		midBalanceEquationsBlock[midDepth],
		fixedFluxesEquationsBlock[fixedReactions, 100]]


(* Declaration of symbols and stoichiometry equations *)

basicEquationsBlock[midDepth_Integer] :=
	formatLines[
		{"midmix \"mix compartment-specific MIDs into fitted metabolite MIDs\"",
		"objective \"defines the L2 error objective function f\"",
		"rangeobj \"objective function for confidence intervals\"",
		"massbaleq(i) \"mass balance equation\"",
		"objbound \"objective bound for range problem\"",
		"midsum(expmt,e) \"the MIDs for each EMU e must sum to 1\"",
		"mesum(ms) \"the measurement-EMU mixing coefficients sum to 1\"",
		"emubalance(expmt,e,`1`) \"mass balance for EMUs\"",
		"condense(expmt,e,ep,`2`) \"convolution of MIDs\" ;",
		"midmix(expmt, ms, `1`) ..",
		"  xfit(expmt, ms, `1`) =e= sum((e)$msEmu(ms,e), measEmu(ms,e)*x(expmt, e, `1`)) ; ",
		"rangeobj .. ",
		"  fr =e= sum(j, rn(j)*v(j)) / (sum(j, rd(j)*v(j)) + rc) ; ",
		"objbound .. ",
		"  fbound =g= f ;",
		"massbaleq(massbal(i)) .. ",
		"  sum(j, S(i,j)*v(j)) =e= 0 ;"}, {massIndex["k", midDepth], massIndex["kpp", midDepth]}]


(* ::Text:: *)
(*Without isotopes*)


basicEquationsBlock[] :=
	{
		"objective \"defines the L2 error objective function f\"",
		"rangeobj \"objective function for the range problem\"",
		"massbaleq(i) \"mass balance equation\"",
		"objbound \"objective bound for range problem\" ;",
		"rangeobj .. ",
		"  fr =e= sum(j, rn(j)*v(j)) / (sum(j, rd(j)*v(j)) + rc) ; ",
		"objbound .. ",
		"  fbound =g= f ;",
		"massbaleq(massbal(i)) .. ",
		"  sum(j, S(i,j)*v(j)) =e= 0 ;"}


fluxObjectiveEquationsBlock[] :=
	{"objective .. ",
	 "  f =e= sum( msFlux, sqr((sum((j)$measFluxReact(msFlux, j), measFluxReact(msFlux, j)*v(j)) - vobs(msFlux)) / vsd(msFlux))) ;"}


(* ::Text:: *)
(*MID and flux objective*)
(*Fit measured MIDs and measured fluxes*)


midFluxObjectiveEquationsBlock[midDepth_Integer] :=
	formatLines[
		{"objective .. ",
		(* measured fluxes *)
		"  f =e= sum( msFlux, sqr((sum((j)$measFluxReact(msFlux, j), measFluxReact(msFlux, j)*v(j)) - vobs(msFlux)) / vsd(msFlux))) + ",
		(* measured MIDs *)
		"    sum((expmt, ms, `1`)$xobs(expmt, ms, `1`), sqr((xfit(expmt, ms, `1`) - xobs(expmt, ms, `1`)) / midsd(expmt, ms)) ) ;"},
		{massIndex["k", midDepth]}]


(* ::Text:: *)
(*MID objective*)


midObjectiveEquationsBlock[midDepth_Integer] :=
	{"objective .. ",
	 "  f =e= sum((expmt, ms," <> massIndex["k", midDepth] <>")$xobs(expmt, ms," 
		<> massIndex["k", midDepth] <>"), " <>
	 "sqr( (sum((e)$msEmu(ms,e), measEmu(ms,e)*x(expmt, e," <> massIndex["k", midDepth] <>"))"
		<> " - xobs(expmt, ms," <> massIndex["k", midDepth] <>")) / midsd(expmt,ms)) ) ;"}		 


(* ::Text:: *)
(*Equation requiring a set of fluxes to sum to a constant. Note that this may be incompatible with measured reactions!*)


fixedFluxesEquationsBlock[fixedReactions:{_String..}, total_] :=
	{"fixedfluxes .. " <> ListToString[
		Table["v('" <> fixedReactions[[i]] <> "') ", {i, Length[fixedReactions]}], "+ "]
		<> "=e= " <> ToString[total] <> " ;"}


(* ::Text:: *)
(*Dummy for case with no fixed fluxes*)


fixedFluxesEquationsBlock[{}, total_] :=
	{"fixedfluxes .. 1 =e= 1 ;"}


(* ::Text:: *)
(*MID balance equations*)


midBalanceEquationsBlock[midDepth_Integer] :=
	{"midsum(expmt,e) .. ",
	 "  sum((" <> massIndex["k", midDepth] <>")$(" <> massIndexConstraint["e", "k", midDepth] <> 
		"), x(expmt,e," <> massIndex["k", midDepth] <>")) =e= 1 ;",
	 "mesum(ms) .. ",
	 "  sum(e$msEmu(ms,e), measEmu(ms,e)) =e= 1 ;",
	 "emubalance(expmt, intemu(e)," <> massIndex["k", midDepth] <>")$(" <> 
		massIndexConstraint["e", "k", midDepth] <> ") .. ",
	 "  sum((ep, j)$er(ep, e, j), er(ep, e, j)*v(j)*x(expmt, ep," <> massIndex["k", midDepth] <>")) + ",
	 "  sum((ep, epp, j)$ecr(ep, epp, e, j), ecr(ep, epp, e, j)*v(j)*z(expmt, ep,epp,"
		<> massIndex["k", midDepth] <>"))",
	 "  =e= ( sum((ep,j)$er(ep, e, j), er(ep, e, j)*v(j)) + ",
	 "  sum((ep, epp, j)$ecr(ep, epp, e, j), ecr(ep, epp, e, j)*v(j)) )"<>
		" *x(expmt, e," <> massIndex["k", midDepth] <>") ; ",
	 "condense(expmt, e,ep," <> massIndex["kpp", midDepth] <>")$(cond(e, ep) and (" <>
		massIndexConstraintCond[midDepth] <> ")) .. ",
	 "  z(expmt, e,ep," <> massIndex["kpp", midDepth] <>") =e= sum((" <>
			massIndex["k", midDepth] <> "," <> massIndex["kp", midDepth] <>")$(" <>
			ListToString[Table[
				massFunctionIndex["kpp", k] <> " = " <> massFunctionIndex["k", k] <>
				" + " <> massFunctionIndex["kp", k],
				{k, midDepth}], " and "] <>
			"), x(expmt, e," <> massIndex["k", midDepth] <>")*x(expmt, ep," <> massIndex["kp", midDepth] <>")) ;"}


(* ::Subsection::Closed:: *)
(*Initial point block*)


MakeInitialPointBlock[emuEq:{_EMUEquation..}, emap_EMUFluxMap, initFlux_FluxState,
			expName:{_String...}, fluxSol:{_EMUSolution..}, rl:{_String..}] :=
	Block[{v0, ei, ec, X0, Z0},
		v0 = Join[emap[initFlux], BoundaryRelease[initFlux]];
		ei = Flatten[InternalEMUs /@ emuEq];
		ec = Flatten[CondensationEMUs /@ emuEq];
		Join[
			{"** MakeInitialPointBlock[] **"},
			Table[
				"v.l('" <> rl[[j]] <> "') = " <> ToString[CForm[v0[[j]]]] <> " ;",
				{j, Length[v0]}],
			Flatten[Table[
				(* loop over experiments *)
				X0 = fluxSol[[k]] /@ ei;
				Z0 = fluxSol[[k]] /@ ec;
				{Table[
					"x.l('" <> expName[[k]] <> "', '" <> 
						emuString[ei[[i]]] <> "'," <>  massIndexVector[t-1] <>  ") = "  <>
						ToString[CForm[X0[[i, Sequence @@ t]]]] <> " ;",
					{i, Length[ei]},
					{t, Tuples[Map[Range[#+1]&, Dimensions[ei[[i]]]]]}],
				 Table[
					"z.l('" <> expName[[k]] <> "', '"  <>
						emuString[ec[[i,1]]] <>  "','" <> emuString[ec[[i,2]]] <>  "'," <>
						massIndexVector[t-1] <>  ") = " <>
						ToString[CForm[Z0[[i, Sequence @@ t]]]] <> " ;",
					{i, Length[ec]},
					{t, Tuples[Map[Range[#+1]&, Dimensions[ec[[i]]]]]}]},
				{k, Length[expName]}]]]]		


(* ::Text:: *)
(*Without isotopes.*)


MakeInitialPointBlock[emap_EMUFluxMap, initFlux_FluxState, rl:{_String..}] :=
	Block[{v0},
		v0 = Join[emap[initFlux], BoundaryRelease[initFlux]];
		Join[
			{"** MakeInitialPointBlock[] **"},
			Table[
				"v.l('" <> rl[[j]] <> "') = " <> ToString[CForm[v0[[j]]]] <> " ;",
				{j, Length[v0]}]]]


(* ::Subsection::Closed:: *)
(*Models block*)


(* ::Text:: *)
(*This defines the models used*)


(* ::Text:: *)
(*Using isotopes*)


MakeIsotopeModelsBlock[] :=
	{"** MakeModelsBlock[] ** ",
	 "Model",
	 "  emuopt / midmix, objective, massbaleq, mesum, midsum, emubalance, condense /",
	 "  fluxci / rangeobj, objbound, midmix, objective, massbaleq, mesum, midsum, emubalance, condense / ;"}


(* ::Text:: *)
(*Without isotopes. Should rename the GAMS problems here.*)


MakeFluxModelsBlock[] :=
	{"** MakeModelsBlock[] ** ",
	 "Model",
	 "  emuopt / objective, massbaleq /",
	 "  fluxci / rangeobj, objbound, objective, massbaleq / ;"}


(* ::Text:: *)
(*With isotopes, without flux data (sink equaton), including ratio CI problem*)


MakeSinkModelsBlock[] :=
	{"** MakeModelsBlock[] ** ",
	 "Model",
	 "  emuopt / objective, massbaleq, mesum, midsum, emubalance, condense, fixedfluxes /",
	 "  fluxci / rangeobj, objbound, objective, massbaleq, mesum, midsum, emubalance, condense, fixedfluxes /",
	 "  fluxratioci / ratiociobj, objbound, objective, massbaleq, mesum, midsum, emubalance, condense, fixedfluxes / ;"}


(* ::Subsection::Closed:: *)
(*Optimization block*)


(* ::Text:: *)
(*This is the code for solving the flux optimization problem and writing solutions to file.*)


miIndexList[midDepth_Integer] :=
	ListToString[Table["k" <> ToString[k] <> ".tl", {k, midDepth}], ","]


MakeOptimizationBlock[modelName_String, midDepth_Integer] :=
	formatLines[
	{"** MakeOptimizationBlock[] **",
	 "Display 'optimizing fluxes ...' ;",
	 "Solve emuopt using NLP minimizing f ;",
	 "** Save optimum values **",
	 "vopt.l(j) = v.l(j) ;",
	 "xopt.l(expmt, e, `1`) = x.l(expmt, e,`1`) ;",
	 "zopt.l(expmt, e, ep, `2`) =  z.l(expmt, e, ep, `2`) ;",
	 "** Write objective value **",
	 "file objval /" <> modelName <> "_objval.tsv/ ;",
	 "objval.pc = 6 ;",
	 "put objval ; put f.l:0:10/ ; putclose objval;",
	 "** Write optimal fluxes **",
	 "file fluxopt /" <> modelName <> "_fluxopt.tsv/ ;",
	 "fluxopt.pc = 6 ;",
	 "put fluxopt ; loop(j, put j.tl, v.l(j):0:10/); putclose fluxopt;",
	 "** Write optimal MIDs **",
	 "file isotopes /" <> modelName <> "_isotopes.tsv/ ;",
	 "isotopes.pc = 6 ;",
	 "put isotopes ; loop((expmt, e, `1`)$(" 
		<> massIndexConstraint["e", "k", midDepth] <> "), "
		<> "put expmt.tl, e.tl, " <> miIndexList[midDepth] <> ","
		<> "x.l(expmt,e, `1`):0:10/); "
		<> "putclose isotopes;",
	 "** Write EMU mixing coefficients **",
	 "file emumixing /" <> modelName <> "_emumixing.tsv/ ;",
	 "emumixing.pc = 6 ;",
	 "put emumixing ; loop((ms, e)$msEmu(ms,e), "
		<> "put ms.tl, e.tl, measEmu.l(ms,e):0:10/); "
		<> "putclose emumixing;",
	 "** Write fitted MIDs **",
	 "file fitmids /" <> modelName <> "_fitmids.tsv/ ;",
	 "fitmids.pc = 6 ;",
	 "put fitmids ; loop((expmt, ms, `1`)$xobs(expmt, ms, `1`), "
		<> "put expmt.tl, ms.tl, " <> miIndexList[midDepth] <> ","
		<> "xfit.l(expmt,ms,`1`):0:10/); "
		<> "putclose fitmids;",
	 "** Write marginal MID values (derivatives) **",
	 "file midmarginal /" <> modelName <> "_midmarginal.tsv/ ;",
	 "midmarginal.pc = 6 ;",
	 "put midmarginal ; loop((expmt, ms, `1`)$xobs(expmt, ms, `1`), "
		<> "put expmt.tl, ms.tl, " <> miIndexList[midDepth] <> ","
		<> "midmix.m(expmt, ms, `1`):0:10/); "
		<> "putclose midmarginal;"
	},
	{massIndex["k", midDepth], massIndex["kpp", midDepth]}]


(* ::Text:: *)
(*Without isotopes*)


MakeOptimizationBlock[modelName_String] :=
	{"** MakeOptimizationBlock[] **",
	 "Display 'optimizing fluxes ...' ;",
	 "Solve emuopt using NLP minimizing f ;",
	 "vopt.l(j) = v.l(j) ;",
	 "file fluxopt /" <> modelName <> "_fluxopt.tsv/ ;",
	 "fluxopt.pc = 6 ;",
	 "put fluxopt ; loop(j, put j.tl, v.l(j):0:10/); putclose fluxopt;" }


(* ::Subsection::Closed:: *)
(*Flux confidence intervals block*)


(* ::Text:: *)
(*This uses a loop to solve the flux range optimization problem for an indicated set of reactions,*)
(*and writes solutions to file.*)
(*The parameter chisq controls the optimum bound; this is a chi-square distribution quantile, for examle F(0.95) for a 5% CI*)


(* ::Text:: *)
(*This version from a list of Reaction, generates CI objectives for all net fluxes*)


MakeNetFluxesCIBlock[modelName_String, midDepth_Integer, rl:{_Reaction..}, chisq_] :=
	MakeFluxCIBlock[modelName, midDepth, ReactionID /@ rl,
		(* numerator from reaction coefficients *)
		Table[
			If[ReversibleQ[r],
				{{ReactionID[r] <> "_f", 1}, {ReactionID[r] <> "_r", -1}},
				{{ReactionID[r] <> "_f", 1}}],
			{r, rl}],
		(* denominator is constant 1 *)
		Table[{}, {r, rl}], Table[1, {r, rl}],
		chisq]


(* ::Text:: *)
(*Create range (confidence interval) objectives. We specify a ratio of linear combination for each range objective and an optional constant in the numerator. In the special case of defining a range objective over a single variable X, we put X in the numerator and 1 in the denominator.*)


MakeFluxCIBlock[
		modelName_String, midDepth_Integer,
		rangeNames:{_String..}, rangeNumerator:{{{_String, _Integer}..}..},
		rangeDenominator:{{{_String, _Integer}...}..}, rangeConst:{_?NumericQ..},
		chisq_] :=
	Join[
	{"** MakeFluxCIBlock[] **",
	 "Set roname \"range objective name\" / " <> ListToString[rangeNames, ", "] <>" / ;",
	 "Parameter rangenum(roname, j) \"range numerator coefficients\" ;"},
	Join @@ Table[
		StringTemplate["  rangenum('`1`', '`2`') = `3` ;"][rangeNames[[i]], rn[[1]], rn[[2]]],
		{i, Length[rangeNames]}, {rn, rangeNumerator[[i]]}],
	{"Parameter rangeden(roname, j) \"range denominator coefficients\" ;"},
	Join @@ Table[
		If[Length[rangeDenominator[[i]]] > 0,
			Table[
				StringTemplate["  rangeden('`1`', '`2`') = `3` ;"][rangeNames[[i]], rd[[1]], rd[[2]]],
				{rd, rangeDenominator[[i]]}],
			(* avoid GAMS empty parameter set error by setting an arbitrary element to zero *)
			{StringTemplate["  rangeden('`1`', '`2`') = 0 ;"][rangeNames[[i]], rangeNumerator[[i, 1, 1]]]}],
		{i, Length[rangeNames]}],
	{"Parameter rangeconst(roname) \"range denominator constant\" ;"},
	Table[
		StringTemplate["  rangeconst('`1`') = `2` ;"][rangeNames[[i]], rangeConst[[i]]],
		{i, Length[rangeNames]}],
	{"fbound.fx = f.l + " <> ToString[chisq] <> " ;",
	 (* upper range problems *)
	 "file rangehi /" <> modelName <> "_rangehi.tsv/ ;",
	 "rangehi.pc = 6 ;",
	 "rangehi.pw = 32767 ;",
	 "file fluxhi /" <> modelName <> "_fluxhi.tsv/ ;",
	 "fluxhi.pc = 6 ;",
	 "fluxhi.pw = 32767 ;",
	 (* loop over range problems *)
	 "loop(roname,",
	(* initialize to MFA optimum *)
	 "  v.l(j) = vopt.l(j) ;",
	 If[midDepth > 0,
		StringTemplate["  x.l(expmt, e, `1`) = xopt.l(expmt, e, `1`) ;"][massIndex["k", midDepth]],
		""],
	 (* set objective coefficients *)
	 "  rn.fx(j) = rangenum(roname, j) ;",
	 "  rd.fx(j) = rangeden(roname, j) ;",
	 "  rc.fx = rangeconst(roname) ;",
	 "  Solve fluxci using NLP maximizing fr ;",
	 "  put rangehi ;",
	 "  put roname.tl, fr.l / ;",
	 "  put fluxhi ;",
	 "  put roname.tl ;",
	 "  loop(j, put v.l(j):0:10 ; ) ;", 
	 "  put / ;",
	 ") ;",
	 "putclose rangehi ;",
	 "putclose fluxhi ;",
	(* lower range problems *)
	 "file rangelow /" <> modelName <> "_rangelow.tsv/ ;",
	 "rangelow.pc = 6 ;",
	 "rangelow.pw = 32767 ;",
	 "file fluxlow /" <> modelName <> "_fluxlow.tsv/ ;",
	 "fluxlow.pc = 6 ;",
	 "fluxlow.pw = 32767 ;",
	(* loop over range problems *)
	 "loop(roname,",
	(* initialize to MFA optimum *)
	 "  v.l(j) = vopt.l(j) ;",
	 If[midDepth > 0,
		StringTemplate["  x.l(expmt, e, `1`) = xopt.l(expmt, e, `1`) ;"][massIndex["k", midDepth]],
		""],
	 (* set objective coefficients *)
	 "  rn.fx(j) = rangenum(roname, j) ;",
	 "  rd.fx(j) = rangeden(roname, j) ;",
	 "  rc.fx = rangeconst(roname) ;",
	 "  Solve fluxci using NLP minimizing fr ;",
	 "  put rangelow ;",
	 "  put roname.tl, fr.l / ;",
	 "  put fluxlow ;",
	 "  put roname.tl ;",
	 "  loop(j, put v.l(j):0:10 ; ) ;", 
	 "  put / ;",
	 ") ;",
	 "putclose rangelow ;",
	 "putclose fluxlow ;"
	 }]


(* ::Subsection::Closed:: *)
(*Ratio confidence intervals block*)


(* ::Text:: *)
(*This uses a loop to minimize/maximize flux ratios while keeping the optimum bound <  fopt + chisq, a chi square quantil at the desired confidence level (see Antoniewicz et al Metab Eng 2006, eq 19)*)
(*We initialize from the optimum point vopt, xopt*)


MakeFluxRatioCIBlock[modelName_String, chisq_] :=
	{"** MakeFluxRatioCIBlock[] **",
	 "fbound.fx = f.l + " <> ToString[chisq] <> " ;",
	(* upper range problems *)
	 "file ratiohi /" <> modelName <> "_ratiohi.tsv/ ;",
	 "ratiohi.pc = 6;",
	 "put ratiohi;",
	 "loop((ir,jr)$(S(ir,jr) > 0),",
	 "  display 'Upper CI ir = ', ir, ', jr = ', jr ;",
	 "  v.l(j) = vopt.l(j) ;",
	 "  x.l(expmt, e," <> massIndex["k", midDepth] <>") = " <>
		"xopt.l(expmt, e, " <> massIndex["k", midDepth] <>") ;",
	 "  ci.fx(i,j) = 0;",
	 "  ci.fx(ir,jr) = 1;",
	 "  Solve fluxratioci using NLP maximizing fci ;",
	 "  put ir.tl, jr.tl, fci.l:0:10/ ;",
	 ") ;",
	 "putclose ratiohi;",
	(* lower range problems *)
	 "file ratiolow /" <> modelName <> "_ratiolow.tsv/ ;",
	 "ratiolow.pc = 6;",
	 "put ratiolow;",
	 "loop((ir,jr)$(S(ir,jr) > 0),",
	 "  display 'Lower CI ir = ', ir, ', jr = ', jr ;",
	 "  v.l(j) = vopt.l(j) ;",
	 "  x.l(expmt, e," <> massIndex["k", midDepth] <>") = " <>
		"xopt.l(expmt, e," <> massIndex["k", midDepth] <>") ;",
	 "  ci.fx(i,j) = 0;",
	 "  ci.fx(ir,jr) = 1;",
	 "  Solve fluxratioci using NLP minimizing fci ;",
	 "  put ir.tl, jr.tl, fci.l:0:10/ ;",
	 ") ;",
	 "putclose ratiolow;"}


(* ::Subsection:: *)
(*Complete export function*)


(* ::Text:: *)
(*This joins the various blocks created by the other functions and saves as a text file.*)
(*If the FluxCI option is given, it must be a list of Reaction for which to get net flux confidence intervals.*)


Options[GAMSExport] = {
	"FluxCI" -> {}, "FluxRatioCI" -> {}, "ConfidenceLevel" -> 0.9, "TestCase" -> False};


GAMSExport[gamsDir_String, modelName_String, M_MetabolicNetwork,
			emuReact:{_EMUReaction...}, emuEq:{_EMUEquation..}, 
			substrateMIDs:{{_Rule..}..},
			measFluxNames:{_String..}, measReaction:{_String..}, measFluxMean_?VectorQ, measFluxSD_?VectorQ, fluxMixing_SparseArray,
			expName:{_String..}, measName:{_String..}, measMID_List, midsd_List,
			measEMU:{_EMU..}, emuMixing_SparseArray, 
			initFlux_FluxState, OptionsPattern[]] :=
	Block[{B, emap, rl, fluxSol, \[Epsilon]f = N[10^(-6)], midDepth, chisq},
		(* calculate balanced network *)
		B = BoundaryNetwork[M];
		rl = gamsReactionList[M, B];
		emap = EMUFluxMap[M];
		midDepth = Length[Dimensions[First[InternalEMUs[First[emuEq]]]]];
		(* simulate MIDs at initial fluxes, for each experiment *)
		fluxSol = Table[EMUSimulate[emuEq, s, emap[initFlux]], {s, substrateMIDs}];
		Export[FileNameJoin[{gamsDir, modelName <> ".gms"}],
			Join[
				MakeHeader[],
				MakeStoichiometryBlock[M, B],
				MakeEMUBlock[emuEq],
				MakeMeasurementsBlock[
					measFluxNames, measReaction, measFluxMean, measFluxSD, fluxMixing,
					expName, measName, measMID, midsd, measEMU, emuMixing],
				MakeEMUReactionsBlock[emuReact, rl],
				MakeEMUCondensationsBlock[emuReact, emuEq, rl],
				MakeVariablesBlock[emuEq, expName, substrateMIDs, 0.001, 100000],
				MakeEquationsBlock[midDepth],
				MakeInitialPointBlock[emuEq, emap, initFlux, expName, fluxSol, rl],
				MakeIsotopeModelsBlock[],
				MakeOptimizationBlock[modelName, midDepth],
				If[Length[OptionValue["FluxCI"]] > 0,
					(* determine confidence intervals on net fluxes *)
					MakeNetFluxesCIBlock[modelName, midDepth, OptionValue["FluxCI"],
						Quantile[ChiSquareDistribution[1], OptionValue["ConfidenceLevel"]]],
					{}]
				],
			"Lines"]]


(* ::Text:: *)
(*A version for confidence intervals on ratios*)


GAMSExportRatioCI[gamsDir_String, modelName_String, M_MetabolicNetwork,
			emuReact:{_EMUReaction...}, emuEq:{_EMUEquation..}, 
			substrateMIDs:{{_Rule..}..},
			measFluxNames:{_String..}, measReaction:{_String..}, measFluxMean_?VectorQ, measFluxSD_?VectorQ, fluxMixing_SparseArray,
			expName:{_String..}, measName:{_String..}, measMID_List, midsd_List,
			measEMU:{_EMU..}, emuMixing_SparseArray, 
			initFlux_FluxState,
			ratioNames_, ratioNumerator_, ratioDenominator_, ratioConstant_, chisq_] :=
	Block[{B, emap, rl, fluxSol, \[Epsilon]f = N[10^(-6)], midDepth},
		(* calculate balanced network *)
		B = BoundaryNetwork[M];
		rl = gamsReactionList[M, B];
		emap = EMUFluxMap[M];
		midDepth = Length[Dimensions[First[InternalEMUs[First[emuEq]]]]];
		(* simulate MIDs at initial fluxes, for each experiment *)
		fluxSol = Table[EMUSimulate[emuEq, s, emap[initFlux]], {s, substrateMIDs}];
		Export[FileNameJoin[{gamsDir, modelName <> ".gms"}],
			Join[
				MakeHeader[],
				MakeStoichiometryBlock[M, B],
				MakeEMUBlock[emuEq],
				MakeMeasurementsBlock[
					measFluxNames, measReaction, measFluxMean, measFluxSD, fluxMixing,
					expName, measName, measMID, midsd, measEMU, emuMixing],
				MakeEMUReactionsBlock[emuReact, rl],
				MakeEMUCondensationsBlock[emuReact, emuEq, rl],
				MakeVariablesBlock[emuEq, expName, substrateMIDs, 0.001, 100000],
				MakeEquationsBlock[midDepth],
				MakeInitialPointBlock[emuEq, emap, initFlux, expName, fluxSol, rl],
				MakeIsotopeModelsBlock[],
				MakeOptimizationBlock[modelName, midDepth],
				(* determine confidence interval on rational functions *)
				MakeFluxCIBlock[
					modelName, midDepth,
					ratioNames, ratioNumerator, ratioDenominator, ratioConstant,
					chisq]
				],
			"Lines"]]


(* ::Text:: *)
(*Version without 13C data (fit fluxes only)*)


GAMSExport[gamsDir_String, modelName_String, M_MetabolicNetwork,
			measFluxNames:{_String..}, measReaction:{_String..}, measFluxMean_?VectorQ, measFluxSD_?VectorQ, fluxMixing_SparseArray,
			initFlux_FluxState, OptionsPattern[]] :=
	Block[{B, emap, rl, \[Epsilon]f = N[10^(-6)], chisq},
		(* calculate balanced network *)
		B = BoundaryNetwork[M];
		rl = gamsReactionList[M, B];
		emap = EMUFluxMap[M];
		Export[FileNameJoin[{gamsDir, modelName <> ".gms"}],
			Join[
				MakeHeader[],
				MakeStoichiometryBlock[M, B],
				MakeMeasurementsBlock[
					measFluxNames, measReaction, measFluxMean, measFluxSD, fluxMixing],
				MakeVariablesBlock[0.001, 100000],
				MakeEquationsBlock[],
				MakeInitialPointBlock[emap, initFlux, rl],
				MakeFluxModelsBlock[],
				MakeOptimizationBlock[modelName],
				If[Length[OptionValue["FluxCI"]] > 0,
					MakeNetFluxesCIBlock[modelName, 0, OptionValue["FluxCI"],
						Quantile[ChiSquareDistribution[1], OptionValue["ConfidenceLevel"]]],
					{}]
				],
			"Lines"]]


(* ::Section::Closed:: *)
(*Find balanced flux ratios*)


(* ::Text:: *)
(*This is used to generate "interesting" flux distributions, such that the ratios for selected metabolites (typically those used for tracing) are as close as possible to the ideal 1 / n*)


(* ::Text:: *)
(*This uses a sink reaction (sum of all output fluxes = constant) and no flux measurements. Results from solving this problem will be scaled to total output = 100 (arbitrary).*)


GAMSExportRatioBalanceProblem[modelName_String, M_MetabolicNetwork,
			initFlux_FluxState, selectedMet:{_Metabolite..}, OptionsPattern[]] :=
	Block[{B, emap, rl, \[Epsilon]f = N[10^(-6)], sinkFluxNames},
		(* calculate balanced network *)
		B = BoundaryNetwork[M];
		rl = gamsReactionList[M, B];
		emap = EMUFluxMap[M];
		(* "sink" (output) fluxes *)
		sinkFluxNames = BoundaryFluxNames[OutputMetabolites[M]];
		Export[modelName <> ".gms",
			Join[
				MakeHeader[],
				MakeStoichiometryBlock[M, B],
				{"Set ratiomet(i) \"metabolites to optimize ratios for\" / " 
					<> ListToString[MetaboliteShortName /@ selectedMet,", "] <> " / ;"},
				MakeVariablesBlock[0.001, 10000],
				{"Equations"},
				balancedRatioEquationsBlock[],
				sinkEquationsBlock[sinkFluxNames],
				balancedRatioObjectiveBlock[],
				MakeInitialPointBlock[emap, initFlux, rl],
				balancedRatioModelsBlock[],
				balancedRatioOptimizationBlock[modelName]],
			"Lines"]]


(* ::Text:: *)
(*The objective here attempts to drive all ratios close to 1 / n, where n is the number of ratios into the the same metabolite*)


balancedRatioEquationsBlock[] :=
	{"objective \"defines the L2 error objective function f\"",
	 "massbaleq(i) \"mass balance equation\"",
	 "sinkEquation \"list of reactions that must sum to 1\" ; ",
	 "massbaleq(massbal(i)) .. ",
	 "  sum(j, S(i,j)*v(j)) =e= 0 ;"}


balancedRatioObjectiveBlock[] :=
	{"objective .. ",
	 "f =e= sum((ratiomet(i),j)$(S(i,j) > 0), power(" <>
		"S(i,j)*v(j) / sum(jp$(S(i,jp) > 0), S(i,jp)*v(jp)) " <>
		" - 1 / sum(jp$(S(i,jp) > 0), S(i,jp)), 2)) ; "}


balancedRatioModelsBlock[] :=
	{"** balancedRatioModelsBlock[] ** ",
	 "Model",
	 "  balratiosopt / objective, massbaleq, sinkEquation / ;"}


balancedRatioOptimizationBlock[modelName_String] :=
	{"** balancedRatioModelsBlock[] **",
	 "Solve balratiosopt using NLP minimizing f ;",
	 "vopt.l(j) = v.l(j) ;",
	 "file fluxopt /" <> modelName <> "_fluxopt.tsv/ ;",
	 "fluxopt.pc = 6 ;",
	 "put fluxopt ; loop(j, put j.tl, v.l(j):0:10/); putclose fluxopt;" }


(* ::Section::Closed:: *)
(*Import GAMS result files*)


(* ::Subsection::Closed:: *)
(*ImportGAMSFluxes*)


(* ::Text:: *)
(*This yields a FluxState object*)


fluxOptFileName[gamsDir_String, modelName_String] := 
	FileNameJoin[{gamsDir, modelName <> "_fluxopt.tsv"}]


(* ::Text:: *)
(*Import a list of rules assigning flux values to reaction names*)


importFluxRules[gamsDir_String, modelName_String] :=
	Dispatch[Join[
		Map[(Rule@@#)&, Import[fluxOptFileName[gamsDir, modelName]]],
		(* assign zero to reverse (_r) fluxes for irreversible reactions *)
		{_String -> 0.}]]


ImportGAMSFluxes[gamsDir_String, modelName_String, M_MetabolicNetwork] :=
	ImportGAMSFluxes[gamsDir, modelName, M, EMUFluxMap[M]]


ImportGAMSFluxes[gamsDir_String, modelName_String, M_MetabolicNetwork, emap_EMUFluxMap] :=
	ImportGAMSFluxes[importFluxRules[gamsDir, modelName], M, emap]


ImportGAMSFluxes[repl_, M_MetabolicNetwork, emap_EMUFluxMap] :=
	Block[{B, nu, nr},
		B = BoundaryNetwork[M];
		(* list of uptake/fwd/rev/release fluxes *)
		nu = Length[InputMetabolites[M]];
		nr = Length[OutputMetabolites[M]];
		FluxState @@ ({Table[StringJoin[ReactionID[r], "_f"], {r, Reactions[M]}],
			 Table[StringJoin[ReactionID[r], "_r"], {r, Reactions[M]}],
			 Take[BoundaryFluxNames[B], nu],
			 Take[BoundaryFluxNames[B], -nr]} /. repl)]


(* ::Subsection::Closed:: *)
(*ImportGAMSIsotopes*)


(* ::Text:: *)
(*This yields a list of MIDs, one for each EMU listed*)


ImportGAMSIsotopes[gamsDir_String, modelName_String, expNames:{_String..}, emus:{_EMU..}] :=
	Block[{t,n,es,emuGams, xflat},
		t = Import[FileNameJoin[{gamsDir, modelName <> "_isotopes.tsv"}]];
		(* arrange isotope fractions into MID vectors *)
		t = GatherBy[t, Take[#, 2]&];
		emuGams = t[[All,1,{1,2}]]; (* list of {exp, emu} pairs *)
		t = Drop[t,None,None,2];
		(* order MIDs by isotope *)
		t = Table[Last /@ SortBy[x,Most], {x,t}];
		(* map to given list of experiments and EMUs *)
		Table[
			xflat = Thread[{exp, emuString /@ emus}] /. Dispatch[Thread[Rule[emuGams, t]]];
			Table[
				ArrayReshape[xflat[[i]], Dimensions[emus[[i]]]+1],
				{i, Length[emus]}],
			{exp, expNames}]]


(* ::Subsection::Closed:: *)
(*ImportGAMSEmuMixing*)


(* ::Text:: *)
(*This is the estimated mixing matrix for measurements--targetEMUs.*)


ImportGAMSEmuMixing[gamsDir_String, modelName_String, measName:{_String...}, measEMU:{_EMU..}] :=
	Block[{t,n,es,emuGams, xflat},
		t = Import[FileNameJoin[{gamsDir, modelName <> "_emumixing.tsv"}]];
		SparseArray[
			Thread[Rule[
				Thread[{ReindexList[t[[All,1]], measName],
						ReindexList[t[[All,2]], emuString /@ measEMU]}],
				t[[All,3]] ]],
			{Length[measName], Length[measEMU]}]]


(* ::Subsection::Closed:: *)
(*ImportGAMSFluxCI*)


(* ::Text:: *)
(*Import a flux vector at each CI boundary*)


fluxLowFileName[gamsDir_String, modelName_String] := 
	FileNameJoin[{gamsDir, modelName<>"_fluxlow.tsv"}]


fluxHiFileName[gamsDir_String, modelName_String] := 
	FileNameJoin[{gamsDir, modelName<>"_fluxhi.tsv"}]


importFluxNames[gamsDir_String, modelName_String] :=
	First /@ Import[fluxOptFileName[gamsDir, modelName]]


ImportGAMSFluxCI[gamsDir_String, modelName_String, M_MetabolicNetwork] :=
	ImportGAMSFluxCI[gamsDir, modelName, M, EMUFluxMap[M]]


ImportGAMSFluxCI[gamsDir_String, modelName_String, M_MetabolicNetwork, emap_EMUFluxMap] :=
	Block[{names, fluxes},
		names = importFluxNames[gamsDir, modelName];
		Transpose[Table[
			fluxes = Rest /@ Import[fn[gamsDir, modelName]];
			Table[
				ImportGAMSFluxes[
					Join[Thread[Rule[names, f]], {_String -> 0.}],
					M, emap],
				{f, fluxes}],
			{fn, {fluxLowFileName, fluxHiFileName}}]]]


(* ::Subsection::Closed:: *)
(*ImportGAMSFluxRatioCI*)


(* ::Text:: *)
(*This yields a sparse matrix of the same dimensions as Stoichiometry[M], with positive entries for flux ratios derived from forward reactions, and negative for reverse reactions.*)


ImportGAMSFluxRatioCI[dirName_String, modelName_String, M_MetabolicNetwork] :=
	{ImportGAMSFluxRatioCI[FileNameJoin[{dirName, modelName<>"_ratiolow.tsv"}], M],
	 ImportGAMSFluxRatioCI[FileNameJoin[{dirName, modelName<>"_ratiohi.tsv"}], M]}


ImportGAMSFluxRatioCI[fileName_String, M_MetabolicNetwork] :=
	Block[{met, rf, rr},
		met = MetaboliteShortName /@ Metabolites[M];
		rf = Table[StringJoin[ReactionID[r], "_f"], {r, Reactions[M]}];
		rr = Table[StringJoin[ReactionID[r], "_r"], {r, Reactions[M]}];
		convertGamsFluxRatios[met, rf, rr,
			Transpose[Import[fileName, "TSV"]]]]


convertGamsFluxRatios[met:{_String..}, rf:{_String..}, rr:{_String..}, {ml_, rl_, ratios_}] :=
	SparseArray[
		Join[
			Cases[
				Transpose[{ReindexList[ml, met], ReindexList[rl, rf], ratios}],
				{i_Integer, j_Integer, x_} -> Rule[{i,j}, x]],
			Cases[
				Transpose[{ReindexList[ml, met], ReindexList[rl, rr], ratios}],
				{i_Integer, j_Integer, x_} -> Rule[{i,j}, -x]]]]


(* ::Section::Closed:: *)
(*End package*)


End[];

EndPackage[];
