
(* Some functions for working with isotope and isotopomer distributions *)

BeginPackage["NilssonLab`IsotopeDistributions`", {"NilssonLab`Utilities`"}];


FragmentIsotopomerDistribution::usage = "FragmentIsotopomerDistribution[d, f] generates the distribution of fragment d of a molecules with isotopomer distribution d."
EmpiricalIsotopomerDistribution::usage = "EmpiricalIsotopomerDistribution[{{isotopomer, frequency}, ...}] specifies an isotopomer distribution as a list ofspecific isotopomers and their frequencies";
EmpiricalMID::usage = "This encapsulates a specific mass isotopomer distribution, that is, an array of MI fractions.";
MIDNormalizedQ;
MIDNormalize;
NaturalIsotopomerDistribution;
NaturalMID;
MixtureIsotopomerDistribution;
TracerPurityIsotopomerDistribution;
TracerPurityMID;
MixtureMID;
MarginalizeToMID::usage = "MarginalizeToMID[dist] where dist is an isotopomer distritbution determines the corresponding specific MID.";
MIDEnrichment;
ProductIsotopomerDistribution::usage = "";
ImportMID;
ImportMixtureID;
ImportMIDTable::usage = "ImportMIDTable[fileName] imports MIDs from a tab-separated peak area file."
NACorrectionMatrix;
NACorrect;


Begin["NilssonLab`IsotopeDistributions`Private`"];

(* Isotopomer distributions (IDs) are given as lists {{isotopomer, fraction} ... }
   where isotopomer completely specifies the (positional) isotopomer state across
   the elements involved.  For example, {{1,0,1} is a isotopomer of one element
   and three atoms, of which atoms 1 and 3 are heavy. *)

(* Mass isotopomer distributions (MIDs) are specified as an array of MI fractions,
   whose dimensions depend on the number of elements *)


(* Empirical isotopomer distributions are given explicitly as a list of probabilities (fractions)
   for each isotopomer or mass isotopomers. *)

EmpiricalIsotopomerDistribution /:
MarginalizeToMID[EmpiricalIsotopomerDistribution[ed_List]] :=
	Block[{d, mp},
		(* MID dimensions *)
		d = Length /@ ed[[1,1]];
		mp = Transpose[{Map[Total, ed[[All,1]], {2}], ed[[All,2]]}];
		mp = GatherBy[mp, First];
		EmpiricalMID[
			Replace[
				Outer[List, Sequence@@Table[Range[n+1],{n,d}]],
				Append[Thread[Rule[mp[[All,1,1]] + 1, Total /@ mp[[All, All, 2]]]], _ -> 0],
				{Length[d]}]]]

(* This would be more efficient if the EmpiricalIsotopomerDistribution was stored as a replacement list*)

EmpiricalIsotopomerDistribution /:
PDF[EmpiricalIsotopomerDistribution[ed_List], x:{{(0|1)...}..}] :=
	Replace[x, Append[Map[Rule@@#&, ed], _ -> 0]]

(* Given an empirical isotopomer distribution ed, find the isotopomer distribution of a fragment
   specified by atom numbers f for each element, e.g. {{2,3}, {1}} for the fragment consisting of
   atoms 2,3 of element 1 and atom 1 of element 2. *)

FragmentIsotopomerDistribution[EmpiricalIsotopomerDistribution[ed_List], f_List] :=
	EmpiricalIsotopomerDistribution[
		Table[
			{Table[Part[ed[[i,1,j]], f[[j]]], {j,Length[f]}], ed[[i,2]]},
			{i, Length[ed]}]]			


(* EmpiricalMID encapsulates a specific (explicit) mass isotopomer distribution, that is,
   an array of MI fractions. The array must have n dimensions for n elements
   (e.g. for carbon only, a vector; for carbon-nitrogen, a matrix).
   The data distribution need not be normalized to 1, so we can store peak areas; see MIDNormalize *)

Format[EmpiricalMID[data_List]] :=
	"<MID dim = " <> ListToString[Dimensions[data],"x"] <> ">"

Total[EmpiricalMID[data_List]] ^:= Total[Flatten[data]]

Dimensions[EmpiricalMID[data_List]] ^:= Dimensions[data]

(* Let Part access the encapsulated array directly *)

Part[EmpiricalMID[data_List], index__] ^:= Part[data, index]

(* Let Normal extract the encapsulated data array *)

Normal[EmpiricalMID[data_List]] ^:= data


(* PDF yields individual MI fractions. Here nl is a list of mass isotopomers per element,
   for example {0,2} for M+0 in the first element (say, carbon) and M+2 in the second element 
   (say, nitrogen) *)

EmpiricalMID /:
PDF[EmpiricalMID[data_List], nl:{_Integer..}] := Extract[data, nl + 1] / Total[Flatten[data]]


(* Test if an MID sums to 1. *)

MIDNormalizedQ[mid_EmpiricalMID] := Chop[Total[mid]] == 1


(* Normalize an MIDs by dividing all data points by their total. If the total is zero,
   we assign an unlabeled MID {1, 0, 0, ... } *)

MIDNormalize[EmpiricalMID[mid_List]] :=
	Block[{t},
		t = Total[Flatten[mid]];
		If[Positive[t],
			EmpiricalMID[mid / t],
			EmpiricalMID[
				Normal[SparseArray[Table[1, {Depth[mid]-1}] -> 1, Dimensions[mid]]]]]]


(* Isotopic enrichment, only defined for single-element MIDs *)

MIDEnrichment[e_EmpiricalMID] := MIDEnrichment[Normal[MIDNormalize[e]]]

MIDEnrichment[x_?VectorQ] := Range[0, Length[x]-1].x / (Length[x]-1)


(* A natural distribution for compounds with dimensions al and isotope probabilities pl
   for each element, sampled randomly and independently. Note that one or more dimensions
   can be empty (zero atoms). The length of each  xl[[i]] must equal the corresponding al[[i]] *)

NaturalIsotopomerDistribution /:
PDF[NaturalIsotopomerDistribution[al:{_Integer..}, pl_List], xl:{{(0|1)...}..}] :=
	Block[{n1, n0},
		n1 = Total /@ xl; n0 = al - n1;
		Times @@ Table[
			Product[pl[[i]]^n1[[i]] * (1 - pl[[i]])^n0[[i]], {i, Length[pl]}],
			{i, Length[xl]}]] /; (Length[al] == Length[xl])

NaturalIsotopomerDistribution /:
PDF[NaturalIsotopomerDistribution[_, _], {}] := 1


(* For a given fragment specified by a list of atoms numbers (only fragment size matters in this case) *)

MarginalizeToMID[
		NaturalIsotopomerDistribution[al:{_Integer..}, pl_List]] :=
	NaturalMID[al, pl] /; (Length[al] == Length[pl])

FragmentIsotopomerDistribution[NaturalIsotopomerDistribution[al_List, pl_List], f_List] :=
	NaturalIsotopomerDistribution[Length /@ f, pl]

(* Natural mass isotopomer distribution for a compound with dimensions given by
   list nl and list of isotope probabilities pl.  Since elements are independent, this is
   a product of the distributions of the individual elements. *)

NaturalMID /: PDF[NaturalMID[nl_List, pl_List], xl:{_Integer..}] := 
	Product[
		PDF[BinomialDistribution[nl[[i]], pl[[i]]], xl[[i]]],
		{i, Length[pl]}] /; (Length[pl] == Length[xl])

Normal[dist:NaturalMID[nl_List, pl_List]] ^:= 
	Map[PDF[dist, #]&, Outer[List, Sequence @@ Table[Range[0,n], {n, nl}]], {Length[nl]}]



(* ::Text:: *)
(*This models a compound for which isotopomer t (tracer) has probability pt and other isotopomers follow a natural isotopomer distribution with heavy isotope fractions p1.*)



(* ::Subsection:: *)
(*TracerPurityIsotopomerDistribution*)


(* ::Text:: *)
(*Random model, with independent atoms. For a tracer isotopomer t, with heavy probability ptl for the labeled atoms ("purity"), one  for each element; and (natural) heavy probability p1l for the unlabeled atoms, for each element.*)
(*This is a product distribution. In the special case where all atoms either labeled or unlabeled, we obtain a natural distribution.*)


TracerPurityIsotopomerDistribution /:
PDF[TracerPurityIsotopomerDistribution[t:{{(0|1)...}..}, ptl_List, p1l_List], x:{{(0|1)...}..}] :=
	Times @@ Table[
		Times @@ (Transpose[{t[[i]],x[[i]]}] /. {
			{0,0} -> 1-p1l[[i]], {0,1} -> p1l[[i]], {1,0}-> 1-ptl[[i]], {1,1}->ptl[[i]]}),
		{i, Length[x]}]  /; (Length[ptl] == Length[x])


(* ::Subsubsection::Closed:: *)
(*Marginalize to form MID*)


MarginalizeToMID[
		TracerPurityIsotopomerDistribution[t:{{(0|1)...}..}, ptl_List, p1l_List]] :=
	TracerPurityMID[t, ptl, p1l]


(* ::Subsubsection::Closed:: *)
(*Fragment distribution*)


FragmentIsotopomerDistribution[
		TracerPurityIsotopomerDistribution[t:{{(0|1)...}..}, ptl_List, p1l_List],
		f_List] := 
	TracerPurityIsotopomerDistribution[
		Table[Part[t[[i]], f[[i]]],{i, Length[t]}],
		ptl, p1l]


(* ::Text:: *)
(* Model tracer MIDs as a convolution of two binomials, one for the tracer purity process
   and one for the natural process for non-labeled atoms, since we do not know if the heavy
   isotopomers originate from the tracer or the natural process. *)

TracerPurityMID /:
PDF[TracerPurityMID[t:{{(0|1)...}..}, ptl_List, p1l_List], x:{_Integer..}] :=
	Product[
		Sum[
			PDF[NaturalMID[{Total[t[[i]]]}, {ptl[[i]]}], {y}] * 
				PDF[NaturalMID[{Total[1-t[[i]]]}, {p1l[[i]]}], {x[[i]] - y}],
			{y, 0, x[[i]]}],
		{i, Length[x]}]

Normal[dist:TracerPurityMID[t:{{(0|1)...}..}, ptl_List, p1l_List]] ^:= 
	Map[PDF[dist, #]&, Outer[List, Sequence @@ Table[Range[0,n], {n, Length /@ t}]], {Length[t]}]


(* Form a new isotopomer distribution as a linear mixture of a list of distributions
   with given mixture coefficients. We assume all distribution have the same domain
   (metabolite/fragment size). *)

MixtureIsotopomerDistribution /:
PDF[MixtureIsotopomerDistribution[dist_List, coeff_?VectorQ], x:{{(0|1)...}..}] :=
	Map[PDF[#, x]&, dist].coeff


MarginalizeToMID[MixtureIsotopomerDistribution[dist_List, coeff_?VectorQ]] :=
	MixtureMID[MarginalizeToMID /@ dist, coeff]

FragmentIsotopomerDistribution[MixtureIsotopomerDistribution[dist_List, coeff_?VectorQ], f_List] := 
	MixtureIsotopomerDistribution[
		Map[FragmentIsotopomerDistribution[#, f]&, dist], coeff]

(* Syntactic sugar for the case of 2 distributions, for compatibility with older version *)

MixtureIsotopomerDistribution[{id1_, a1_?NumericQ}, {id2_, a2_?NumericQ}] :=
	MixtureIsotopomerDistribution[{id1, id2}, {a1, a2}]

MixtureIsotopomerDistribution[{id_}, {1.0|1}] := id


(* Linear mixtures of MIDs *)

MixtureMID /:
PDF[MixtureMID[mids_List, coeff_?VectorQ], x:{_Integer..}] :=
	Map[PDF[#, x]&, mids].coeff

Normal[MixtureMID[mids_List, coeff_?VectorQ]] ^:= 
	Map[Normal, mids].coeff

MixtureMID[{mid1_, a1_?NumericQ}, {mid2_, a2_?NumericQ}] :=
	MixtureMID[{mid1, mid2}, {a1, a2}]

MixtureMID[{mid_}, {1.0|1}] := mid


(* Natural abundance correction *)

NACorrectionMatrix[n_,p_] :=
	Transpose[Table[
		Join[Table[0,{i}],PDF[BinomialDistribution[n-i,p],Range[0,n-i]]],
		{i,0,n}]]

NACorrect[x_?VectorQ, p_] :=
	LinearSolve[NACorrectionMatrix[Length[x]-1, p], x]


End[];

EndPackage[];
