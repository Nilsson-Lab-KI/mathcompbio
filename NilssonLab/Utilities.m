

BeginPackage["NilssonLab`Utilities`"];


ListToString::usage = "ListToString[list, sep] encodes list as a delimited string, with sep defaulting to |";
ReindexList::usage = "ReindexList[x_List, y_List] finds the position in y of each element of x. Elements of x not found in y are set to Missing, and Missing values in x are skipped. ReindexList[x, y, k] performs the replace at level k.";
PositivePart::usage = "PositivePart[x] gives the positive part of a number (or list of numbers) x.";
NegativePart::usage = "NegativePart[x] gives the negative part of a number (or list of numbers) x.";
NonZeroElements::usage = "NonZeroElements[a] returns a list of the non-zero elements of a SparseArray a";
NonZeroIndices::usage = "NonZeroIndices[a] returns a list of indices of the non-zero elements of a SparseArray a";
SparseArrayFlatten;
SparseIdentity;
SparseZero;
SparseReplace;
ListReshape;


Begin["NilssonLab`Utilities`Private`"];


(* Encode a list as a delimited string *)


ListToString[x_List] := ListToString[x, "|"]

ListToString[x_List/;Length[x]>0, sep_String] :=
	StringJoin @@ Table[ToString[z] <> sep, {z,Most[x]}] <> ToString[Last[x]]

ListToString[{}, sep_String] := ""


(* ReindexList *)

ReindexList[x_List, y_List] := ReindexList[x, y, 1]

ReindexList[x_List, y_List, level_Integer] :=
	Replace[x, Dispatch[Join[
		{Missing -> Missing},
		Thread[Rule[y, Range[Length[y]]]],
		{_ -> Missing}]], {level}]



(* Positive / negative part of arrays *)

PositivePart[X_SparseArray] :=
	Block[{r},
		r = Most[ArrayRules[X]];
		SparseArray[Pick[r, r[[All,2]], _?Positive], Dimensions[X]]]

PositivePart[x_?NumericQ] := UnitStep[x]*x

PositivePart[\[Infinity]] := \[Infinity]

PositivePart[-\[Infinity]] := 0

PositivePart[x_List] := PositivePart /@ x


NegativePart[X_SparseArray] :=
	Block[{r},
		r = Most[ArrayRules[X]];
		-SparseArray[Pick[r, r[[All,2]], _?Negative], Dimensions[X]]]

NegativePart[x_?NumericQ] := -UnitStep[-x]*x

NegativePart[\[Infinity]] := 0

NegativePart[-\[Infinity]] := \[Infinity]

NegativePart[x_List] := NegativePart /@ x


(* Non-zero elements and indices of sparse arrays *)


NonZeroElements[x_SparseArray] :=
	Flatten[Most[ArrayRules[x]][[All,2]]]

NonZeroElements[x_?VectorQ] :=
	NonZeroElements[SparseArray[x]]


NonZeroIndices[x_SparseArray] :=
	Most[ArrayRules[x]][[All,1]]

NonZeroIndices[x_?VectorQ] :=
	NonZeroIndices[SparseArray[x]]


(* ArrayFlatten crashes with empty sparse matrices, use this instead *)

SparseArrayFlatten[A_] :=
	Join @@ Table[
		Transpose[Join @@ (Transpose/@a)],
		{a, A}]


(*Sparse array-compatible replace function*)

SparseReplace[A_SparseArray, r_] :=
	SparseArray[ArrayRules[A] /. r, Dimensions[A]]

SparseReplace[{}, r_] := {}


(* Sparse identity matrix *)

SparseIdentity[n_Integer] := SparseArray[{i_, i_} -> 1, {n,n}]


(* Sparse zero (empty) matrix *)

SparseZero[m_Integer,n_Integer] := SparseArray[{}, {m, n}]


(* A recursive version of ArrayReshape that works on arbitrary lists, assuming that
   Flatten[a] is the same length as the vector x *)

ListReshape[x_List, a_List/;Depth[a]>2] :=
	Block[{n, xs},
		n = Replace[a, {l_List :> Length[Flatten[l]],_->1}, 1];
		xs = SplitVector[x, n];
		Table[ListReshape[xs[[i]],a[[i]]],{i,Length[a]}]]


End[];

EndPackage[];
