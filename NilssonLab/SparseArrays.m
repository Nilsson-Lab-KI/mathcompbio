
(* Some utilities for working with SparseArray objects. *)

BeginPackage["NilssonLab`SparseArrays`", {"NilssonLab`Utilities`"}];

SparseTotal::usage = "SparseTotal[A] behaves like the built-in Total[], but works on sparse arrays";
TransposeToFirst::usage = "TransposeToFirst[A, k] transposes a SparseArray A so that dimension k becomes the first dimension, and all other dimensions remain in the same sequence";
TransposeToLast::usage = "TransposeToLast[A, k] transposes a SparseArray A so that dimension k becomes the last dimension, and all other dimensions remain in the same sequence";
DotAt::usage = "DotAt[A, k, B, m] returns the dot product of dimension k of SparseArray A with dimension m of SparseArray B. Dot[A, k, v] returns the product of dimensions k of A with a vector v.";
TimesAt::usage = "TimesAt[A, k, B, m] returns the elementwise product of dimension(s) k of A with dimensions(s) m of B.";


Begin["NilssonLab`SparseArrays`Private`"];


(* Transpose an array so that dimension k in the becomes dimension 1,
   and all other dimensions remain in the same sequence. Note that this is not specific to sparse arrays *)

TransposeToFirst[A_?ArrayQ, 1] := A

TransposeToFirst[A_?ArrayQ, k_Integer] :=
	Transpose[A, Join[Range[k-1]+1, {1}, Range[k+2, ArrayDepth[A]]]]

(* Generalize to a list of dimensions; dimensions kl are put first, then follow all other dimensions *)

TransposeToFirst[A_?ArrayQ, kl:{_Integer..}] :=
	Transpose[A, dimVector[kl, ArrayDepth[A]]]

dimVector[kl:{_Integer..}, n_Integer] :=
	Block[{j = 1, l = Length[kl]+1},
		Table[If[MemberQ[kl, i], j++, l++], {i, n}]]


(* Transpose a SparseArray so that dimension k in the becomes the last dimension,
   and all other dimensions remain in the same sequence *)

TransposeToLast[A_?ArrayQ, k_Integer] :=
	Transpose[A, Join[Range[k-1], {ArrayDepth[A]}, Range[k, ArrayDepth[A]-1]]]


(* Total sum of sparse matrix elements *)

SparseTotal[A_SparseArray] := SparseArray[Table[1, {Length[A]}]].A

(*Total down to level n*)

SparseTotal[A_SparseArray, n_Integer /; n > 1] := 
	SparseTotal[SparseTotal[A], n-1]

SparseTotal[A_SparseArray, 1] := SparseTotal[A]

(*Total at level n*)

SparseTotal[A_SparseArray, {n_Integer}] := 
	SparseTotal[TransposeToFirst[A, n]]


(* Dot product at a particular dimension *)

DotAt[A_SparseArray, k_Integer, B_SparseArray, m_Integer] :=
	Dot[TransposeToLast[A,k], TransposeToFirst[B,m]]

DotAt[A_SparseArray, k_Integer, v_?VectorQ] :=
	Dot[TransposeToLast[A,k], v]


(* Elementwise product. Here B must match dimensions of A completely. *)

TimesAt[A_SparseArray, {}, B_SparseArray] := A

TimesAt[A_SparseArray, kl:{_Integer..}, B_SparseArray] :=
	Transpose[
		TransposeToFirst[A, kl]*B,
		Join[kl, Complement[Range[ArrayDepth[A]], kl]]]

TimesAt[A_SparseArray, k_Integer, B_SparseArray] :=
	TimesAt[A, {k}, B]


End[];

EndPackage[];
