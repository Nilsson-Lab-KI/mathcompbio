

BeginPackage["NilssonLab`Utilities`"];


ListToString::usage = "ListToString[list, sep] encodes list as a delimited string, with sep defaulting to |";
ReindexList::usage = "ReindexList[x_List, y_List] finds the position in y of each element of x. Elements of x not found in y are set to Missing, and Missing values in x are skipped. ReindexList[x, y, k] performs the replace at level k.";


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


End[];

EndPackage[];
