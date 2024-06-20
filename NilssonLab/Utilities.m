

BeginPackage["NilssonLab`Utilities`"];


ListToString::usage = "ListToString[list, sep] encodes list as a delimited string, with sep defaulting to |";


Begin["NilssonLab`Utilities`Private`"];


(* Encode a list as a delimited string *)


ListToString[x_List] := ListToString[x, "|"]

ListToString[x_List/;Length[x]>0, sep_String] :=
	StringJoin @@ Table[ToString[z] <> sep, {z,Most[x]}] <> ToString[Last[x]]

ListToString[{}, sep_String] := ""


End[];

EndPackage[];
