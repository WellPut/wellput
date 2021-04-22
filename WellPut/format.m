(* ::Package:: *)


BeginPackage["WellPut`"]


Begin["`Private`"]


(* charges *)

Format[Q[a_],TraditionalForm] := Subscript["Q", a];

(* wilson coefficients *)

(* performSum::usage = "performSum[...] performs the sum over the Wilson coefficients"; *)

(* masses *)

Format[x[a_,b_],TraditionalForm] := Subsuperscript["x", b, a];

Format[m[a_],TraditionalForm] := Subscript["m", a];

(* replaceX::usage = "replaceX[...] replaces all instances of x_a^b by the mass ratios"; *)

(* sums *)

Format[sum[a_,{f0,fields__}],TraditionalForm] :=
  Row[{Underscript["\[Sum]", Row[{fields}]], "(", If[Length[{a}]<2, HoldForm[a], HoldForm[Times[a]]], ")"}];
Format[sum[a_,{f1_,fields___}],TraditionalForm] :=
  Row[{Underscript["\[Sum]", Row[{f1,fields}]], "(", If[Length[{a}]<2, HoldForm[a], HoldForm[Times[a]]], ")"}];

(* couplings *)


Format[g["vff",c_,v_,fa_,fb_],TraditionalForm] := Subsuperscript[g, Row[{v,bar[fa],fb}], c];
Format[g["sff",c_,s_,fa_,fb_],TraditionalForm] := Subsuperscript[y, Row[{s,bar[fa],fb}], c];
Format[g["vvv",a_,b_,c_],TraditionalForm] := Subscript[g, Row[{a,b,c}]];
Format[g["vss",a_,b_,c_],TraditionalForm] := Subscript[g, Row[{a,b,c}]];
Format[g["vvs",a_,b_,c_],TraditionalForm] := Subscript[g, Row[{a,b,c}]];


Format[bar[i_],TraditionalForm] := OverBar[i];

(* Fields *)



(* Loop Functions *)

Format[F["VAZ" ,a__,_  ],TraditionalForm] := Subsuperscript["F", "V"  , "\[Gamma]Z"        ][a];
Format[F["VBZ" ,c_, a__],TraditionalForm] := Subsuperscript["F", "V"  , ToString[c]<>",BZ" ][a];
Format[F["VB'Z",c_, a__],TraditionalForm] := Subsuperscript["F", "V"  , ToString[c]<>",B'Z"][a];
Format[F["VZ"  ,a__    ],TraditionalForm] := Subsuperscript["F", "V"  , "Z"                ][a];
Format[F["V'Z" ,a__    ],TraditionalForm] := Subsuperscript["F", "V'" , "Z"                ][a];
Format[F["V''Z",a__    ],TraditionalForm] := Subsuperscript["F", "V''", "Z"                ][a];

Format[F["VSB" ,a__    ],TraditionalForm] := Subsuperscript["F", "VS" , "B"                ][a];
Format[F["VSZ" ,a__    ],TraditionalForm] := Subsuperscript["F", "VS" , "Z"                ][a];
Format[F["VS'Z",a__    ],TraditionalForm] := Subsuperscript["F", "VS'", "Z"                ][a];

Format[F["SA"  ,a_,_   ],TraditionalForm] := Subsuperscript["F", "S"  , "\[Gamma]"         ][a];
Format[F["SB"  ,a__    ],TraditionalForm] := Subsuperscript["F", "S"  , "B"                ][a];
Format[F["SZ"  ,a__    ],TraditionalForm] := Subsuperscript["F", "S"  , "Z"                ][a];
Format[F["S'Z" ,a__    ],TraditionalForm] := Subsuperscript["F", "S'" , "Z"                ][a];
Format[F["S''Z",a__    ],TraditionalForm] := Subsuperscript["F", "S''", "Z"                ][a];

Format[F["dS",  a_,_   ], TraditionalForm] := Subsuperscript["F", "S",   "d"                ][a];
Format[F["dS'", a_,_   ], TraditionalForm] := Subsuperscript["F", "S'",  "d"                ][a];
Format[F["dV",  a__,_  ], TraditionalForm] := Subsuperscript["F", "V",   "d"                ][a];
Format[F["dV'", a_,_   ], TraditionalForm] := Subsuperscript["F", "V'",  "d"                ][a];


End[]

EndPackage[]
