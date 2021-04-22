(* ::Package:: *)

BeginPackage["WellPut`"]

Begin["`Private`"]

wellPutInfo["Startup"] :=
  Print[
    Style["This Package is based on the work ",15],
		Style[Hyperlink["[2104.?????]","https://arxiv.org/abs/2104.?????"],15],"\n",
    Style["Type \"wellPutInfo[]\" for a description of all available functions.",15]];

wellPutInfo[] :=
  Print[
    getc::usage, "\n\n",
    F::usage, "\n\n",
    replaceFunctions::usage, "\n\n",
    g::usage, "\n\n",
    L::usage, "\n",
    R::usage, "\n",
    bar::usage, "\n\n",
    el::usage, "\n",
    Q::usage, "\n",
    m::usage, "\n",
    x::usage, "\n",
    mu::usage, "\n",
    sum::usage];

x[a_,a_] := 1;

x[a_,b_] /; (m[a] /. Lmassrule /. Fmassrule) == 0 := 0;

replaceX[arg_] := arg /. x[a_,b_] :> repx[a,b];

repx[a_,b_] := m[a]^2/m[b]^2 /. Lmassrule /. Fmassrule /. Vmassrule;

bar[L] := R;
bar[R] := L;

bar[bar[a_]] := a;
x[a_, bar[b_]] := x[a, b];
x[bar[a_], b_] := x[a, b];

Clear[getc];
Options[getc] := {
  PerformSum -> True,
  ReplaceCharges -> True,
  ReplaceCouplings -> True,
  ReplaceFunctions -> False,
  ReplaceZeroCouplings -> True,
  CouplingRules -> {},
  Leptons -> {},
  Quarks -> {},
  Scalars -> {},
  ZBosons -> {},
  Vectors -> {},
  Externals -> {}
  };

getc[7, {i_, j_}, opts : OptionsPattern[]] := getc[{L}, {i, j}, opts];

getc[{cq_}, {i_, j_}, OptionsPattern[]] := Block[
        {res},
        setModel @@ ({#, OptionValue[ToExpression@#]} & /@ {"Leptons", "Quarks", "Scalars","ZBosons","Vectors", "Externals"});
        res = dipolesum[{i, j}, cq];
        If[OptionValue[PerformSum] == True,
           res = performSum[res]];
        If[OptionValue[ReplaceZeroCouplings] == True,
           res = res /. Join[
	     Cases[OptionValue[CouplingRules], _[_, _[0, _]]],
             Cases[OptionValue[CouplingRules], _[_, 0]]]];
        If[OptionValue[ReplaceCouplings] == True,
           res = res /. OptionValue[CouplingRules]];
        If[OptionValue[ReplaceFunctions] == True,
           res = replaceFunctions[res]];
        If[OptionValue[ReplaceCharges] == True,
           res = res /. chargerule];
        res];

getc[9,{i_,j_,l_}, opts:OptionsPattern[]] := (getc[{L, R}, {i, j, l}, opts] + getc[{L, L}, {i, j, l}, {opts}])/2;

getc[10,{i_,j_,l_}, opts:OptionsPattern[]] := (getc[{L, R}, {i, j, l}, opts] - getc[{L, L}, {i, j, l}, {opts}])/2;

getc[{cq_,cl_},{i_,j_,l_}, OptionsPattern[]] := Block[
  {res},
  setModel@@({#,OptionValue[ToExpression@#]}& /@ {"Leptons","Quarks","Scalars","ZBosons","Vectors","Externals"});
  res = mygetc[{i,j,l},cq,cl];
  If[OptionValue[PerformSum] == True,
     res = performSum[res]];
  If[OptionValue[ReplaceZeroCouplings] == True,
     res = res /. Join[
       Cases[OptionValue[CouplingRules],_[_,_[0,_]]],
       Cases[OptionValue[CouplingRules],_[_,0]]]];
  If[OptionValue[ReplaceCouplings] == True,
     res = res /. OptionValue[CouplingRules]];
  If[OptionValue[ReplaceFunctions] == True,
     res = replaceFunctions[res]];
  If[OptionValue[ReplaceCharges] == True,
     res = res /. chargerule];
  res];

dipolesum[{i_, j_}, cq_] :=
  sum[(1/m[s1]^2)*g["sff", bar[cq], bar[s1], i, f1] *
      (m[f1] * g["sff", bar[cq], s1, f1, j] *  
       F["dS",x[f1, s1], {f1, s1}] + 
       m[j] * g["sff", cq, s1, f1, j] * 
       F["dS'",x[f1, s1], {f1, s1}]), {f1, s1}] +
  sum[(1/m[v1]^2) *
      g["vff", cq, bar[v1], i, f1] *
      m[j] * g["vff", cq, v1, f1, j] * 
      F["dV",x[f0, v1], x[f1, v1], {f1, v1}],
      {f0, f1, v1}] + 
  sum[(1/m[v1]^2) *
      g["vff", cq, bar[v1], i, f1] *
      m[f1] * g["vff", bar[cq], v1, f1, j] *
      F["dV'",x[f1, v1], {f1, v1}],
      {f1, v1}];

mygetc[{i_,j_,l_},cq_,cl_] :=
  sum[(1/m[v1]^2) *
      g["vff",cq,bar[v1],i,f1] *
      g["vff",cq,v1,f1,j] *
      el^2 Q[l] F["VAZ",x[f0,v1],x[f1,v1],{f1,v1}],
      {f0,f1,v1}] +
  sum[(1/m[v1]^2) *
      g["vff",cq, bar[v2], i, f1] *
      g["vff",cq, v1, f1, j] *
      g["vff",cl,bar[v1],l,f3] *
      g["vff",cl,v2,f3,l] *
      F["VBZ",cl,x[f0,v1],x[f1,v1],x[v1,v2],x[f3,v1]],
      {f0,f1,f3,v1,v2}] +
  sum[(1/m[v1]^2) *
      g["vff",cq, bar[v2], i, f1] *
      g["vff",cq, v1, f1, j] *
      g["vff",cl,v2,l,f3] *
      g["vff",cl,bar[v1],f3,l] *
      F["VB'Z",cl,x[f0,v1],x[f1,v1],x[v1,v2],x[f3,v1]],
      {f0,f1,f3,v1,v2}] +
  sum[(1/m[z1]^2) *
      g["vff",cl, z1, l, l] *
      g["vff",cq, bar[v2], i, f1] *
      g["vff",cq, v1, f1, j] *
      g["vvv",z1, bar[v1], v2] *
      F["V''Z",x[f0,v1],x[f1,v1],x[v1,v2]],
      {f0,f1,v1,v2,z1}] +
  sum[(1/m[z1]^2) *
      g["vff",cl, z1, l, l] *
      g["vff",cq, bar[v1], i, f2] *
      g["vff",cq, v1, f1, j] *
      g["vff",cq, z1, f2, f1] *
      F["VZ",x[f1,v1],x[f2,v1]],
      {f1,f2,v1,z1}] +
  sum[(1/m[z1]^2) *
      g["vff",cl, z1, l, l] *
      g["vff",cq, bar[v1], i, f2] *
      g["vff",cq, v1, f1, j] *
      g["vff",bar[cq], z1, f2, f1] *
      F["V'Z",x[f1,v1],x[f2,v1]],
      {f1,f2,v1,z1}] +
  sum[(1/m[s1]^2) *
      g["sff",bar[cq], bar[s1], i, f1] *
      g["sff",cq, s1, f1, j] *
      el^2 * Q[l] *
      F["SA",x[f1,s1],{f1,s1}],
      {f1,s1}] +
  sum[(1/m[s1]^2) *
      g["sff",bar[cq], bar[s2], i, f1] *
      g["sff",cq, s1, f1, j] *
      (g["sff",bar[cl],bar[s1],l,f3] g["sff",cl,s2,f3,l] -
	  g["sff",bar[cl],s2,l,f3] g["sff",cl,bar[s1],f3,l] ) *
      F["SB",x[f1,s1],x[s1,s2],x[f3,s1]],
      {f1,f3,s1,s2}]+
   sum[(1/m[z1]^2) *
      g["vff",cl, z1, l, l] * 
      g["sff",bar[cq], bar[s1], i, f1] * g["sff",cq, s1, f1, j]*
      g["vff",cq, z1, j, j]*F["SZ",x[f1,s1],x[f1,s1]],{s1,f1,z1}]+
  sum[(1/m[z1]^2) *
      g["vff",cl, z1, l, l] * 
      g["sff",bar[cq], bar[s1], i, f1] *
      g["sff",cq, s2, f1, j] * 
      g["vss",z1,s1,bar[s2]] *
      F["SZ",x[f1,s1],x[f1,s2]],{s1,s2,f1,z1}]+
   sum[(1/m[z1]^2) *
      g["vff",cl, z1, l, l] * 
      g["sff",bar[cq], bar[s1], i, f2] * g["sff",cq, s1, f1, j]* 
      (g["vff",cl, z1, f2, f1]* F["S'Z",x[f1,s1],x[f2,s1]]+
      g["vff",bar[cl], z1, f2, f1]* F["S''Z",x[f1,s1],x[f2,s1]]),
      {s1,f1,f2,z1}] +
   sum[(1/m[s1]^2) *
      ( g["vff",cq, bar[v1], i, f1] * g["sff",cq, s1, f1, j]+
      g["sff",bar[cq], bar[s1], i, f1] * g["sff",cq, v1, f1, j])*
      ( g["sff",bar[cq], bar[s1], l, f3] * g["vff",cq, v1, f3, l]+
      g["vff", cq, bar[v1], l, f3] * g["sff",cq, s1, f3, l])*
      F["VSB", x[f1,s1],x[s1,v1],x[f3,s1]], {s1,v1,f1,f3}] +
   sum[(1/m[z1]^2) *
       g["vff",cl, z1, l, l] * (
	 g["vff",cq, bar[v1], i, f1] *
	 g["sff",cq, s1, f1, j] *
	 g["vvs",z1,v1,bar[s1]]*
	 F["VSZ", x[f1,s1], x[f1,v1]] + 
	 g["sff",bar[cq], bar[s1], i, f1] *
	 g["vff",cq, v1, f1, j] *
	 g["vvs",z1, bar[v1], s1] *
	 F["VS'Z", x[f1,s1], x[f1,v1]]),
       {s1,v1,f1,z1}];

setModel[] := Block[
  {},
  chargerule = Union[
    Echargerule,
    Lchargerule,
    Fchargerule,
    Vchargerule,
    Zchargerule,
    Schargerule];
  allfields = Union[
    Elist,
    Llist,
    Flist,
    Vlist,
    Slist,
    Zlist]; ];

setModel[{type_String,fieldInfo_}, a___] := Block[
  {},
  setFields[type,fieldInfo];
  setModel[a]];

setFields["Leptons"  ,fieldInfo_] := 
  {Llist, Lnmb, Lchargerule, Lmassrule} = setFields[fieldInfo];
setFields["Quarks"   ,fieldInfo_] := 
  {Flist, Fnmb, Fchargerule, Fmassrule} = setFields[fieldInfo];
setFields["Vectors"  ,fieldInfo_] := 
  {Vlist, Vnmb, Vchargerule, Vmassrule} = setFields[fieldInfo];
setFields["Scalars"  ,fieldInfo_] := 
  {Slist, Snmb, Schargerule, Smassrule} = setFields[fieldInfo];
setFields["ZBosons"  ,fieldInfo_] := 
  {Zlist, Znmb, Zchargerule, Zmassrule} = setFields[fieldInfo];
setFields["Externals",fieldInfo_] := 
  {Elist, Echargerule, Emassrule} = setFields[fieldInfo][[{1,3,4}]];

setFields[{{}}] := {{},0,{},{}};

setFields[fieldinfo_] := Block[
  {Plist, Pnmb, Pchargerule, Pmassrule},
  Plist = Transpose[fieldinfo][[1]];
  Pnmb = Length@Plist;
  Pchargerule = MapThread[Rule,{Q/@Plist,Transpose[fieldinfo][[3]]}];
  Pmassrule = MapThread[Rule,{m/@Plist,Transpose[fieldinfo][[2]]}];
  {Plist, Pnmb, Pchargerule, Pmassrule}];

(* perform field summations *)

chargeconservationrule = {
  g["vff"|"sff",_,a_,b_,c_] :> 0 /;
  ((Q[a]-Q[b]+Q[c] /. chargerule) =!= 0) &&
  (And @@ (MemberQ[allfields,#] & /@ {a,b,c})),
  g["vvv"|"vvs"|"vss",a_,b_,c_] :> 0 /;
  ((Q[a]+Q[b]+Q[c] /. chargerule) =!= 0) &&
  (And @@ (MemberQ[allfields,#] & /@ {a,b,c}))};


performSum[arg_] := arg /. sum -> dosum /. chargeconservationrule;

Qfn = (f1|f2);
Qln = (f3);
Qvn = (v1|v2);
Qzn = (z1);
Qsn = (s1|s2);

dosum[arg_,{},___] := arg;

dosum[arg_,{f0,fields__}] :=
  dosum[
    If[
      (m[Flist[[1]]] /. Fmassrule) == 0,
      arg /. x[f0,_] -> 0,
      arg /. x[f0,_] -> x[Flist[[1]],Vlist[[iv1]]]],
    {fields},"GIM"];

dosum[arg_,{fn:Qfn,fields___}       ] := dosum[Sum[arg /. fn :> Flist[[i]],{i,  Fnmb}],{fields}      ];
dosum[arg_,{fn:Qfn,fields___},"GIM" ] := dosum[Sum[arg /. fn :> Flist[[i]],{i,2,Fnmb}],{fields},"GIM"];

dosum[arg_,{ln:Qln,fields___},gim___] := dosum[Sum[arg /. ln :> Llist[[i]],{i,  Lnmb}],{fields}, gim];
dosum[arg_,{vn:Qvn,fields___},gim___] := dosum[Sum[arg /. vn :> Vlist[[i]],{i,  Vnmb}],{fields}, gim];
dosum[arg_,{sn:Qsn,fields___},gim___] := dosum[Sum[arg /. sn :> Slist[[i]],{i,  Snmb}],{fields}, gim];
dosum[arg_,{zn:Qzn,fields___},gim___] := dosum[Sum[arg /. zn :> Zlist[[i]],{i,  Znmb}],{fields}, gim];

(* general properties of the loop functions *)

F["VAZ",0,0,__] := 0;

F["VZ",x_,x_] := 0;
F["VZ",x_,y_] /; Not[OrderedQ[{x,y}]] := F["VZ",y,x];

F["VBZ",_,x_,x_,__] := 0;
F["V'Z",x_,x_] := 0;
F["V'Z",0 ,x_] := 0;
F["V'Z",x_,0 ] := 0;

F["V''Z",x_,y_,1] := 0;

(* analytic form of loop functions *)

replaceFunctions[arg_] :=
  arg /. { F["VAZ" ,a__] :> FVAZ[a],
	   F["VBZ" ,a__] :> FVBZ[a],
	   F["VB'Z",a__] :> FVBpZ[a],
	   F["VZ"  ,a__] :> FVZ[a],
	   F["V'Z" ,a__] :> FVpZ[a],
	   F["V''Z",a__] :> FVppZ[a],
	   F["VSB" ,a__] :> FVSB[a],
	   F["VSZ" ,a__] :> FVSZ[a],
	   F["VS'Z",a__] :> FVSpZ[a],
	   F["SA"  ,a__] :> FSA[a],
	   F["SB"  ,a__] :> FSB[a],
	   F["SZ"  ,a__] :> FSZ[a],
	   F["S'Z" ,a__] :> FSpZ[a],
	   F["S''Z",a__] :> FSppZ[a],
	   F["dS",  a__] :> FdS[a],
	   F["dS'", a__] :> FdSp[a],
	   F["dV",  a__] :> FdV[a],
	   F["dV'", a__] :> FdVp[a]}

fvaz[x_,{f_,v_}] :=
  ((1 + 7*x*(-3 + 2*x))/(12*(-1 + x)^3) + 
   ((-4 + (16 - 9*x)*x)*Log[x])/(6*(-1 + x)^4))*Q[f] + 
  ((-37 + x*(87 + 2*x*(-16 + 3*(-3 + x)*x)))/(12*(-1 + x)^3) + 
   (x*(6 + x*(-3 + 2*x)*(5 + 4*x))*Log[x])/(6*(-1 + x)^4))*Q[v];

FVAZ[x0_,x_,{f_,v_}] /; x0 =!= 0 && x =!= 0 :=
  fvaz[x,{f,v}] - fvaz[x0,{f,v}];

FVAZ[0,x_,{f_,v_}] /; x =!= 0 :=
  ((8 - 7*(-3 + x)*(-2 + x)*x)/(12*(-1 + x)^3) + 
   ((-4 + (16 - 9*x)*x)*Log[x])/(6*(-1 + x)^4) + 
   (2*Log[mu^2/m[v]^2])/3)*Q[f] + 
  ((x*(-24 + x*(79 + x*(-55 + 6*x))))/(12*(-1 + x)^3) + 
   (x*(6 + x*(-3 + 2*x)*(5 + 4*x))*Log[x])/(6*(-1 + x)^4))*Q[v];

FVAZ[0,0,{f_,v_}] := 0;

FVBZ[L,x0_,x_,y_,z_] /;
x0 =!= 0 && x =!= 0 && y =!= 1 && z =!=0 :=
  (2*x*y - 2*x0*y - 3/(-1 + x*y) + 3/(-1 + x0*y))/4 + 
  (x*y*(3*x^2*y*(-2 + x + x*y) - (-1 + x)*(4 + x*y*(-2 + x*y))*z)*Log[x])/
  (4*(-1 + x)*(-1 + x*y)^2*(x - z)) + 
  (x0*y*(-3*x0^2*y*(-2 + x0 + x0*y) + (-1 + x0)*(4 + x0*y*(-2 + x0*y))*z)*Log[x0])/
  (4*(-1 + x0)*(-1 + x0*y)^2*(x0 - z)) - 
  (3*(x - x0)*y^2*(1 - x0 + z - y*z + x*(-1 + x0*y*(2 - y + (-1 + y)*y*z)))*Log[y])/
  (4*(-1 + y)*(-1 + x*y)^2*(-1 + x0*y)^2*(-1 + y*z)) + 
  ((x - x0)*y*z^2*(4 + (-4 + y*(-4 + z))*z)*Log[z])/
  (4*(x - z)*(-1 + z)*(-x0 + z)*(-1 + y*z));

FVBZ[L,0,x_,y_,z_] /; x =!= 0 && y =!= 1 && z =!=0 :=
  Limit[FVBZ[L,x0,x,y,z],x0->0];

FVBZ[L,x0_,0,y_,z_] /; y =!= 1 && z =!=0 :=
  Limit[FVBZ[L,x0,x,y,z],x->0];

FVBZ[L,x0_,x_,1,z_] /; z =!=0 :=
  Limit[FVBZ[L,x0,x,y,z],y->1];

FVBZ[L,x0_,x_,y_,0] :=
  Limit[FVBZ[L,x0,x,y,z],z->0];
  
FVBZ[R,x0_,x_,y_,z_] /;
x0 =!= 0 && x =!= 0 && y =!= 1 && z =!=0 :=
  1/4 (2 x y-2 x0 y-3/(-1+x y)+3/(-1+x0 y))+(x y (3 x (4+x y (-6+x+x y))-(-1+x) (4+x y (-2+x y)) z) Log[x])/(4 (-1+x) (-1+x y)^2 (x-z))+
  (x0 y (-3 x0 (4+x0 y (-6+x0+x0 y))+(-1+x0) (4+x0 y (-2+x0 y)) z) Log[x0])/(4 (-1+x0) (-1+x0 y)^2 (x0-z))-
  (3 (x-x0) y^2 (-3+x0 (-1+4 y)+z-y z+x (-1+y (4+x0 (2-5 y+(-1+y) y z)))) Log[y])/(4 (-1+y) (-1+x y)^2 (-1+x0 y)^2 (-1+y z))+
  ((x-x0) y (-4+z) z^2 (-4+y z) Log[z])/(4 (x-z) (-1+z) (-x0+z) (-1+y z));
  
FVBZ[R,0,x_,y_,z_] /; x =!= 0 && y =!= 1 && z =!=0 :=
  Limit[FVBZ[R,x0,x,y,z],x0->0];

FVBZ[R,x0_,0,y_,z_] /; y =!= 1 && z =!=0 :=
  Limit[FVBZ[R,x0,x,y,z],x->0];

FVBZ[R,x0_,x_,1,z_] /; z =!=0 :=
  Limit[FVBZ[R,x0,x,y,z],y->1];

FVBZ[R,x0_,x_,y_,0] :=
  Limit[FVBZ[R,x0,x,y,z],z->0];
  
FVBpZ[L,x0_,x_,y_,z_] /;
x0 =!= 0 && x =!= 0 && y =!= 1 && z =!=0 :=
  1/4 (2 x y-2 x0 y-3/(-1+x y)+3/(-1+x0 y))+(x y (3 x (4+x y (-6+x+x y))-(-1+x) (4+x y (-2+x y)) z) Log[x])/(4 (-1+x) (-1+x y)^2 (x-z))+
  (x0 y (-3 x0 (4+x0 y (-6+x0+x0 y))+(-1+x0) (4+x0 y (-2+x0 y)) z) Log[x0])/(4 (-1+x0) (-1+x0 y)^2 (x0-z))-
  (3 (x-x0) y^2 (-3+x0 (-1+4 y)+z-y z+x (-1+y (4+x0 (2-5 y+(-1+y) y z)))) Log[y])/(4 (-1+y) (-1+x y)^2 (-1+x0 y)^2 (-1+y z))+
  ((x-x0) y (-4+z) z^2 (-4+y z) Log[z])/(4 (x-z) (-1+z) (-x0+z) (-1+y z));
  
FVBpZ[L,0,x_,y_,z_] /; x =!= 0 && y =!= 1 && z =!=0 :=
  Limit[FVBpZ[L,x0,x,y,z],x0->0];

FVBpZ[L,x0_,0,y_,z_] /; y =!= 1 && z =!=0 :=
  Limit[FVBpZ[L,x0,x,y,z],x->0];

FVBpZ[L,x0_,x_,1,z_] /; z =!=0 :=
  Limit[FVBpZ[L,x0,x,y,z],y->1];

FVBpZ[L,x0_,x_,y_,0] :=
  Limit[FVBpZ[L,x0,x,y,z],z->0];

FVBpZ[R,x0_,x_,y_,z_] /;
x0 =!= 0 && x =!= 0 && y =!= 1 && z =!=0 :=
  (2*x*y - 2*x0*y - 3/(-1 + x*y) + 3/(-1 + x0*y))/4 + 
  (x*y*(3*x^2*y*(-2 + x + x*y) - (-1 + x)*(4 + x*y*(-2 + x*y))*z)*Log[x])/
  (4*(-1 + x)*(-1 + x*y)^2*(x - z)) + 
  (x0*y*(-3*x0^2*y*(-2 + x0 + x0*y) + (-1 + x0)*(4 + x0*y*(-2 + x0*y))*z)*Log[x0])/
  (4*(-1 + x0)*(-1 + x0*y)^2*(x0 - z)) - 
  (3*(x - x0)*y^2*(1 - x0 + z - y*z + x*(-1 + x0*y*(2 - y + (-1 + y)*y*z)))*Log[y])/
  (4*(-1 + y)*(-1 + x*y)^2*(-1 + x0*y)^2*(-1 + y*z)) + 
  ((x - x0)*y*z^2*(4 + (-4 + y*(-4 + z))*z)*Log[z])/
  (4*(x - z)*(-1 + z)*(-x0 + z)*(-1 + y*z));

FVBpZ[R,0,x_,y_,z_] /; x =!= 0 && y =!= 1 && z =!=0 :=
  Limit[FVBpZ[R,x0,x,y,z],x0->0];

FVBpZ[R,x0_,0,y_,z_] /; y =!= 1 && z =!=0 :=
  Limit[FVBpZ[R,x0,x,y,z],x->0];

FVBpZ[R,x0_,x_,1,z_] /; z =!=0 :=
  Limit[FVBpZ[R,x0,x,y,z],y->1];

FVBpZ[R,x0_,x_,y_,0] :=
  Limit[FVBZ[L,x0,x,y,z],z->0];  
   
FVZ[x1_,x2_] /; x1 =!= 0 && x2 =!= 0 && x1 =!= x2 :=
  - (x1 + x2) / 2 +
  (x1*x2*Log[x1])/(x1 - x2) +
  (x1*x2*Log[x2])/(-x1 + x2);

FVZ[0,x2_] := Limit[FVZ[x1,x2],x1->0];

FVZ[x1_,0] := Limit[FVZ[x1,x2],x2->0];

FVpZ[x1_,x2_] /; x1 =!= 0 && x2 =!= 0 && x1 =!= x2 :=
  ((-4 + x2)*Sqrt[x1*x2])/(2*(-1 + x2)) -
  ((-4 + x1)*x1^(3/2)*Sqrt[x2]*Log[x1])/(2*(-1 + x1)*(x1 - x2)) + 
  (Sqrt[x1*x2]*(-3*x2^2 + x1*(4 + (-2 + x2)*x2))*Log[x2])/(2*(x1 - x2)*(-1 + x2)^2);

FVppZ[x0_,x_,y_]/; x0=!=0 && y=!=1 := (x (5+y+2 x (-1+y) y))/(-4+4 x y)-(x0 (5+y+2 x0 (-1+y) y))/(-4+4 x0 y)-(x (-1+y) (4+x (-10+3 x (-1+y)) y) Log[x])/(4 (-1+x) (-1+x y)^2)+(x0 (-1+y) (4+x0 (-10+3 x0 (-1+y)) y) Log[x0])/(4 (-1+x0) (-1+x0 y)^2)+(y ((x (5+y+x (3-y (5+4 y))))/(-1+x y)^2+(x0 (-5-y+x0 (-3+y (5+4 y))))/(-1+x0 y)^2) Log[y])/(4 (-1+y));
FVppZ[0,x_,y_]/; y=!=1 := Limit[FVppZ[x0,x,y],x0->0];

FVppZ[x0_,x_,1]/; x0=!=0 := 0;

FVppZ[0,x_,1]:= 0;


FVSB[x_,y_,z_]/; z=!=0 && y=!=1 :=(x^(3/2) (-4+x y) Sqrt[z] Log[x])/(4(-1+x) (-1+x y) (x-z))-(3 Sqrt[x] y Sqrt[z] Log[y])/(4(-1+y) (-1+x y) (-1+y z))-(Sqrt[x] z^(3/2) (-4+y z) Log[z])/(4 (x-z) (-1+z) (-1+y z));

FVSB[x_,1,z_]/; z=!=0 := Limit[FVSB[x,y,z],y->1];

FVSB[x_,y_,0]/; y=!=1 := 0;

FVSB[x_,1,0] := 0;

FVSZ[xs_,xv_]:=((5-4 xv) Sqrt[xv])/(4 (-1+xv))-(Sqrt[xv] (-4 xs+xv) Log[xs])/(4(-1+xs) (xs-xv))+(xv^(3/2) (-3+xs+2 xv) Log[xv])/(4 (-1+xv)^2 (-xs+xv));

FVSpZ[xs_,xv_]:=((1-2 xs) Sqrt[xv])/(4 (-1+xs))+(xs (-3+4 xs-xv) Sqrt[xv] Log[xs])/(4 (-1+xs)^2 (xs-xv))-(3 xs Sqrt[xv] Log[xv])/(4 (xs-xv) (-1+xv));

FSA[x_,{f_,s_}] :=
  Q[s]((2-7 x+11 x^2)/(36 (-1+x)^3)-(x^3 Log[x])/(6 (-1+x)^4))+
  Q[f]((16-29 x+7 x^2)/(36 (-1+x)^3)+((-2+3 x) Log[x])/(6 (-1+x)^4));
  
FSB[x_,y_,z_]/; z=!=0 && y=!=1 :=(x^2 y Log[x])/(4 (-1+x)(-1+x y)(x-z))+(y Log[y])/(4 (-1+y)(-1+x y)(-1+y z))+(y z^2 Log[z])/(4 (-1+z)(-x+z)(-1+y z));

FSB[x_,y_,0]/; y=!=1 := Limit[FSB[x,y,z], z->0];

FSB[x_,1,z_]/; z=!=0 := Limit[FSB[x,y,z], y->1];

FSB[x_,1,0]:= Limit[Limit[FSB[x,y,z], y->1],z->0];

FSZ[x1_,x2_]/; x1=!=x2 :=(1-2 x2)/(2 (x2-1))-(x2 Log[x1])/(2 (x1-1)(x1-x2))+((x1-1) x2 Log[x2])/(2 (x1-x2) (x2-1)^2);

FSZ[x1_,x1_] := Limit[FSZ[x1,x2],x2->x1];

FSpZ[x1_,x2_]/; x1=!=x2:=Sqrt[x1 x2]/(x1-x2) ((x1 Log[x1])/(x1-1) -(x2 Log[x2])/(x2-1));

FSpZ[x1_,x1_]:= Limit[FSpZ[x1,x2],x2->x1];

FSppZ[x1_,x2_]/; x1=!=x2 :=x2/(2 (x2-1))-(x1^2 Log[x1])/(2 (x1-1) (x1-x2))+(x2 (x1 (-2+x2)+x2) Log[x2])/(2 (x1-x2) (-1+x2)^2);

FSppZ[x1_,x1_]:= 0;


FdS[x_, {f_, s_}] /; x =!= 0 := 
  Q[f] ((-3 + x)/(4 (-1 + x)^2) + Log[x]/(2 (-1 + x)^3)) + 
   Q[s] ((1 + x)/(4 (-1 + x)^2) - (x Log[x])/(2 (-1 + x)^3));

FdSp[x_, {f_, s_}] /; x =!= 0 := 
  Q[f] ((-2 - 5 x + x^2)/(24 (-1 + x)^3) + (x Log[x])/(
      4 (-1 + x)^4)) + 
   Q[s] ((-1 + 5 x + 2 x^2)/(24 (-1 + x)^3) - (x^2 Log[x])/(
      4 (-1 + x)^4));

fdV[x_, {f_, v_}] := 
  Q[f] ((1 - 5 x - 2 x^2)/(8 (-1 + x)^3) + (3 x^2 Log[x])/(
      4 (-1 + x)^4)) + 
   Q[v] ((2 - 7 x + 11 x^2)/(8 (-1 + x)^3) - (3 x^3 Log[x])/(
      4 (-1 + x)^4));

FdV[x0_, x_, {f_, v_}] /; x0 =!= 0 && x =!= 0 := 
  fdV[x, {f, v}] - fdV[x0, {f, v}];

FdV[0, x_, {f_, v_}] := 
  fdV[x, {f, v}] - Limit[fdV[x0, {f, v}], x0 -> 0];

FdVp[x_, {f_, v_}] /; x =!= 0 := 
  Q[f] ((4 + x + x^2)/(4 (-1 + x)^2) - (3 x Log[x])/(2 (-1 + x)^3)) - 
   Q[v] (-((4 - 11 x + x^2)/(4 (-1 + x)^2)) - (3 x^2 Log[x])/(
      2 (-1 + x)^3));

End[]

EndPackage[]

wellPutInfo["Startup"];

SetOptions[getc,
	   "Leptons" -> {{\[Nu],0,0}},
	   "Quarks" -> {{u,0,2/3},{t, Null, 2/3}},
	   "Scalars" -> {{}},
	   "ZBosons" -> {{Z, Null, 0}},
	   "Vectors" -> {{W, Null, 1}},
	   "Externals" -> {{s,0,-1/3},{b,mb,-1/3},{\[Nu],0,0},{\[Mu],0,-1}},
	   CouplingRules -> {
	     g["vff",L,W,t,b]->-(el/Sqrt[2]/SW) Vtb,
	     g["vff",L,bar[W],s,t]->-(el/Sqrt[2]/SW) Vtsc,
	     g["vff",L,bar[W],e|\[Mu],\[Nu]]->-(el/Sqrt[2]/SW),
	     g["vff",L,W,\[Nu],e|\[Mu]]->-(el/Sqrt[2]/SW),
	     g["vff",R,W,__]->0,
	     g["vff",_,Z,a_,b_]:>0/;a=!=b} ];

getc[{L, L}, {s, b, \[Mu]}, {PerformSum -> False}];
