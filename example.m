(* ::Package:: *)

(* ::Chapter:: *)
(*Weak Effective Lagrangian Library for Perturbative Unitary Theories*)


(*$DirectDMDirectory="<path/to/folder/containing/WellPut>";*)
$WellPutDirectory=NotebookDirectory[];
AppendTo[$Path,$WellPutDirectory];


<<WellPut`


(* ::Section::Closed:: *)
(*List of Available Functions and Examples *)


(* ::Input:: *)
(*wellPutInfo[];*)


(* ::Subsubsection:: *)
(*Gauge Independent Loop Functions*)


(* ::Text:: *)
(*There are two \[Gamma]-Penguin contributions. The one that incorporates internal vector bosons is combined with the Z-Bosons,*)
(*while the other \[Gamma]-Penguin involves only internal scalars and fermions. *)


photonFunctions={F["VAZ",Subscript[x, 0],x,{f,v}],F["SA",x,{f,s}]};


(* ::Text:: *)
(*The loop - function depends on the mass ratios of  the internal particles, and in case of the combined *)
(*\[Gamma] - Z - Penguin also on the mass of the fermion that has been eliminated using a GIM like mechanism.*)
(*In addition there is an implicit dependence on the charge of the internal fermion f and the internal vector / scalar (v / s).*)
(*The explicit from can be found using the "replaceFunctions" function:*)


replaceFunctions[photonFunctions]


(* ::Text:: *)
(*They can be formatted in a nicer way using*)


TraditionalForm[photonFunctions]


(* ::Text:: *)
(*The combined contributions of boxes and Z-Bosons are functions of the chirality (and *)


(* ::Text:: *)
(*The remaining functions only depend on the mass ratios. They are:*)


remainingFunctions={
F["VBZ",L,Subscript[x, 0],x,y,z],F["VB'Z",L,Subscript[x, 0],x,y,z],
F["VBZ",R,Subscript[x, 0],x,y,z],F["VB'Z",R,Subscript[x, 0],x,y,z],
F["VZ",x,y],F["V'Z",x,y],F["V''Z",Subscript[x, 0],x,y],
F["VSB",x,y,z],F["VSZ",x,y],F["VS'Z",x,y],
F["SB",x,y,z],F["SZ",x,y],F["S'Z",x,y],F["S''Z",x,y]};


TraditionalForm[remainingFunctions]


replaceFunctions[remainingFunctions]


(* ::Section::Closed:: *)
(*Standard Model*)


Format[Vtsc,TraditionalForm]:=Subsuperscript[V,ts,"*"];
Format[Vtb,TraditionalForm]:=Subscript[V,tb];


(* ::Text:: *)
(*In the Standard Model we only have two contributions for C9 given by two gauge independent functions:*)


Expand@getc[9,{s,b,\[Mu]}]


(* ::Text:: *)
(*Here we already substituted the standard model couplings via the "CoupingRules" option of get c. *)
(*Switching the replacement off and using traditional form:*)


TraditionalForm@Expand@getc[9,{s,b,\[Mu]},ReplaceCouplings->False]


(* ::Text:: *)
(*Note that above we have used the fact that the Z couples diagonally. Switching  this off we find:*)


TraditionalForm@Expand@getc[9,{s,b,\[Mu]},ReplaceZeroCouplings->False,
"ReplaceCouplings"->False]


(* ::Text:: *)
(*The analytic functions are given via the "ReplaceFunctions" option.  Note that we also replace the charges appearing in the*)
(*loop functions that involve photon couplings, since the "ReplaceCharges" option is set to true.*)


SetOptions[getc,ReplaceFunctions->True];
res7=Collect[getc[7,{s,b}],Log[_],Simplify]/.{SW->Sin[W]};
res9=Collect[getc[9,{s,b,\[Mu]}],{SW,Log[_]},Simplify]/.{SW->Sin[W]};
res10=Collect[getc[10,{s,b,\[Mu]}],{SW,Log[_]},Simplify]/.{SW->Sin[W]};
SetOptions[getc,ReplaceFunctions->False];


TraditionalForm[TableForm[{res7,res9,res10}]]


(* ::Section::Closed:: *)
(*2HDM Type II*)


<<WellPut`


(* ::Text:: *)
(*Let us only look at the contribution of charged scalars to C_7 in the 2HDM type II*)


THDMCouplings = {
g["sff",L,Hp,t,b] ->m[t]/tb Vtb el/(SW Sqrt[2] m[W]),
g["sff",R,Hp,t,b] ->m[b]tb Vtb el/(SW Sqrt[2] m[W]),
g["sff",L,bar[Hp],s,t]->m[s]tb Vtsc el/(SW Sqrt[2]m [W]),
g["sff",R,bar[Hp],s,t]->m[t]/tb Vtsc el/(SW Sqrt[2]m [W]),
g["sff",L|R,Hp,\[Nu],\[Mu]] ->0,
g["sff",L,Hp,u,b] ->0,
g["sff",R,bar[Hp],s,u]->0,
g["vff",L,Z,t,t]->(-1/2+2/3 SW^2)(el/(SW CW)),g["vff",L,Z,b,b]->(1/2-1/3 SW^2)(el/(SW CW)),
g["vff",L,Z,\[Mu],\[Mu]]->(1/2-SW^2)(el/(SW CW)),
g["vff",L,Z,\[Nu],\[Nu]]->(-1/2)(el/(SW CW)),
g["vff",R,Z,t,t]->(2/3SW^2)(el/(SW CW)),
g["vff",R,Z,\[Mu],\[Mu]]->(-SW^2)(el/(SW CW)),
g["vss",Z,Hp,bar[Hp]]-> (-1/2+SW^2)  (el/(SW CW))};


SetOptions[getc,{
Leptons->{{\[Nu],0,0}},
Scalars->{{Hp,Null,1}},
Vectors->{{}},
Quarks->{{u,0,2/3},{t,Null,2/3}},
Externals->{{s,0,-1/3},{b,mb,-1/3},{\[Nu],0,0},{\[Mu],0,-1}},
ZBosons->{{Z,Null,0}},
CouplingRules->THDMCouplings}];


norm=y (m[Hp]^2 SW^2 2)/(m[t]^2 el^2Vtb Vtsc m[b]);
Collect[Expand[getc[7,{s,b},ReplaceFunctions->True]norm]/.x[__]->y,tb,
	Collect[#,Log[_],Simplify]&]


(* ::Section::Closed:: *)
(*Standard Model with 3 Up Quarks*)


<<WellPut`


(* ::Text:: *)
(*When we also include a massless charm quark in our calculation, we still get the same result for the top contribution:*)


Format[Vtsc,TraditionalForm]:=Subsuperscript[V,ts,"*"];
Format[Vtb,TraditionalForm]:=Subscript[V,tb];getc[9,{s,b,\[Mu]},
	Quarks->{{u,0,2/3},{c,0,2/3},{t,Null,2/3}}]//TraditionalForm


(* ::Section::Closed:: *)
(*U(1)' Model*)


<<WellPut`


(* ::Subsection:: *)
(*Define Contributing Particles and Coupling Constants*)


SetOptions[getc,
{Quarks->{{u,0,2/3},{t,Null,2/3},{T,Null,2/3}},
ZBosons->{{Z,,0},{Zp,,0}},
CouplingRules->U1pcouplings}];


mix[L,t]=Cos[l];
mix[L,T]=Sin[l];
mix[R,t]=Cos[r];
mix[R,T]=Sin[r];
U1pcouplings=
{g["vff",L,W,t:(t|T),b]:>-mix[L,t](el/Sqrt[2]/SW) Vtb,
g["vff",L,bar[W],s,t:(t|T)]:>-mix[L,t] (el/Sqrt[2]/SW) Vtsc,

g["vff",L,Z,t,T]->-(el/(2CW SW))Cos[l]Sin[l],
g["vff",L,Z,T,t]->-(el/(2CW SW))Cos[l]Sin[l],

g["vff",c_,Zp,t1:(t|T),t2:(t|T)]:>gt qp mix[c,t1]mix[c,t2],

g["vff",L,bar[W],e|\[Mu],\[Nu]]->-(el/Sqrt[2]/SW),
g["vff",L,W,\[Nu],e|\[Mu]]->-(el/Sqrt[2]/SW),

g["vff",L,Z,\[Mu],\[Mu]]->el (-1/2+SW^2)/(CW SW),
g["vff",R,Z,\[Mu],\[Mu]]->el (SW^2)/(CW SW),

g["vff",L,Zp,\[Mu],\[Mu]]->-gt (qlpv-qlpa),
g["vff",R,Zp,\[Mu],\[Mu]]->-gt (qlpv+qlpa),

g["vff",R,W,__]->0,
g["vff",L,Z|Zp,u,t|T]->0,
g["vff",L,Z|Zp,t|T,u]->0,
g["vff",R,Z,a_,b_]:>0/;a=!=b}/.
{SW->Sin[w],CW->Cos[w]};


Format[Vtsc,TraditionalForm]:=Subsuperscript[V,ts,"*"];
Format[qp,TraditionalForm]:=q';
Format[qlpa,TraditionalForm]:=Subscript[q',10];
Format[qlpv,TraditionalForm]:=Subscript[q',9];
Format[el,TraditionalForm]:=e;
Format[Vtb,TraditionalForm]:=Subscript[V,tb];
Format[gt,TraditionalForm]:=\!\(\*OverscriptBox[\(g\), \(~\)]\);


(* ::Subsection:: *)
(*Write Wilson Coefficients in terms of Gauge Independent Functions*)


deltac10=Expand[getc[10,{s,b,\[Mu]}]-(getc[10,{s,b,\[Mu]}]/.l->0)];
deltac9=Expand[getc[9,{s,b,\[Mu]}]-(getc[9,{s,b,\[Mu]}]/.l->0)];


Collect[deltac9,F[__],Simplify]//TraditionalForm


Collect[deltac10,F[__],Simplify]//TraditionalForm


(* ::Subsection:: *)
(*Study the decoupling behaviour*)


(* ::Text:: *)
(*In the limit of a large T mass m_T we can expand in the small parameter \[Epsilon] = m_t / m_T:*)


decoupling={l->Tan[r]\[Epsilon],x[t,W]->x,x[T,W]->x/\[Epsilon]^2};


(* ::Text:: *)
(*Next we use the loop functions in the new physics contributions to c9 and c10 and apply the above decoupling relations:*)


deltac9=Expand[getc[9,{s,b,\[Mu]},ReplaceFunctions->True]-\
	(getc[9,{s,b,\[Mu]},ReplaceFunctions->True]/.l->0)];
deltac9=deltac9/.decoupling/.Log[a_/\[Epsilon]^2]->Log[a]-2 Log[\[Epsilon]];
deltac10=Expand[getc[10,{s,b,\[Mu]},ReplaceFunctions->True]-\
	(getc[10,{s,b,\[Mu]},ReplaceFunctions->True]/.l->0)];
deltac10=deltac10/.decoupling/.Log[a_/\[Epsilon]^2]->Log[a]-2 Log[\[Epsilon]];


(* ::Text:: *)
(*One can also look at c7 even though the effect will be suppressed by \[Epsilon]^2*)


deltac7=Expand[getc[7,{s,b},ReplaceFunctions->True]-
(getc[7,{s,b},ReplaceFunctions->True]/.l->0)];
deltac7=deltac7/.decoupling/.Log[a_/\[Epsilon]^2]->Log[a]-2 Log[\[Epsilon]];


(* ::Text:: *)
(*The new physics contributions of c9 has three terms that are suppressed by the square of MZ', MW and MZ respectively.*)


{c9zp,c9w,c9z}=List@@Collect[deltac9,m[_]]/.m[Z]->m[W]/Cos[w];


(* ::Text:: *)
(*The latter two decouple when we send mT -> infinity:*)


Collect[Normal[Series[c9w+c9z,{\[Epsilon],0,0}]],Log[_],Simplify[#,{x>0,\[Epsilon]>0}]&]/.U1pcouplings


(* ::Text:: *)
(*While the Z' term gives the dominant new physics contribution in the limit mT -> infinity*)


exc9=Collect[Normal[Series[c9zp,{\[Epsilon],0,0}]],Log[_],
	FullSimplify[#,{x>0,\[Epsilon]>0}]&]/.\[Epsilon]->Subscript[m,t]/Subscript[m,T]


(* ::Text:: *)
(*Similiar for c10 and using that MW = MZ cos(W).*)


exc10=Collect[Normal[Series[deltac10/.m[Z]->m[W]/Cos[w],{\[Epsilon],0,0}]],Log[_],
	FullSimplify[#,{x>0,\[Epsilon]>0}]&]/.\[Epsilon]->Subscript[m,t]/Subscript[m,T]


normc9=x gt^2 el^2 qp qlpv (Vtb Vtsc)Sin[r]^2/(Sin[w]m[Zp])^2/2;
normc10=x gt^2 el^2 qp qlpa (Vtb Vtsc)Sin[r]^2/(Sin[w]m[Zp])^2/2;


Collect[FullSimplify[exc9/normc9],{Sec[_],Log[_]},FullSimplify]//TraditionalForm


Collect[FullSimplify[exc10/normc10],{Sec[_],Log[_]},FullSimplify]//TraditionalForm


exc7=Collect[Normal[Series[deltac7/.m[Z]->m[W]/Cos[w],{\[Epsilon],0,2}]],Log[_],FullSimplify[#,{x>0,\[Epsilon]>0}]&]/.\[Epsilon]->Subscript[m,t]/Subscript[m,T]


(* ::Section::Closed:: *)
(*U(1)L\[Mu] - L\[Tau] Model*)


<<WellPut`


(* ::Text:: *)
(*Note that we have the charge conjugated d~ inside the loop*)


(* ::Subsection:: *)
(*Define Contributing Particles and Coupling Constants*)


U1mutauCouplings = {
g["vff",_,Zp,(Np|Nm),(Np|Nm)]->-gx Qn/2,
g["sff",R,bar[dtc],d:(s|b),Nm_]->-1/Sqrt[2] y[d],
g["sff",L,dtc,Nm_,d:(s|b)]->-1/Sqrt[2] Conjugate[y[d]],
g["sff",L,bar[dtc],d:(s|b),Nm_]->0,
g["sff",R,dtc,Nm_,d:(s|b)]->0,
g["vff",L_,Z,Nm|Np,_]->0,
g["vss",Zp,dtc,bar[dtc]]->-gx Qn,
g["vff",L_, Zp, \[Mu], \[Mu]]->-1 gx, 
g["vss",Z,dtc,bar[dtc]]->-(1-2/3 SW^2)el/(2SW CW),
g["vff",L,Z,b,b]->(1-2/3 SW^2)(el/(2SW CW)),
g["vff",L,Z,\[Mu],\[Mu]]->(1-2 SW^2)(el/(2SW CW)),
g["vff",R,Z,\[Mu],\[Mu]]->-(SW^2)(el/(SW CW)),
g["vff", L_, Zp, b, b]-> 0 (*Qn gx*),
g["sff", R_, bar[dtc], \[Mu], \[Nu]]->0, (*dtc-colored scalar*)
g["sff", R_, dtc, \[Mu], \[Nu]]->0};


SetOptions[getc,{
Leptons->{{\[Nu],0,0}},
Scalars->{{dtc,Null,1/3}},
Vectors->{{}},
Quarks->{{Np, Null,0},{Nm,Null,0}},
Externals->{{s,0,-1/3},{b,mb,-1/3},{\[Nu],0,0},{\[Mu],0,-1}},
ZBosons->{{Z,Null,0},{Zp,Null, 0}},
CouplingRules->U1mutauCouplings}];


Format[Zp,TraditionalForm]:=Z';
Format[Nm,TraditionalForm]:=Subscript[N,"-"];
Format[Np,TraditionalForm]:=Subscript[N,"+"];
Format[dtc,TraditionalForm]:=Superscript[\!\(\*OverscriptBox[\(d\), \(~\)]\),"*"];


(* ::Subsection:: *)
(*Write Wilson Coefficients in terms of Gauge Independent Functions*)


resc9=Collect[getc[9,{s,b,\[Mu]}],m[_],Simplify]//TraditionalForm


resc10=Collect[getc[10,{s,b,\[Mu]}],m[_],Simplify]//TraditionalForm


resc7=Collect[getc[7,{s,b}],m[_],Simplify]//TraditionalForm


(* ::Subsection:: *)
(*Result Baek*)


resc9=Collect[getc[9,{s,b,\[Mu]},ReplaceFunctions->True],m[_],Simplify];
resc9=Coefficient[resc9,m[Zp],-2]/.{x[Nm,dtc]->x1,x[Np,dtc]->x2};


resc7=Collect[getc[7,{s,b},ReplaceFunctions->True],m[_],Simplify]\
	/.{x[Nm,dtc]->x1,x[Np,dtc]->x2};


j[x_]:=(x Log[x])/(x-1)
k[x_]:=(x^2Log[x])/(x-1)
Y[x_]:=((-3x^2+4x-1+2x^2 Log[x])/(8(x-1)^2)) 
J1[x_]:=(1-6x+3x^2+2x^3-6x^2 Log[x])/(12(1-x)^4);


j[x1_,x2_]:=(j[x1]-j[x2])/(x1-x2) 
k[x1_,x2_]:=((k[x1]-k[x2])/(x1-x2))  


knownResC7[xm_,xp_]:=J1[xm]+J1[xp];


knownResC9[xm_,xp_]:=1/4+Sqrt[xm xp]j[xm,xp]-1/2 k[xm,xp]+Y[xm]+Y[xp]


Collect[1/2 Conjugate[y[b]]y[s] Qn gx^2knownResC9[x2,x1]-resc9,Log[_],FullSimplify[#,{x2>x1>0}]&]


((Conjugate[y[b]] y[s] m[b])/(12 m[dtc]^2 ) knownResC7[x1,x2]) - resc7//Simplify


(* ::Section::Closed:: *)
(*Model with vector like fermions and neutral scalars*)


<<WellPut`


(* ::Subsection:: *)
(*Define Contributing Particles and Coupling Constants*)


Format[\[Psi]l]:=Subscript[\[Psi], L];
Format[\[Phi]l]:=Subscript[\[Phi], L];
Format[\[Phi]h]:=Subscript[\[Phi], H];
Format[\[Psi]q]:=Subscript[\[Psi], Q];


c[hc[a_]]:=a;
Format[hc[arg_]]:=SuperStar[arg];
ModelCouplings={
g["sff",R,\[CapitalPhi]:(\[Phi]l|\[Phi]h),f:(b),\[Psi]_]->
hc[g["sff",L,bar[\[CapitalPhi]],\[Psi],f]],
g["sff",L,(\[Phi]l|\[Phi]h),_,f2_]->0,
g["sff",R,bar[\[CapitalPhi]_],_,(b|s|\[Mu])]->0,
g["sff",R,\[CapitalPhi]_,_,(b|s|\[Mu])]->0,
g["vff",L,Z,\[Psi]q,\[Psi]q]->(1-2/3 sW^2)(g/(2cW)),
g["vff",R,Z,\[Psi]q,\[Psi]q]->(1-2/3 sW^2)(g/(2cW)),
g["vff",L,Z,b,b]->(1-2/3 sW^2)(g/(2cW)),
g["vff",L,Z,\[Mu],\[Mu]]->(1-2 sW^2)(g/(2cW)),
g["vff",R,Z,\[Mu],\[Mu]]->-(sW^2)(g/cW)};


SetOptions[getc,{
Leptons->{{\[Psi]l,Null,-1}},
Scalars->{{bar[\[Phi]l],Null,0},{bar[\[Phi]h],Null,0}},
Vectors->{{}},
Quarks->{{\[Psi]q, Null,-1/3}},
Externals->{{s,0,-1/3},{b,mb,-1/3},{\[Nu],0,0},{\[Mu],0,-1}},
ZBosons->{{Z,Null,0}},
CouplingRules->ModelCouplings}];


(* ::Subsection:: *)
(*Compare with Grinstein, Pokorski & Ross *)


myc9=getc[9,{s,b,\[Mu]}];
myc10=getc[10,{s,b,\[Mu]}];


(* ::Text:: *)
(*Keep the box contributions only*)


c9b=myc9/.F[a_,b__]:>0/;StringMatchQ[a,RegularExpression["S.*[ZA]"]];
c10b=myc10/.F[a_,b__]:>0/;StringMatchQ[a,RegularExpression["S.*[ZA]"]];


(* ::Text:: *)
(*The replacement rule for F["SB",...] ->G[...]*)


Collect[(4/y\[Phi]l)replaceFunctions@F["SB",x[\[Psi]q,\[Phi]l],x[\[Phi]l,\[Phi]h],x[\[Psi]l,\[Phi]l]]/.{
	x[\[Psi]l,\[Phi]l]->m\[Psi]l^2/m\[Phi]l^2,x[\[Psi]q,\[Phi]l]->m\[Psi]q^2/m\[Phi]l^2,
	x[\[Psi]l,\[Phi]h]->m\[Psi]l^2/m\[Phi]h^2,x[\[Psi]q,\[Phi]h]->m\[Psi]q^2/m\[Phi]h^2,
	x[\[Phi]l,\[Phi]h]->m\[Phi]l^2/m\[Phi]h^2,x[\[Phi]h,\[Phi]l]->m\[Phi]h^2/m\[Phi]l^2}/.{Log[a_^2/b_^2]:>2Log[a]-2Log[b]}/.{
Log[m\[Phi]l]->1/2Log[y\[Phi]l]+Log[m\[Psi]q],Log[m\[Phi]h]->1/2Log[y\[Phi]h]+Log[m\[Psi]q],Log[m\[Psi]l]->1/2Log[y\[Psi]l]+Log[m\[Psi]q]
}/.{m\[Phi]l->Sqrt[y\[Phi]l]*m\[Psi]q,m\[Phi]h->Sqrt[y\[Phi]h]*m\[Psi]q,m\[Psi]l->Sqrt[y\[Psi]l]m\[Psi]q},
_Log,FullSimplify[#,Assumptions->{y\[Phi]l>0,y\[Psi]l>0,y\[Phi]h>0,m\[Psi]q>0}]&]


(* ::Text:: *)
(*This gives: F["SB",x[\[Psi]q,\[Phi]l],x[\[Phi]l,\[Phi]h],x[\[Psi]l,\[Phi]l]] / m\[Phi]l^2 -> 4 * G[yl,y\[Phi]L,y\[Phi]H] / m\[Psi]q^2*)


c9c=c9b/.{
F["SB",x[\[Psi]q,\[Phi]l],x[\[Phi]l,\[Phi]h],x[\[Psi]l,\[Phi]l]]:>m[bar[\[Phi]l]]^2/m[\[Psi]q]^2/4G[yl,y\[Phi]L,y\[Phi]H],
F["SB",x[\[Psi]q,\[Phi]h],x[\[Phi]h,\[Phi]l],x[\[Psi]l,\[Phi]h]]:>m[bar[\[Phi]h]]^2/m[\[Psi]q]^2/4G[yl,y\[Phi]L,y\[Phi]H],
F["SB",x[\[Psi]q,\[Phi]l],1,x[\[Psi]l,\[Phi]l]]:>m[bar[\[Phi]l]]^2/m[\[Psi]q]^2/4G[yl,y\[Phi]L,y\[Phi]L],
F["SB",x[\[Psi]q,\[Phi]h],1,x[\[Psi]l,\[Phi]h]]:>m[bar[\[Phi]h]]^2/m[\[Psi]q]^2/4G[yl,y\[Phi]H,y\[Phi]H],
g["sff",R,\[Phi]l,\[Mu],\[Psi]l]:>Subscript[\[CapitalGamma],\[Mu]]*VE2\[Mu]\[Conjugate],
hc[g["sff",R,\[Phi]l,\[Mu],\[Psi]l]]:>Subscript[\[CapitalGamma],\[Mu]]\[Conjugate]*VE2\[Mu],
g["sff",R,\[Phi]h,\[Mu],\[Psi]l]:>-Subscript[\[CapitalGamma],\[Mu]]*VE2\[Mu]\[Conjugate],
hc[g["sff",R,\[Phi]h,\[Mu],\[Psi]l]]:>-Subscript[\[CapitalGamma],\[Mu]]\[Conjugate]*VE2\[Mu],
g["sff",R,\[Phi]l,s,\[Psi]q]->Subscript[\[CapitalGamma],b]*V3s\[Conjugate]+Subscript[OverTilde[\[CapitalGamma]],s]V2s\[Conjugate],
g["sff",R,\[Phi]h,s,\[Psi]q]->sign (Subscript[\[CapitalGamma],b]*V3s\[Conjugate]-Subscript[OverTilde[\[CapitalGamma]],s]V2s\[Conjugate]),
g["sff",L,bar[\[Phi]l],\[Psi]q,b]->Subscript[\[CapitalGamma],b]\[Conjugate]*V3b+Subscript[OverTilde[\[CapitalGamma]],s]\[Conjugate]V2b,
g["sff",L,bar[\[Phi]h],\[Psi]q,b]->sign (Subscript[\[CapitalGamma],b]\[Conjugate]*V3b-Subscript[OverTilde[\[CapitalGamma]],s]\[Conjugate]V2b)
}//Simplify


Collect[4\[Pi]/(\[Alpha]*16\[Pi]^2)c9c/.Conjugate[a_]:>SuperStar[a]/.sign->-1,_G,FullSimplify]/.{
VE2\[Mu]*SuperStar[VE2\[Mu]]->Abs[VE2\[Mu]]^2
}
