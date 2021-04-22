(* ::Package:: *)


BeginPackage["WellPut`"]


Q::usage = "Q[a] represents the charge of particle a in units of the electric charge";

el::usage = "el represents the electric coupling constants";

wellPutInfo::usage =
"Describes all available functions";

getc::usage =
"getc[op_,ef_] generates the Wilson coefficient of the operator \"op\" for the external fields \"ef\".
Currently coefficients for dipole and current-current operators can be generated.

For current-current operators
  - op = 9 gives C_9
  - op = 10 gives C_10
  - op = {\[Sigma]1,\[Sigma]2} gives C_\[Sigma]1\[Sigma]2,
         where \[Sigma]1 and \[Sigma]2 denote the chirality of the quark and lepton current respectively.
  - ef = {di,dj,mu} specifies the flavour of the incoming quark dj, the outgoing quark di and the lepton pair.
         The symbol used should be members of the particle identifier of the Externals defined below.

For dipole operators
  - op = 7 gives D_L (or C_7\[Gamma])
  - op = \[Sigma] gives D_\[Sigma].
  - ef = {di,dj} specifies the flavour of the incoming quark dj, the outgoing quark di.
         The symbols used should be members of the particle identifier of the Externals defined below.

Possible Options are:
  PerformSum -> True, 
  ReplaceCharges -> True,
  ReplaceCouplings -> True,
  ReplaceFunctions -> False,
  ReplaceZeroCouplings -> True,
  CouplingRules -> {g[vff,L,W,t,b]->-(el/Sqrt[2]/SW) Vtb, ... },
  Leptons -> {{\[Nu],0,0}}, 
  Quarks -> {{u,0,2/3},{t, Null, 2/3}},
  Scalars -> {{}},
  ZBosons -> {{Z, Null, 0}},
  Vectors -> {{W, Null, 1}},
  Externals -> {{s,0,-1/3},{b,mb,-1/3},{\[Nu],0,0},{\[Mu],0,-1}}

Leptons, Quarks, Scalars, ZBosons and Vectors specify the particles that appear as
internal particles in the evaluation of the Wilson coefficients. Externals are possible external particles.
All particles are given as a list of the form {{particleid1, massinfo1, charge1}, {particleid2, massinfo2, charge2}, ...}.
E.g. Quarks = {{u, 0, 2/3}, {c, 0, 2/3}, {t, Null, 2/3}} would imply that 3 quarks appear in the sum,
where u, c and t are symbols that identify the particles. The mass of u and c will be treated as zero
in the matching calculation, while the mass of the top quark is kept. u, c and t all have charge 2/3.
For Quarks, a GIM like mechanism will be performed on the first fermion on the list.
The charges of the particles have to such that they contribute in the diagrams in arxiv:0421.???? .

PerformSum performs the sums over internal particles if set to True.

If ReplaceCoupling is set True, the replacement rule specified in CouplingRules
will be applied to the generated expression.

CouplingRules will be replacement rules of generic coupling constants.

If ReplaceZeroCouplings is set True couplings that will be replaced to zero
in the rule of CouplingRules will be removed from the generic expression after
particle summation has been performed.

The general expressions and the loop functions contain charges these will be replaced by the
charges specified in the other options if ReplaceCharges is set to True.";

x::usage = "x[a,b] represents the mass ratio m[a]^2/m[b]^2";

m::usage = "m[a] represents the mass of particle a";

mu::usage = "mu is the scale of dimensional regularisation";

sum::usage "sum[...,{f1, f2, ...}] represents the sum over f1, f2, ... .
If the sum specification starts with f0, a GIM like mechanism is performed";

replaceFunctions::usage =
"replaceFunctions[arg_] replaces all loop functions of the form F[id_,x_,...,particlid___] by their analytitic expression.";

F::usage =
"F[identifier_, massratios__, particleid___List] represents the loop function as defined in the paper ... .
- massratios are a sequence of x_a^b = m[b]^2/m[a]^2.
- particleid is of the form {fermionid, bosonid} of particles inside the loop function.
  This is only present for loopfunctions that have photon couplings (\"VAZ\", \"SA\",
  \"dV\", \"dV'\", \"dS\", \"dS'\"), and the identifiers are used to fix the fermion charges.
- loop functions identifiers for the neutral current oeprators are:
  \"VAZ\", \"VBZ\", \"VB'Z\", \"VZ\", \"V'Z\", \"V''Z\",
  \"VSB\", \"VSZ\", \"VS'Z\",
  \"SA\", \"SB\", \"SZ\", \"S'Z\", \"S''Z\"
- loop functions identifiers for the dipole operators are:
  \"dV\", \"dV'\", \"dS\", \"dS'\"";

g::usage =
"g[id_,(\[Sigma]_,)fields__] represents the coupling of the theory.
 All fields given in the order of the identifier (id), which specifies the type of coupling:
 - \"vff\" is a vector-fermion1-fermion2 coupling, where
   \[Sigma] is the chirality, the vector boson and fermion2 are incoming, fermion1 outgoing.
 - \"sff\" is a scalar-fermion1-fermion2 coupling, where
   \[Sigma] is the chirality, the scalar and fermion2 are incoming, fermion1 outgoing.
 - \"vvv\" is the triple vector coupling and all fields are incoming.
 - \"vvs\" is the coupling of two vectors with one scalar, where all fields are incoming.
 - \"vss\" is the coupling of one vectors with two scalar, where all fields are incoming.";

L::usage = "L represents left handed coupling";
R::usage = "R represents right handed coupling";

bar::usage =
"bar[f] denotes the charge conjugation of field f, while bar[L] = R and bar[R] = L are flipped chiralities";

s1::usage = "represents scalars in the loop";
s2::usage = "represents scalars in the loop";
v1::usage = "represents vector bosons in the loop";
v2::usage = "represents vector bosons in the loop";
z1::usage = "represents neutral vector bosons (Z's) that are not in the loop";
f0::usage = "represents the fermion that is elimenated using GIM like unitarity";
f1::usage = "represents fermions coupling to the external quarks";
f2::usage = "represents fermions coupling to the external quarks";
f3::usage = "represents fermions coupling to the external leptons";


PerformSum::usage = "PerformSum is an option for getc that specifies whether the sum over loop fields is performed.";
ReplaceCharges::usage = "ReplaceCharges is an option for getc that, if set to True, will replace any charges that appear in the loop functions by their values sepcified in CouplingRules.";
ReplaceFunctions::usage = "ReplaceFunctions is an option for getc that replaces the symbolic representation F[id,x,...,particle] by their explicit analytic expressions.";
ReplaceCouplings::usage = "ReplaceCouplings  is an option for getc that specifies the couplings of the new states.";
ReplaceZeroCouplings::usage = "ReplaceZeroCouplings is an option for getc that, if set to True, ...";
CouplingRules::usage = "CouplingRules is an option for getc that specifies the list of couplings of the loop fields.";
Leptons::usage = "Leptons is an option for getc that specifies color-neutral loop fields. It should be a list of the form {{\[Psi],Mass,Charge},...}.";
Quarks::usage = "Leptons is an option for getc that specifies colored loop fields. It should be a list of the form {{Q,Mass,Charge},...}.";
Scalars::usage = "Leptons is an option for getc that specifies color-neutral loop fields. It should be a list of the form {{\[Phi],Mass,Charge},...}.";
ZBosons::usage = "Leptons is an option for getc that specifies neutral spin-1 loop fields. It should be a list of the form {{Z,Mass,Charge},...}.";
Vectors::usage = "Vectors is an option for getc that specifies charged spin-1 loop fields. It should be a list of the form {{W,Mass,Charge},...}.";
Externals::usage = "Leptons is an option for getc that specifies external fields. It should be a list of the form {{Field,Mass,Charge},...}.";

EndPackage[]

