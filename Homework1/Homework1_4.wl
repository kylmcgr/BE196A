(* ::Package:: *)

(* ::Section:: *)
(*Define CRN notation and data types*)


<<Notation`


Notation[\!\(\*
TagBox[
RowBox[{"r_", 
OverscriptBox["\[LongRightArrow]", "k_"], "p_"}],
"NotationTemplateTag"]\) \[DoubleLongLeftArrow] \!\(\*
TagBox[
RowBox[{"rxn", "[", 
RowBox[{"r_", ",", "p_", ",", "k_"}], "]"}],
"NotationTemplateTag"]\)]


rxn::about="Represents an irreversible reaction. eg. rxn[a+b,c,1]";
revrxn::about="Represents a reversible reaction. eg. revrxn[a+b,c,1,1]";
conc::about="Initial concentration: conc[{x,y},10]";


revrxn[r_,p_,k1_,k2_]:=Sequence[rxn[r,p,k1],rxn[p,r,k2]]


SetAttributes[{gate,rxn,revrxn}, HoldAll]


Seq:=Sequence (*to use instead of /.List->Sequence in definition of gates*)


conc[xs_List,c_]:=(conc[#,c]&/@xs)/.List->Sequence


(* ::Section:: *)
(*Define CRN functions*)


(* ::Text:: *)
(*Converts a system of gates with conc statements to NDSolve statement:*)


gatesysToNDSolve[gatesys_,endtime_]:=
Module[{
inputspecs=Cases[gatesys,conc[x_,c_]:>{x,c}],rsys=DeleteCases[gatesys,conc[___]]},
rxnsysToNDSolve[rsys,inputspecs,endtime]]


(* ::Text:: *)
(*Gets the set of all species for a reaction system:*)


spcsInRxnsys[rsys_]:=Cases[Cases[rsys,rxn[r_,p_,_]:>Seq[r,p]]/.Times|Plus->Seq,_Symbol|_Symbol[__]]//Union


(* ::Text:: *)
(*Converts a reaction system to a list of odes:*)


rxnsysToOdesys[rsys_]:=
Module[{spcs=spcsInRxnsys[rsys], rrates, spccoeffs,odes},
(*create terms for the rates of each reaction*)
rrates = rsys/.rxn[r_,_,k_]:>(k (r/.{Times[b_,s_]:>s^b,Plus->Times}));
(*for each species, get list of net coefficient for each reaction*)
spccoeffs=Coefficient[rsys/.rxn[r_,p_,_]:>p-r,#]& /@ spcs;
(*create ode for each species*)
odes=MapThread[Function[{spc,coeffs},HoldForm[D[spc,t]]==Total[coeffs*rrates]],{spcs, spccoeffs}];
(*change all species conc to be functions of t*)
odes/.s_/;MemberQ[spcs,s]:>s[t]]//ReleaseHold;


(* ::Text:: *)
(*Converts a reaction system to an NDSolve statement, and solves:*)


rxnsysToNDSolve[rsys_,inputspecs_, endtime_]:=
Module[{spcs=spcsInRxnsys[rsys],odesys=rxnsysToOdesys[rsys], NDSolveEqns},
NDSolveEqns = Join[odesys,(#[0]==Cases[inputspecs,{#,c0_}:>c0]/.{{c0_}:> c0,{}->0})&/@spcs];
Quiet[NDSolve[NDSolveEqns, spcs, {t,0,endtime},
AccuracyGoal->MachinePrecision,PrecisionGoal->MachinePrecision,WorkingPrecision->MachinePrecision,MaxSteps->Infinity],{NDSolve::"precw"}][[1]]]


(* ::Text:: *)
(*Generates list of species matching formal pattern:*)


spcsPattern[rsys_,pattern_]:=Cases[spcsInRxnsys[rsys],pattern]


(* ::Text:: *)
(*Generates list of species matching a pattern (eg "g$*" to match all implicit gates; can also do RegularExpression["o..d.\$.*"]):*)


spcsPatternString[rsys_,pattern_]:=Select[spcsInRxnsys[rsys],StringMatchQ[ToString[#],pattern]&]


(* ::Text:: *)
(*Generates set of initial values for given list of species:*)


initsSpcs[spcs_,value_]:={#,value}&/@spcs


(* ::Text:: *)
(*Generates set of initial values for species matching formal pattern:*)


initsSpcsPattern[rsys_,pattern_,value_]:={#,value}&/@spcsPattern[rsys,pattern]


(* ::Text:: *)
(*Generates set of initial values for species matching string pattern:*)


initsSpcsPatternString[rsys_,pattern_,value_]:={#,value}&/@spcsPatternString[rsys,pattern]


(* ::Section:: *)
(*Define seesaw parameters*)


(* ::Text:: *)
(*Rate constants:*)


kf = 2*10^6; (* fast strand displacement rate with extended toehold, unit: /M/s *)
ks = 5*10^4; (* slow strand displacement rate with universal toehold, unit: /M/s *)
krf = 26; (* fast strand dissociation rate with universal toehold, unit: /s *)
krs = 1.3; (* slow strand dissociation rate with extended toehold, unit: /s *)
kl = 10; (* strand displacement leak rate, unit: /M/s *)


(* ::Text:: *)
(*Standard concentration 1x:*)


c = 50*10^-9; (* unit: M *)


(* ::Section:: *)
(*Define seesaw functions*)


(* ::Text:: *)
(*Translates a seesaw gate into a list of reactions:*)


seesaw[x_,l_List,r_List]:={
(* Toehold exchange reactions *)
Outer[revrxn[w[#1,x]+g[x,w[x,#2]],g[w[#1,x],x]+w[x,#2],ks,ks]&,l,r],
(* Thresholding reactions *)
Outer[rxn[w[#,x]+th[w[#,x],x],waste,kf]&,l],
Outer[rxn[w[x,#]+th[x,w[x,#]],waste,kf]&,r],
(* Universal toehold binding reactions *)
Outer[revrxn[g[w[#,x],x]+W,g[w[#,x],x,W],kf,krf]&,l],
Outer[revrxn[g[x,w[x,#]]+W,g[W,x,w[x,#]],kf,krf]&,r],
Outer[revrxn[th[w[#,x],x]+W,th[W,w[#,x],x],kf,krs]&,l],
Outer[revrxn[th[x,w[x,#]]+W,th[x,w[x,#],W],kf,krs]&,r],
Outer[revrxn[w[#,x]+G,w[#,x,G],kf,krf]&,l],
Outer[revrxn[w[#,x]+TH,w[#,x,TH],kf,krs]&,l],
(* Leak reactions *)
Outer[rxn[g[w[#1,x],x]+w[#2,x],g[w[#2,x],x]+w[#1,x],kl]&,l,l],
Outer[rxn[g[x,w[x,#1]]+w[x,#2],g[x,w[x,#2]]+w[x,#1],kl]&,r,r]
}/.List->Sequence


(* ::Text:: *)
(*Translates a reporter into a list of reactions and the initial concentration of the reporter:*)


reporter[x_,l_]:={
rxn[w[l,x]+Rep[x],Fluor[x],ks],
revrxn[Rep[x]+W, Rep[x,W],kf,krf],
revrxn[w[l,x]+G,w[l,x,G],kf,krf],
revrxn[w[l,x]+TH,w[l,x,TH],kf,krs],
conc[Rep[x],1.5*c]
}/.List->Sequence


(* ::Text:: *)
(*Translates logic OR operation into a list of seesaw gates and the concentrations of initial species:*)


seesawOR[x1_,x2_,l_List,r_List]:=
Module[{f},
{seesaw[x1,l,{x2}],
seesaw[x2,{x1},{r,f}],
conc[g[x1,w[x1,x2]],Length[l]*c],
Outer[conc[g[x2,w[x2,#]],1*c]&,r],
conc[w[x2,f],2*Length[r]*c],
conc[th[w[x1,x2],x2],1.1*0.6*c]
}]/.List->Sequence


(* ::Text:: *)
(*Translates logic AND operation into a list of seesaw gates and the concentrations of initial species:*)


seesawAND[x1_,x2_,l_List,r_List]:=
Module[{f},
{seesaw[x1,l,{x2}],
seesaw[x2,{x1},{r,f}],
conc[g[x1,w[x1,x2]],Length[l]*c],
Outer[conc[g[x2,w[x2,#]],1*c]&,r],
conc[w[x2,f],2*Length[r]*c],
conc[th[w[x1,x2],x2],1.1*(Length[l]-1+0.2)*c]
}]/.List->Sequence


(* ::Text:: *)
(*Translates inputs with fan-out into a list of seesaw gates and the concentrations of initial species:*)


inputfanout[x_,l_,r_List]:=
Module[{f},
{seesaw[x,{l},{r,f}],
Outer[conc[g[x,w[x,#]],1*c]&,r],
conc[w[x,f],2*Length[r]*c],
conc[th[w[l,x],x],1.1*0.2*c]
}]/.List->Sequence


(* ::Section:: *)
(*Define input value and simulation time*)


OFF = 0.1;    
ON = 0.9;
inputs=Flatten[Table[{{l},{k},{j},{i}}, {i, 0.1,0.9,0.8}, {j, 0.1,0.9,0.8}, {k, 0.1,0.9,0.8}, {l, 0.1,0.9,0.8}],3];
(*input = {{ON},{ON},{ON},{ON}}*)
time = 10; (* unit: hour *)


(* ::Section:: *)
(*Simulation*)


SIMcircuits=Table[Table[
gatesys={
seesawOR[20,21,{7,11},{38,40}],
seesawAND[22,23,{5,9},{36,42}],
seesawOR[28,29,{13,19},{40}],
seesawAND[30,31,{15,17},{42}],
seesawAND[32,33,{13,17},{48}],
seesawOR[34,35,{15,19},{50}],
seesawAND[36,37,{15,19,23},{46}],
seesawOR[38,39,{13,17,21},{44}],
seesawAND[40,41,{21,29},{46}],
seesawOR[42,43,{23,31},{44}],
seesawAND[44,45,{39,43},{52}],
seesawOR[46,47,{37,41},{54}],
inputfanout[13,12,{28,32,38}],
inputfanout[15,14,{30,34,36}],
inputfanout[17,16,{30,32,38}],
inputfanout[19,18,{28,34,36}],
reporter[48,33],
reporter[50,35],
reporter[52,45],
reporter[54,47],
conc[w[5,22],(1-x1)*c],
conc[w[7,20],x1*c],
conc[w[9,22],(1-x2)*c],
conc[w[11,20],x2*c],
conc[w[12,13],(1-x3)*c],
conc[w[14,15],x3*c],
conc[w[16,17],(1-x4)*c],
conc[w[18,19],x4*c],
conc[W,56*c],
conc[G,58*c],
conc[TH,1.1*12.6*c]
};
sol=gatesysToNDSolve[gatesys,time*60*60]//ReleaseHold//Flatten;
{Fluor[48][t*60*60]/c/.sol,Fluor[50][t*60*60]/c/.sol,Fluor[52][t*60*60]/c/.sol,Fluor[54][t*60*60]/c/.sol},
{x1,inputs[[i]][[1]]},{x2,inputs[[i]][[2]]},{x3,inputs[[i]][[3]]},{x4,inputs[[i]][[4]]}
],{i, 1, Length[inputs]}];


(* ::Section:: *)
(*Plot*)


plots=Table[Plot[Evaluate[SIMcircuits[[i]]//Flatten],{t,0,time},
Frame->True,FrameLabel->{"Time (hours)","Output"},
PlotStyle->{Thick},
PlotLegends->LineLegend[Automatic,{"y2^0","y2^1","y1^0","y1^1"}],
LabelStyle->Directive[FontSize->12,FontFamily->"Helvetica"],
PlotLabel->Flatten[Reverse[inputs[[i]],1]]==i-1,
PlotRange->All,AspectRatio->1/1.32,ImageSize->250],{i, 1, Length[SIMcircuits]}]


(* ::Section:: *)
(*The output values are given by the different curves. y2^0 and y1^0 (blue and green curves) indicate a 0 at that binary position and y2^1 and y1^1 (yellow and red curves) indicate a 1 at that binary position. If the curve increases in the graph, then that number is present at that position. The output is calculated from the binary number y2y1. So yellow and green lines that increase give the binary number 10 which is 2.*)
