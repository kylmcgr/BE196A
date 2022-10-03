(* ::Package:: *)

(* ::Text:: *)
(*DNA strand displacement (DSD) diagram package is developed by Lulu Qian. *)
(*Version 1.0 is released on 06/20/2015.*)
(*Copyright 2015. *)
(*http://qianlab.caltech.edu/DSDdiagram/*)


(* define domain colors *)
domainColor=Flatten[Table[{
	Table[Lighter[ColorData["DarkRainbow"][i],j],{i,0,3/4,1/4}],
	Table[Lighter[ColorData["DarkTerrain"][i],j],{i,1/4,3/4,1/2}]},
{j,0,0.8,0.2}],2];
domainColor=Append[domainColor,Black]; (* the last color is used for edges with no specific identity *)


Options[diagramDSD]={DNASequence->False,VectorFormat->True,DistinctComplement->True,
LineBasepair->False,GreyBasepair->False,DotBasepair->False,
DomainSeparator->False,SquiggleToehold->False,MaxToeholdLength->8,
ReverseOrient->False,DomainLabelSize->24,DomainLabel->True,LineThickness->3};


diagramDSD[molecule_,OptionsPattern[]]:=
Module[{domains,Nd,label,seq,length,orient,color,style,pos,comp,basepair,hairpin,
pcomp,ko,kd,kc,lc,kp,forward,nick5,nick3,xmin,xmax,ymin,ymax,name,a1,a2,l,r,d,cen,lr,i,j,
DomainObject,TextObject,plotColor},
domains=StringSplit[molecule," "];
Nd=Sum[If[domains[[i]]=="+" || StringTake[domains[[i]],1]=="@",0,1],{i,1,Length[domains]}]; (* number of domains *)
label=Table["",{i,1,Nd}];
seq=Table["",{i,1,Nd}];
length=Table[0,{i,1,Nd}];
orient=Table[0,{i,1,Nd}];
color=Table[-1,{i,1,Nd}];
style=Table[1,{i,1,Nd}]; (* 1: solid line, 2: dashed line, 3: solid line with an arrow, 4: dashed line with an arrow *)
pos=Table[{{0,0},{0,0},{0,0},{0,0}},{i,1,Nd}];
comp=Table[{0,""},{i,1,Nd}];
basepair=Table[0,{i,1,Nd}];
hairpin=Table[0,{i,1,Nd}];
If[OptionValue[DNASequence], pcomp=-6, pcomp=-2-(OptionValue[LineThickness]-3)/10];
ko=0; (* current orientation *)
kd=1; (* current domain *)
kc=0; (* current complementary domain that hasn't been paired *)
lc=0; (* last ( in the current strand *)
kp={0,0}; (* current position *)
forward=1; (* forward mode *)
nick5=0; (* nick at the 5' end *)
nick3=0; (* nick at the 3' end *)
xmin=0; xmax=0; ymin=0; ymax=0;  (* image size *)
l=0; r=0; d=0; cen=0; lr=1; (* length, radius, open angle, center, and loop ratio of a hairpin *)
i=1;
While[i<=Length[domains],
Which[
domains[[i]]=="+",
	lc=0;
	If[forward==1,
		kp={-1,-1}; i++,
		While[StringFreeQ[domains[[i]],")"], 
			If[StringTake[domains[[i]],1]=="@",
			ko+=ToExpression[StringTake[domains[[i]],{2,-1}]]]; ko=Mod[ko,360]; i++]; 
		forward=1], 
StringTake[domains[[i]],1]=="@",
	If[forward==1, ko+=ToExpression[StringTake[domains[[i]],{2,-1}]]; ko=Mod[ko,360]; i++,
		ko-=ToExpression[StringTake[domains[[i]],{2,-1}]]; ko=Mod[ko,360]; i--],
True,
	If[StringTake[domains[[i]],1]==":", nick5=1; domains[[i]]=StringTake[domains[[i]],{2,-1}]];
	If[StringTake[domains[[i]],-1]==":", nick3=1; domains[[i]]=StringTake[domains[[i]],{1,-2}]];
	Which[
	StringTake[domains[[i]],-1]=="(",
		kc++; lc=kd;
		label[[kd]]=StringTake[domains[[i]],{1,-2}];
		comp[[kc]]={kd,label[[kd]]};
		basepair[[kd]]=1,
	domains[[i]]==")",
		If[StringTake[comp[[kc,2]],-1]=="*",
			label[[kd]]=StringTake[comp[[kc,2]],{1,-2}],
			label[[kd]]=comp[[kc,2]]<>"*"]; 
		kp=pos[[comp[[kc,1]],2]]+RotationTransform[orient[[comp[[kc,1]]]] Degree][{0,pcomp}];
		ko=orient[[comp[[kc,1]]]]+180; ko=Mod[ko,360];
		If[lc>0,
			l=Sum[length[[j]],{j,lc+1,kd-1}];
			If[OptionValue[DNASequence],If[l<10,lr=10/l*1.6,lr=1.6],If[l<6,lr=6/l]];
			l=l*lr;
			r=Re[r/.FindRoot[r*(2Pi-2ArcSin[-pcomp/(2r)])==l,{r,1}]];
			d=ArcSin[-pcomp/(2r)];
			cen=(pos[[lc,2]]+kp)/2+RotationTransform[orient[[lc]] Degree][{Sqrt[r^2-(pcomp/2)^2],0}];
			If[cen-{r,0}<xmin,xmin=cen-{r,0}];
			If[cen+{r,0}>xmax,xmax=cen+{r,0}];
			If[cen-{0,r}<ymin,ymin=cen-{0,r}];
			If[cen+{0,r}>ymax,ymax=cen+{0,r}];
			For[j=lc+1,j<kd,j++,
				hairpin[[j]]=1;
				pos[[j,1]]=cen;
				If[j==lc+1,
					pos[[j,2]]={orient[[lc]]/360*2Pi+Pi-d,orient[[lc]]/360*2Pi+Pi-d-((2Pi-2d)*length[[j]]*lr/l)},
					pos[[j,2]]={pos[[j-1,2,2]],pos[[j-1,2,2]]-((2Pi-2d)*length[[j]]*lr/l)}];
				pos[[j,3]]=pos[[j,1]]+{Cos[(pos[[j,2,1]]+pos[[j,2,2]])/2]*(r+2),Sin[(pos[[j,2,1]]+pos[[j,2,2]])/2]*(r+2)};
				length[[j]]=r;
				If[pos[[j,3,1]]<xmin,xmin=pos[[j,3,1]]];
				If[pos[[j,3,1]]>xmax,xmax=pos[[j,3,1]]];
				If[pos[[j,3,2]]<ymin,ymin=pos[[j,3,2]]];
				If[pos[[j,3,2]]>ymax,ymax=pos[[j,3,2]]];
			];
			lc=0
		];
		kc--,
	True, 
		label[[kd]]=domains[[i]]];
	If[StringTake[label[[kd]],-1]=="*",
		name=StringTake[label[[kd]],{1,-2}]; 
		If[OptionValue[DistinctComplement],style[[kd]]=2,style[[kd]]=1]; 
		If[OptionValue[ReverseOrient], 
			seq[[kd]]=StringReverse[StringReplace[StringReverse[ToExpression[name]],{"C"->"G","G"->"C","T"->"A","A"->"T"}]],
			seq[[kd]]=StringReverse[StringReplace[ToExpression[name],{"C"->"G","G"->"C","T"->"A","A"->"T"}]]
		],
		If[StringTake[label[[kd]],1]=="_",
			name=StringTake[label[[kd]],{2,-1}], 
			name=label[[kd]]
		]; 
		style[[kd]]=1;
		If[OptionValue[ReverseOrient], 
			seq[[kd]]=StringReverse[ToExpression[name]],
			seq[[kd]]=ToExpression[name]
		]];
	length[[kd]]=StringLength[ToExpression[name]];
	color[[kd]]=cl[ToString[name]];
	orient[[kd]]=ko;
	If[StringTake[domains[[i]],1]!="_" && StringTake[domains[[i]],1]!="@" && domains[[i]]!="+" &&
		((!OptionValue[ReverseOrient] && (i==Length[domains] || domains[[i+1]]=="+" || StringTake[domains[[i+1]],1]=="_")) ||
		(OptionValue[ReverseOrient] && (i==1 || domains[[i-1]]=="+" || StringTake[domains[[i-1]],1]=="_" || 
										((i==2 || domains[[i-2]]=="+") && StringTake[domains[[i-1]],1]=="@")))), 
		style[[kd]]+=2];
	If[kp=={-1,-1}, 
		While[StringFreeQ[domains[[i]],")"], i++]; 
		kp=pos[[comp[[kc,1]],2]]+RotationTransform[orient[[comp[[kc,1]]]] Degree][{0,pcomp}];
		ko=orient[[comp[[kc,1]]]]+180; ko=Mod[ko,360];
		i--;
		forward=0, 
			
		If[StringTake[domains[[i]],1]=="_", 
			If[i==Length[domains] || domains[[i+1]]=="+",
				pos[[kd,1]]=kp+RotationTransform[orient[[kd]] Degree][{0.7,0}]; 
				pos[[kd,2]]=pos[[kd,1]];
				pos[[kd,3]]=pos[[kd,1]]+RotationTransform[orient[[kd]] Degree][{1.5+StringLength[seq[[kd]]]/2,0}],

				pos[[kd,1]]=kp-RotationTransform[orient[[kd]] Degree][{0.7,0}]; 
				pos[[kd,2]]=pos[[kd,1]];
				pos[[kd,3]]=pos[[kd,1]]-RotationTransform[orient[[kd]] Degree][{1.5+StringLength[seq[[kd]]]/2,0}]];
			If[forward==1, i++, i--],

		If[forward==1,
			pos[[kd,1]]=kp;
			pos[[kd,2]]=kp+RotationTransform[orient[[kd]] Degree][{length[[kd]],0}];
			If[EvenQ[orient[[kd]]/90],
				pos[[kd,3]]=(pos[[kd,1]]+pos[[kd,2]])/2+RotationTransform[orient[[kd]] Degree][{0,2}],
				pos[[kd,3]]=(pos[[kd,1]]+pos[[kd,2]])/2+RotationTransform[orient[[kd]] Degree][{0,1.5+StringLength[label[[kd]]]/2}]
			];				
			If[OptionValue[DNASequence], pos[[kd,4]]=(pos[[kd,1]]+pos[[kd,2]])/2+RotationTransform[orient[[kd]] Degree][{0,-2}]];
			pos[[kd,1]]=kp+RotationTransform[orient[[kd]] Degree][{nick5,0}];
			pos[[kd,2]]=kp+RotationTransform[orient[[kd]] Degree][{length[[kd]]-nick3,0}];
			kp=pos[[kd,2]];
			If[nick5==1, domains[[i]]=":"<>domains[[i]]; nick5=0];
			If[nick3==1, domains[[i]]=domains[[i]]<>":"; nick3=0];
			i++,
	
			pos[[kd,2]]=kp;
			pos[[kd,1]]=kp-RotationTransform[orient[[kd]] Degree][{length[[kd]],0}];	
			If[EvenQ[orient[[kd]]/90],
				pos[[kd,3]]=(pos[[kd,1]]+pos[[kd,2]])/2+RotationTransform[orient[[kd]] Degree][{0,2}],
				pos[[kd,3]]=(pos[[kd,1]]+pos[[kd,2]])/2+RotationTransform[orient[[kd]] Degree][{0,1.5+StringLength[label[[kd]]]/2}]
			];			
			If[OptionValue[DNASequence], pos[[kd,4]]=(pos[[kd,1]]+pos[[kd,2]])/2+RotationTransform[orient[[kd]] Degree][{0,-2}]];
			pos[[kd,2]]=kp-RotationTransform[orient[[kd]] Degree][{nick3,0}];
			pos[[kd,1]]=kp-RotationTransform[orient[[kd]] Degree][{length[[kd]]-nick5,0}];
			kp=pos[[kd,1]];
			If[nick5==1, domains[[i]]=":"<>domains[[i]]; nick5=0];
			If[nick3==1, domains[[i]]=domains[[i]]<>":"; nick3=0];
			i--
		]];
		For[j=1,j<=4,j++,
			If[pos[[kd,j,1]]<xmin, xmin=pos[[kd,j,1]]];
			If[pos[[kd,j,1]]>xmax, xmax=pos[[kd,j,1]]];
			If[pos[[kd,j,2]]<ymin, ymin=pos[[kd,j,2]]];
			If[pos[[kd,j,2]]>ymax, ymax=pos[[kd,j,2]]]];
		kd++]
	]
];
Style[Graphics[{
Table[If[OptionValue[GreyBasepair] && basepair[[i]]==1,
		{Lighter[Gray,0.5],Polygon[{pos[[i,1]],pos[[i,2]],
									pos[[i,2]]+RotationTransform[orient[[i]] Degree][{0,pcomp}],
									pos[[i,1]]+RotationTransform[orient[[i]] Degree][{0,pcomp}]}]}],
{i,1,Nd}],
Table[{
	If[OptionValue[LineBasepair] && basepair[[i]]==1,
		Table[{Black,Line[{pos[[i,1]]+RotationTransform[orient[[i]] Degree][{j-0.5,0}]+RotationTransform[orient[[i]] Degree][{0,pcomp*0.2}],
				pos[[i,1]]+RotationTransform[orient[[i]] Degree][{j-0.5,0}]+RotationTransform[orient[[i]] Degree][{0,pcomp*0.8}]}]},
		{j,1,length[[i]]}]],
	If[OptionValue[DotBasepair] && basepair[[i]]==1,
		Table[{Black,Disk[pos[[i,1]]+RotationTransform[orient[[i]] Degree][{j-0.5,0}]+RotationTransform[orient[[i]] Degree][{0,pcomp*0.5}],0.3]},
		{j,1,length[[i]]}]],
	If[OptionValue[DomainSeparator] && style[[i]]<3 && StringTake[label[[i]],1]!="_",
		If[hairpin[[i]]==1,
			{Black,Line[{pos[[i,1]]+{Cos[pos[[i,2,2]]]*length[[i]],Sin[pos[[i,2,2]]]*length[[i]]},
						pos[[i,1]]+{Cos[pos[[i,2,2]]]*(length[[i]]+pcomp/2),Sin[pos[[i,2,2]]]*(length[[i]]+pcomp/2)}}]},
			{Black,Line[{pos[[i,2]],pos[[i,2]]+RotationTransform[orient[[i]] Degree][{0,pcomp*0.5}]}]}]],

	If[StringTake[label[[i]],1]=="_",
		DomainObject=Disk[pos[[i,1]],0.7];
		TextObject=Text[Style[seq[[i]],Black,FontSize->20,FontFamily->"Calibri"],pos[[i,3]]],
		If[OptionValue[SquiggleToehold] && StringLength[seq[[i]]]<=OptionValue[MaxToeholdLength], 
			If[hairpin[[i]]==1,
				If[!OptionValue[ReverseOrient],
				DomainObject=Line[Flatten[{{pos[[i,1]]+{Cos[pos[[i,2,1]]]*r,Sin[pos[[i,2,1]]]*r}}, 
							Table[pos[[i,1]]+{Cos[pos[[i,2,1]]-(j-0.5)*(pos[[i,2,1]]-pos[[i,2,2]])/(lr*StringLength[seq[[i]]])]*(r+0.4*(2*Mod[2j,2]-1)),
								Sin[pos[[i,2,1]]-(j-0.5)*(pos[[i,2,1]]-pos[[i,2,2]])/(lr*StringLength[seq[[i]]])]*(r+0.4*(2*Mod[2j,2]-1))},
								{j,1,lr*StringLength[seq[[i]]]-1,0.5}],
							{pos[[i,1]]+{Cos[pos[[i,2,2]]+(pos[[i,2,1]]-pos[[i,2,2]])/(lr*StringLength[seq[[i]]])]*r,
									Sin[pos[[i,2,2]]+(pos[[i,2,1]]-pos[[i,2,2]])/(lr*StringLength[seq[[i]]])]*r},
							pos[[i,1]]+{Cos[pos[[i,2,2]]]*r,Sin[pos[[i,2,2]]]*r}}},1]],
				DomainObject=Line[Flatten[{{pos[[i,1]]+{Cos[pos[[i,2,1]]]*r,Sin[pos[[i,2,1]]]*r}, 
							pos[[i,1]]+{Cos[pos[[i,2,1]]-(pos[[i,2,1]]-pos[[i,2,2]])/(lr*StringLength[seq[[i]]])]*r,
									Sin[pos[[i,2,1]]-(pos[[i,2,1]]-pos[[i,2,2]])/(lr*StringLength[seq[[i]]])]*r}},
							Table[pos[[i,1]]+{Cos[pos[[i,2,1]]-(j-0.5)*(pos[[i,2,1]]-pos[[i,2,2]])/(lr*StringLength[seq[[i]]])]*(r+0.4*(2*Mod[2j,2]-1)),
								Sin[pos[[i,2,1]]-(j-0.5)*(pos[[i,2,1]]-pos[[i,2,2]])/(lr*StringLength[seq[[i]]])]*(r+0.4*(2*Mod[2j,2]-1))},
								{j,2,lr*StringLength[seq[[i]]],0.5}],
							{pos[[i,1]]+{Cos[pos[[i,2,2]]]*r,Sin[pos[[i,2,2]]]*r}}},1]]],
				If[!OptionValue[ReverseOrient],
				DomainObject=Line[Flatten[{{pos[[i,1]]}, Table[pos[[i,1]]+RotationTransform[orient[[i]] Degree][{j-0.5,0}]+
										RotationTransform[orient[[i]] Degree][{0,0.4*(2*Mod[2j,2]-1)}],{j,1,length[[i]]-1,0.5}],
								{pos[[i,2]]+RotationTransform[orient[[i]] Degree][{-1,0}],pos[[i,2]]}},1]],
				DomainObject=Line[Flatten[{{pos[[i,1]],pos[[i,1]]+RotationTransform[orient[[i]] Degree][{1,0}]},
										Table[pos[[i,1]]+RotationTransform[orient[[i]] Degree][{j-0.5,0}]+
										RotationTransform[orient[[i]] Degree][{0,0.4*(2*Mod[2j,2]-1)}],{j,2,length[[i]],0.5}],{pos[[i,2]]}},1]]
				]
			],
			If[hairpin[[i]]==1,		
			DomainObject=Circle[pos[[i,1]],length[[i]],pos[[i,2]]],	
			DomainObject=Line[{pos[[i,1]],pos[[i,2]]}]]
		];
		TextObject=Text[Style[label[[i]],Black,FontSize->OptionValue[DomainLabelSize],FontFamily->"Calibri"],pos[[i,3]]]];

	If[IntegerQ[color[[i]]],plotColor=domainColor[[color[[i]]]],plotColor=color[[i]]];

	If[OddQ[style[[i]]],
	{plotColor,AbsoluteThickness[OptionValue[LineThickness]],DomainObject,If[OptionValue[DomainLabel],TextObject,{}]},
	{{Lighter[plotColor,0.5],AbsoluteThickness[OptionValue[LineThickness]],DomainObject,If[OptionValue[DomainLabel],TextObject,{}]},
		If[OptionValue[VectorFormat],
	{plotColor,AbsoluteThickness[OptionValue[LineThickness]],AbsoluteDashing[{1,7}],DomainObject},
	{plotColor,AbsoluteThickness[OptionValue[LineThickness]],Dashed,DomainObject}]}],

	If[style[[i]]>2,
		If[!OptionValue[ReverseOrient],
			a1=pos[[i,2]]+RotationTransform[orient[[i]] Degree][{-1,0.5}];
			a2=pos[[i,2]]+RotationTransform[orient[[i]] Degree][{-1,-0.5}],
			a1=pos[[i,1]]+RotationTransform[orient[[i]] Degree][{1,0.5}];
			a2=pos[[i,1]]+RotationTransform[orient[[i]] Degree][{1,-0.5}]
		];
		If[a1[[1]]<xmin, xmin=a1[[1]]];
		If[a1[[1]]>xmax, xmax=a1[[1]]];
		If[a1[[2]]<ymin, ymin=a1[[2]]];
		If[a1[[2]]>ymax, ymax=a1[[2]]];
		If[a2[[1]]<xmin, xmin=a2[[1]]];
		If[a2[[1]]>xmax, xmax=a2[[1]]];
		If[a2[[2]]<ymin, ymin=a2[[2]]];
		If[a2[[2]]>ymax, ymax=a2[[2]]];
		If[!OptionValue[ReverseOrient],	
			{{White,AbsoluteThickness[OptionValue[LineThickness]],Line[{pos[[i,2]]+RotationTransform[orient[[i]] Degree][{-0.1,0}],pos[[i,2]]}]},
			{plotColor,AbsoluteThickness[OptionValue[LineThickness]], Line[{a1,pos[[i,2]],a2}]}},
			{{White,AbsoluteThickness[OptionValue[LineThickness]],Line[{pos[[i,1]]+RotationTransform[orient[[i]] Degree][{-0.1,0}],pos[[i,1]]}]},
			{plotColor,AbsoluteThickness[OptionValue[LineThickness]], Line[{a1,pos[[i,1]],a2}]}}
		]
	],

	If[OptionValue[DNASequence] && StringTake[label[[i]],1]!="_", 
		If[hairpin[[i]]==1,
		Table[Text[Style[StringTake[seq[[i]],{j}],plotColor,FontSize->15,FontFamily->"Courier New"],
			pos[[i,1]]+{Cos[pos[[i,2,1]]-(j-0.5)*(pos[[i,2,1]]-pos[[i,2,2]])/StringLength[seq[[i]]]]*(r-2),
						Sin[pos[[i,2,1]]-(j-0.5)*(pos[[i,2,1]]-pos[[i,2,2]])/StringLength[seq[[i]]]]*(r-2)},
			Automatic,RotationTransform[pos[[i,2,1]]-(j-0.5)*(pos[[i,2,1]]-pos[[i,2,2]])/StringLength[seq[[i]]]-Pi/2][{1,0}]],
		{j,1,StringLength[seq[[i]]]}],
		Text[Style[seq[[i]],plotColor,FontSize->15,FontFamily->"Courier New"],
			pos[[i,4]],Automatic,RotationTransform[orient[[i]] Degree][{1,0}]]]]
	},{i,1,Nd}]},
ImageSize->10*{xmax-xmin+10,ymax-ymin+10},PlotRangePadding->None,ImagePadding->50,ImageMargins->0],
ImageSizeMultipliers->1, ScriptSizeMultipliers->1]
]


AddTurns[molecule_]:=
Module[{domains,Nd,comp,c1,c2,kc,ns,nl,nt,complex,i},
domains=StringSplit[molecule," "];
Nd=Sum[If[domains[[i]]=="+",0,1],{i,1,Length[domains]}]; (* number of domains *)
comp=Table[0,{i,1,Nd}]; (* a stack of complemetary domains *)
c1=0; (* last unpaired ( in the previous strand *)
c2=0; (* last unpaired ( in the current strand *)
kc=0; (* current complementary domain that hasn't been paired *)
ns=1; (* current strand number *)
nl=0; (* number of ( in the current strand *)
nt=0; (* total number of turns *)
For[i=1,i<=Length[domains],i++,
Which[
domains[[i]]=="+", 
	If[c2>c1 && ns>1, domains=Insert[domains,"@",comp[[c2]]+1]; nt++; i++];
	ns++; nl=0; c1=c2,
StringTake[domains[[i]],-1]=="(",
	kc++;
	comp[[kc]]=i;
	If[nl==0 && ns>1, c2=kc-1; domains=Insert[domains,"@",comp[[kc]]]; nt++; i++];
	nl++,
domains[[i]]==")",
	kc--]
];
complex=domains[[1]];
For[i=2,i<=Length[domains],i++,
If[domains[[i]]=="@", domains[[i]]="@"<>ToString[If[nt==1,90,180-360/(nt+1)]]];
complex=complex<>" "<>domains[[i]]
];
complex]


AutoShift[molecule_]:=
Module[{domains,newdomains,Nd,comp,label,pair,nr,kc,nor,complex,i,j},
domains=StringSplit[molecule,"+"];
For[i=1,i<=Length[domains],i++,domains[[i]]=StringSplit[domains[[i]]," "]];
comp=Table[{0,0,""},{i,1,Sum[Length[domains[[i]]],{i,1,Length[domains]}]}]; (* a stack of complemetary domains *)
label=Table["",{i,1,Length[domains]},{j,1,Length[domains[[i]]]}];
pair=Table[{0,0},{i,1,Length[domains]},{j,1,Length[domains[[i]]]}];
nr=Table[0,{i,1,Length[domains]}]; (* number of ) in each strand *)
kc=0; (* current complementary domain that hasn't been paired *)
nor=0; (* the last strand without any ) *)
For[i=1,i<=Length[domains],i++,
	For[j=1,j<=Length[domains[[i]]],j++,
		Which[
		StringTake[domains[[i,j]],-1]=="(",
			kc++;
			label[[i,j]]=StringTake[domains[[i,j]],{1,-2}];
			comp[[kc]]={i,j,label[[i,j]]},
		domains[[i,j]]==")",
			nr[[i]]++;
			pair[[i,j]]={comp[[kc,1]],comp[[kc,2]]};
			pair[[comp[[kc,1]],comp[[kc,2]]]]={i,j};
			If[StringTake[comp[[kc,3]],-1]=="*",
				label[[i,j]]=StringTake[comp[[kc,3]],{1,-2}],
				label[[i,j]]=comp[[kc,3]]<>"*"]; 
			kc--,
		True, 
			label[[i,j]]=domains[[i,j]]]
	];
	If[nr[[i]]==0,nor=i]
];
newdomains=domains;
If[nor>1,
domains=RotateLeft[domains,nor-1];
label=RotateLeft[label,nor-1];
pair=RotateLeft[pair,nor-1];
newdomains=domains;
For[i=1,i<=Length[pair],i++,
	For[j=1,j<=Length[pair[[i]]],j++,
	If[pair[[i,j]]!={0,0},	
	pair[[i,j]]={If[pair[[i,j,1]]>=nor,pair[[i,j,1]]-nor+1,Length[pair]-nor+1+pair[[i,j,1]]],pair[[i,j,2]]}	
	]]
];
kc=0;
For[i=1,i<=Length[domains],i++,
	For[j=1,j<=Length[domains[[i]]],j++,
		Which[
		StringTake[domains[[i,j]],-1]=="(",
			kc++,
		domains[[i,j]]==")",
			kc--;
			If[kc<0,
				newdomains[[i,j]]=label[[i,j]]<>"(";
				newdomains[[pair[[i,j,1]],pair[[i,j,2]]]]=")"]
		]
]]];
complex="";
For[i=1,i<=Length[newdomains],i++,
	For[j=1,j<=Length[newdomains[[i]]],j++,
		If[j>1,complex=complex<>" "<>newdomains[[i,j]],
				complex=complex<>newdomains[[i,j]]]];
	If[i<Length[newdomains],complex=complex<>" + "]];
complex]
