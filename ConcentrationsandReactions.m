function [concentrations, reactions] = ConcentrationsandReactions()
%test test
concentrations{1,1} = 'P680';
concentrations{1,2} = 'P680*';
concentrations{1,3} = 'Phe';
concentrations{1,4} = 'P680+'; 
concentrations{1,5} = 'e-'; 
concentrations{1,6} = 'Phe-'; 
concentrations{1,7} = 'QA';
concentrations{1,8} = 'QA-';
concentrations{1,9} = 'QB'; 
concentrations{1,10} = 'QB-'; 
concentrations{1,11} = 'H+,Stroma'; 
concentrations{1,12} = 'QBH-'; 
concentrations{1,13} = 'QBH2'; 
concentrations{1,14} = 'C02 Ext';
concentrations{1,15} = 'p700*';
concentrations{1,16} = 'A0';
concentrations{1,17} = 'A0-';
concentrations{1,18} = 'A1';
concentrations{1,19} = 'Pc(Cu+)';
concentrations{1,20} = 'P700+';
concentrations{1,21} = 'Pc(Cu2+)';
concentrations{1,22} = 'A1-';
concentrations{1,23} = 'Fd';
concentrations{1,24} = 'Fd-';
concentrations{1,25} = 'NADP+';
concentrations{1,26} = 'NADPH';
concentrations{1,27} = 'ADP';
concentrations{1,28} = 'ATP';
concentrations{1,29} = 'CO2';
concentrations{1,30} = 'H+,Lumen';
concentrations{1,31} = 'PQH2';
concentrations{1,32} = 'PQ';
concentrations{1,33} = 'P700';
concentrations{1,34} = 'P680d';
concentrations{1,35} = 'P680d*';

reactions{1} = 'P680 -> p680*';
reactions{2} = 'P680* + Phe -> P680+ + Phe-'; 
reactions{3} = 'P680+ + e- -> Phe + QA-';
reactions{4} = 'Phe- + QA -> Phe + QA-';
reactions{5} = 'p680+ + QA- -> P680 + QA';
reactions{6} = 'QA- + QB -> QA + QB-';
reactions{7} = 'QA- + QB -> QA + QBH-';
reactions{8} = 'QBH- + H+,Stroma -> QBH2';
reactions{9} = 'QBH2 -> PQH2';
reactions{10} = 'PQH2 + 2PC(CU2++ -> PQ + 2H+,Lumen';
reactions{11} = 'P700 -> p700*';
reactions{12} = 'P700* + A0 -. P700+ + A0-';
reactions{13} = 'A0 + A1 -> A0 + A1-';
reactions{14} = 'PC(Cu+) + P700+ -> A0 + A1-';
reactions{15} = 'P700+ + A1- -> A1 + Fd-';
reactions{16} = 'A1- + Fd -> A1 + Fd-';
reactions{17} = '2Fd- + NADP+ + H+,Stroma -> 2Fd + NADPH';
reactions{18} = 'P700* -> P700';
reactions{19} = 'P680* -> P680';
reactions{20} = 'NADPH -> NADP+ + H+,Stroma';
reactions{21} = 'Fd- + H+,Stroma -> Fd';
reactions{22} = 'ADP + H+,Lumen -> ATP + H+,Stroma';
reactions{23} = 'CO2 + 3ATP+ 2NADPH -> 3ADP + 2NADP+ + H+,Stroma';
reactions{24} = 'PQ -> QB';
reactions{25} = 'ATP -> ADP';
reactions{26} = 'H+,Lumen -> H+, Stroma';
reactions{27} = 'CO2 Ext -> CO2';
reactions{28} = 'PQH@ -> PQ';
reactions{29} = 'NADP+ + H+,Stroma -> NADPH'; 
reactions{30} = 'Fd -> H+,Stroma + Fd-';
reactions{31} = '2H20 -> 4H+,Lumen + 4e- + O2';
reactions{32} = 'P680d -> P680d*';
reactions{33} = 'P680d* -> P680d'; 
reactions{34} = 'H+,Lumen ->';
