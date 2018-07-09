% function find_optimal(analysis_name,nloop)
% calculates optimal calies for kinetic constants and concentrations based
% on experimental data of ANALYSIS


C = xlsread('Coh.csv');
FvFm_exp = C(:,2);

options = optimset('Display', 'iter', 'GradObj', 'off', 'Algorithm', 'interior-point');


k = rand(32,1);
y0 = rand(35,1);
rng('shuffle');
x0 = [y0;k];
%initialize constraint matrices Aeq * x = beq. Each row of Aeq corresponds
%to the coefficients for a linear equation. beq is the right hand side of
%that equation.
Aeq = zeros(10,length(x0));
beq = zeros(10,1);
%set linear constraints for p680 (sum(P680 = 1)
Aeq([1 2 4]) = ones(3,1);
beq(1) = 1;

%set linear equations for Phe
Aeq([3,6]) = ones(2,1);
beq(2) = 1;

%set linear equations for Qa
Aeq([7,8]) = ones(2,1);
beq(3) = 1;

%set linear equations for Qb
Aeq([9 10 12 13]) = ones(4,1);
beq(4) = 1;

xopt = fmincon(@(x) calc_sqerror(x,FvFm_exp),x0,[],[],Aeq,beq,zeros(length(x0),1),[],[],options);

foo = 1


