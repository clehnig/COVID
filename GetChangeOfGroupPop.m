function [dg1,dg2,dg3,dg0] = GetChangeOfGroupPop(group1,group2,group3,group0,I,N)


% This function returns the compartment size changes over a single day. Used in a day loop to update compartment counts

% SEIAR Model SDSU Student Popualation Framework
%-----------------------------------
%       SDSU Compartments 
% ----------------------------------
% -- (on campus) --
% Sc = Susceptible
% Ec = Exposed
% Ic = Symptomatic Infectious
% Ac = Asymptomatic Infectious
% Rc = Recovered
% Nc = total population (Nc=Sc+Ec+Ic+Ac+Rc)

% -- (off campus) --
% So = Susceptible
% Eo = Exposed
% Io = Symptomatic Infectious
% Ao = Asymptomatic Infectious
% Ro = Recovered
% No = total population (No=So+Eo+Io+Ao+Ro)

% -- San Diego County --
% I = Infected 
% N = Population

% -- Parameters (on campus) --
% bc_ic = prob(transmission from Ic to Sc) 
% bc_ac = prob(transmission from Ac to Sc)
% bc_io = prob(transmission from Io to Sc)
% bc_ao = prob(transmission from Ao to Sc)

% -- Parameters (off campus) --
% bo_io = prob(transmission from Io to So) 
% bo_ao = prob(transmission from Ao to So)
% bo_ic = prob(transmission from Ic to So)
% bo_ac = prob(transmission from Ac to So)

% s  = 1/(incubation period) in days 
% u = 1/(infectious period) 
% p  = prob(symptomatic)


bc_ic = .14;
bc_ac = bc_ic;
bc_io = bc_ic*.4;
bc_ao = bc_ic*.4;

bo_io = bc_ic*0.8;
bo_ao = bo_io;
bo_ic = bo_io*.4;
bo_ac = bo_io*.4;

bc = .2;
bo = .2;

s  = 1/5;   % Incubation Period
u = 1/14;   % Infectious Period

p = .9;

Ic=group1(3)+group2(3);
Ac=group1(4)+group2(4);

Io=group3(3)+group0(3);
Ao=group3(4)+group0(4);

Nc = floor(sum(group1) + sum(group2)); 
No = floor(sum(group3) + sum(group0));

    dg1(1) = -group1(1) * ( (bc_ic*Ic + bc_ac*Ac)/Nc + (bc_io*Io + bc_ao*Ao)/No + bc*I/N);
    dg1(2) = abs(dg1(1)) - s*group1(2);
    dg1(3) =   (p)*s*group1(2) - u*group1(3);
    dg1(4) = (1-p)*s*group1(2) - u*group1(4);
    dg1(5) = u * (group1(3) + group1(4));
    
    dg2(1) = -group2(1) * ( (bc_ic*Ic + bc_ac*Ac)/Nc + (bc_io*Io + bc_ao*Ao)/No + bc*I/N);
    dg2(2) = abs(dg2(1)) - s*group2(2);
    dg2(3) =   (p)*s*group2(2) - u*group2(3);
    dg2(4) = (1-p)*s*group2(2) - u*group2(4);
    dg2(5) = u * (group2(3) + group2(4));
    
    
    dg3(1) = -group3(1) * ( (bo_io*Io + bo_ao*Ao)/No + (bo_ic*Ic + bo_ac*Ac)/Nc + bo*I/N);
    dg3(2)  = abs(dg3(1)) - s*group3(2);
    dg3(3)  =   (p)*s*group3(2) - u*group3(3);
    dg3(4)  = (1-p)*s*group3(2) - u*group3(4);
    dg3(5)  = u * (group3(3) + group3(4));
    
    dg0(1) = -group0(1) * ( (bo_io*Io + bo_ao*Ao)/No + (bo_ic*Ic + bo_ac*Ac)/Nc + bo*I/N);
    dg0(2)  = abs(dg0(1)) - s*group0(2);
    dg0(3)  =   (p)*s*group0(2) - u*group0(3);
    dg0(4)  = (1-p)*s*group0(2) - u*group0(4);
    dg0(5)  = u * (group0(3) + group0(4));
    
    
    
end
