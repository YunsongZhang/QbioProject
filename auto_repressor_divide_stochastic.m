%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------Moose code-----------
% ----Yunsong Zhang-------
%-----2015-7-23-----------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code tries to implement a simplest delayed degradation and fire model
% in the 2009 PRL paper by Will Mather et al
% Two reactions are considered:
%   0 --> r  Ka
%   r --> 0  Kd
% However, here we also consider the dilution and random partitioning of
% cell growth and division
clear
close all
clc

V0=1;                    % Initial Volume
T_divide=20;             % period of cell division
Volume=@(t) V0*exp(log(2).*t/T_divide);

% chemical model definition
Model.gamma_r=80;
Model.alpha=300;
Model.C0=10;
Model.tau=1;
Model.beta=0.1;
Model.R0=1;
Model.S=[1,-1];
Model.K=@(x) [Model.alpha*(Model.C0/(Model.C0+x))^2;
               Model.gamma_r*x/(Model.R0+x)+Model.beta*x];
%end of model definition

T_rec=[0];
numR_rec=[0];  
Tmax=200; %2*T_divide;
T_nextdelay=[NaN];

timenow=0;           % real time
time_in_gen=0;       % time in each generation, put to 0 at each division
numR=0;
generation=1;       % number of generation
c=log(2)/T_divide;
cells={};
%%

      rng('shuffle');
      
       [T_rec,numR_rec,con_rec]=stochastic_reaction(numR,Model,Volume,T_divide);
       cells{1}={[T_rec;numR_rec;con_rec]};
       
       numTotal=numR_rec(end);
       
       numR1=floor(numTotal*0.5);
       numR2=numTotal-numR1;
       [T_rec1,numR_rec1,con_rec1]=stochastic_reaction(numR1,Model,Volume,T_divide*(1+0.2*rand()));
       [T_rec2,numR_rec2,con_rec2]=stochastic_reaction(numR2,Model,Volume,T_divide*(1+0.2*rand()));
       cells{2}={[T_rec1;numR_rec1;con_rec1],[T_rec2;numR_rec2;con_rec2]};
      % timenow=timenow+T_divide;
      % generation=generation+1;
   
     figure(1)
     plot(T_rec,numR_rec)
     title('generation 1, number of R');
     shg
     
     figure(2)
     plot(T_rec,con_rec)
     title('generation 1, concentration of R');
     shg
     
     figure(3)
     plot(T_rec1,numR_rec1)
     title('generation 2, number of R,cell1');
     shg
     
     figure(4)
     plot(T_rec1,con_rec1)
     title('generation 2, concentration of R,cell1');
     shg
     
     figure(5)
     plot(T_rec2,numR_rec2)
     title('generation 2, number of R,cell2');
     shg    
     
      figure(6)
     plot(T_rec2,con_rec2)
     title('generation 2, concentration of R,cell2');
     shg       
     
     figure(7)
     plot(T_rec1,con_rec1)
     hold on
     plot(T_rec2,con_rec2)
     legend('cell1','cell2')
     shg
        
   



% figure(1)
% handle=plot(T_rec,numR_rec);
% xlim([min(T_rec),max(T_rec)])
% titlename=sprintf('number of r, T_division=%d',T_divide);
% title(titlename)
% picname=sprintf('./Td%dnum',T_divide);
%saveas(handle,picname,'jpg');



% figure(2)
% handle=plot(T_rec,numR_rec./Volume(T_rec));
% xlim([min(T_rec),max(T_rec)])
% titlename=sprintf('concentration of r, T_division=%d',T_divide);
% title(titlename)
% picname=sprintf('./Td%dcon',T_divide);
%saveas(handle,picname,'jpg');
