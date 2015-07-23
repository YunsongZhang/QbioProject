%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------Moose code-----------
% ----Yunsong Zhang-------
%-----2015-7-22-----------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code tries to implement a simplest delayed degradation and fire model
% in the 2009 PRL paper by Will Mather et al
% Two reactions are considered:
%   0 --> r  Ka
%   r --> 0  Kd

clear
close all
clc

V0=1;                    % Initial Volume
T_divide=10;             % period of cell division
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

timenow=0;
numR=0;
c=log(2)/T_divide;
%%
while timenow<Tmax
    V=Volume(timenow);
    propensity=Model.K(numR/V);
    K_sum=sum(propensity);
    delta_t=log(K_sum/(K_sum+c*log(rand())))/c;
%     T_nextdelay(find(T_nextdelay<timenow))=[];
%     if numR==0
%         timenow=timenow+Model.tau;
%         numR=numR+1;
%         T_rec(end+1)=timenow;
%         numR_rec(end+1)=numR;
%         disp('cont');
%         %timenow
%         %T_nextdelay
%         continue;        
%     end
         
    if timenow+delta_t>1000000 %T_divide
        test=1;
        % cell divides
        %disp('Cell divides');
    elseif timenow+delta_t>min(T_nextdelay) %%&& timenow<min(T_nextdelay)
       test=2; min(T_nextdelay)
        % a delay reaction occurs
        timenow=min(T_nextdelay);
        T_rec(end+1)=timenow;
        numR=numR+1;
        numR_rec(end+1)=numR;
        ind=find(T_nextdelay==min(T_nextdelay));
        T_nextdelay(ind)=[];
        
        
        
        
    elseif rand()<propensity(1)/K_sum
        test=3;
        % a delayed reaction in preparation
        timenow=timenow+delta_t;
        T_nextdelay(end+1)=timenow+Model.tau;
        T_nextdelay=unique(T_nextdelay);
        
        
        
        
    else
        test=4;
        % a reaction occurs immediately       
        timenow=timenow+delta_t;
        numR=numR-1;
        T_rec(end+1)=timenow;
        numR_rec(end+1)=numR;
    end
    
    if delta_t<0
        error('strange delta_t');
    end
%     if T_rec(end)<T_rec(end-1)
%         error('time back!');
%     end
       
        
  
   timenow
        

end

figure(1)
handle=plot(T_rec,numR_rec);
xlim([min(T_rec),max(T_rec)])
titlename=sprintf('number of r, T_division=%d',T_divide);
title(titlename)
picname=sprintf('./Td%dnum',T_divide);
saveas(handle,picname,'jpg');



figure(2)
handle=plot(T_rec,numR_rec./Volume(T_rec));
xlim([min(T_rec),max(T_rec)])
titlename=sprintf('concentration of r, T_division=%d',T_divide);
title(titlename)
picname=sprintf('./Td%dcon',T_divide);
saveas(handle,picname,'jpg');
