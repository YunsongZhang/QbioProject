function [T_rec,numR_rec,con_rec]=stochastic_reaction(InitialNum,Model,Volume,T_divide)

time_in_gen=0;
numR=InitialNum;
c=log(2)/T_divide;
T_nextdelay=[NaN];


T_rec=[0];
numR_rec=[InitialNum];

while time_in_gen<T_divide
    
propensity=Model.K(numR/Volume(time_in_gen));
K_sum=sum(propensity);
delta_t=log(K_sum/(K_sum+c*log(rand())))/c;

        if time_in_gen+delta_t>min(T_nextdelay) 
                % a delay reaction occurs
                time_in_gen=min(T_nextdelay);
               % T_rec(end+1)=time_in_gen;
                numR=numR+1;
               % numR_rec(end+1)=numR;
                ind=find(T_nextdelay==min(T_nextdelay));
                T_nextdelay(ind)=[];

        elseif rand()<propensity(1)/K_sum
                % a delayed reaction in preparation
                time_in_gen=time_in_gen+delta_t;
                T_nextdelay(end+1)=time_in_gen+Model.tau;
                T_nextdelay=unique(T_nextdelay);

        else
                % a reaction occurs immediately       
                time_in_gen=time_in_gen+delta_t;
                numR=numR-1;

        end

        T_rec(end+1)=time_in_gen;
        numR_rec(end+1)=numR;
        
        
end

  
   con_rec=numR_rec./Volume(T_rec);
      
end
