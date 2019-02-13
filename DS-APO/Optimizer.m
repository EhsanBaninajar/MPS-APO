clc
clear all
%% Loading global variables
Global_var
%% objective Function
nsim=1;    % Number of function evaluation
gap_num=15;    % number of gaps
gapsize=15;    %gap size
nb_realz=1;    % nb_realz: number of realizations in each function evaluation
Ttradeoff=10;    %user defined threshold between quality and cpu time (percent)
%% Loading Ti
%Training Image
[Ti,Xsize_Ti,Ysize_Ti,Zsize_Ti,nbvar_Ti,namevar_Ti]=LoadGrid('MultiVariate.sgems');
TiType=[1 1 0]; % "0"for categorical Ti and "1" for continious Ti
gapdatabase
%% SPSA
%parameter initialization
theta=[0.3; 0.1; 0.01; 0.33; 0.33; 0.34 ];% Maximum number of conditioning data, minimum number of replicates, number of multiple grids
thetamin=[0.001; 0.1; 0.01; 0.1; 0.1; 0.1];
thetamax=[0.2; 1; 1; 1; 1; 1;];
p=6;  % Dimension of the parameters search space
%% Optimizer parameters
NSPSA=50;   % Maximum number of iterations
a=2;
A=50;
c=0.2;
alpha=0.602;
gamma=0.101;
stepsizemax=0.01; %to be sure the second step don't degrad the quality a lot
stage=0;
counter=0;
seccounter=1;
FunValue=zeros(NSPSA,p+2);
FunValue_minus=zeros(NSPSA,p+2);
FunValue_plus=zeros(NSPSA,p+2);
FunValue_MSE=zeros(nsim,NSPSA,1);
stepsizestor=zeros(NSPSA,p);
%%
for k=0:NSPSA-1
    ak=a/(k+1+A)^alpha;
    ck=c/(k+1)^gamma;
    delta=2*round(rand(p,1))-1;
    thetaplus=theta+ck*delta;
    thetaminus=theta-ck*delta;
    
    %parameter constraints
    thetaplus=min(thetaplus,thetamax);
    thetaplus=max(thetaplus,thetamin);
    thetaminus=min(thetaminus,thetamax);
    thetaminus=max(thetaminus,thetamin); 
    thetaplus(4:6)=thetaplus(4:6)/sum(thetaplus(4:6));
    thetaminus(4:6)=thetaminus(4:6)/sum(thetaminus(4:6));
    
    %function evaluations
    Funnum=1;
    yvalue=E_fun(theta);
    Funnum=2;
    yminus=E_fun(thetaminus);
    Funnum=3;
    yplus=E_fun(thetaplus);
    Funnum=0;
    
    FunValue(k+1,:)=[theta',MSE_globalmain,Timemain];
    FunValue_minus(k+1,:)=[thetaminus',MSE_globalminus,Timeminus];
    FunValue_plus(k+1,:)=[thetaplus',MSE_globalplus,Timeplus];
    FunValue_MSE(:,k+1,:)=[MSE_realz_global];  %each function evaluation before averaging it into objective function
    FunValue_Time(:,k+1)=[Time_realz_global(:)];  %each time evaluation before averaging it into objective function
    optfigures;
    
    %gradient estimation
    ghat=(yplus-yminus)./(2*ck*delta);
    stepsize=ak*ghat;
    if stage==1
        stepsize=min(stepsize,stepsizemax);
        stepsize=max(stepsize,-stepsizemax);
    end
    theta=theta-stepsize;
    if stage==1 && MSE_globalmain>(funcbeststage0+Ttradeoff*funcbeststage0)
        theta=thetaok;
    else
        thetaok=theta;
    end
	stepsizestor(k+1,:)=[stepsize'];
    
    theta=min(theta,thetamax);
    theta=max(theta,thetamin);
    theta(4:6)=theta(4:6)/sum(theta(4:6));

    %stopping criterion
    if k==0
        func1=MSE_globalmain;
        para1=theta;
    end
    if k>0
        func2=MSE_globalmain;
        para2=theta;
        if func2>func1*0.95
            counter=counter+1;
        else
            counter=0;
        end
        if counter==0
            func1=func2;
        end
    end
    if counter==5 && stage==0
        stage=1;
        counter=0;
        func1=MSE_globalmain;
        para1=theta;
        funcbeststage0=min(FunValue(1:k+1,p+1));
        
        figure(1);
        hold on
        subplot(3,2,1)
        line([k,k],[0,thetamax(1)],'color','red')
        subplot(3,2,3)
        line([k,k],[0,thetamax(2)],'color','red')
        subplot(3,2,5)
        line([k,k],[0,thetamax(3)*100],'color','red')
        subplot(3,2,2)
        line([k,k],[0,thetamax(2)],'color','red')
        subplot(3,2,4)
        line([k,k],[0,thetamax(4)],'color','red')
        subplot(3,2,6)
        line([k,k],[0,thetamax(6)],'color','red')
        figure(2);
        hold on
        line([k,k],[0,MSE_globalmain+0.2],'color','red')
        line([k,k+1],[(Ttradeoff+100)/100*funcbeststage0,(Ttradeoff+100)/100*funcbeststage0],'color','red')
        figure(3);
        hold on
        line([k,k],[0,Timemain+1],'color','red')
    end
    if counter==5 && stage==1
        break
    end    
end
    FunValue(k+2,:)=[theta',MSE_globalmain,Timemain];
    FunValue_minus(k+2,:)=[thetaminus',MSE_globalminus,Timeminus];
    FunValue_plus(k+2,:)=[thetaplus',MSE_globalplus,Timeplus];
    FunValue_MSE(:,k+2,:)=[MSE_realz_global];  %each function evaluation before averaging it into objective function
    FunValue_Time(:,k+2)=[Time_realz_global(:)];  %each time evaluation before averaging it into objective function
    theta
if nsim>1
    figure(4);
    boxplot(FunValue_MSE(:,1:k+1),'Labels',{0:k});
    hold on
    line(1:k+1,FunValue(1:k+1,p+1),'color','black')
    xlim([0, k+2]), ylim([0 inf])
    figure(5);
    boxplot(FunValue_Time(:,1:k+1),'Labels',{0:k});
    hold on
    line(1:k+1,FunValue(1:k+1,p+2),'color','black')
    xlim([0, k+2]), ylim([0 inf])
end
for H=1:5
    filename=['fig',num2str(H),'-','nsim',num2str(nsim),'-TO',num2str(Ttradeoff)];
    savefig(H,filename)
end
save(['nsim',num2str(nsim),'-TO',num2str(Ttradeoff)])