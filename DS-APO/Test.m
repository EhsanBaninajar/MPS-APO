clc
clear all
%% Loading global variables
Global_var
gapsize=10;
%% Loading Ti
%Training Image
[Ti,Xsize_Ti,Ysize_Ti,Zsize_Ti,nbvar_Ti,namevar_Ti]=LoadGrid('MultiVariate.sgems');
gapdatabase
%% Reading the parameters file from source
params=ReadParamsFile('params_read.dat');
%% params values
x=[0.1 0.5 0.2];
nb_realz=params.nb_realz;
params.nb_cells_x=Xsize_Ti;
params.nb_cells_y=Ysize_Ti;
params.nb_cells_z=Zsize_Ti;
nvalue=fix(x(3)*100);
params.t=x(1);%t=x(1)
params.f=x(2);%f=x(2)
params.n=nvalue.*ones(1,nbvar_Ti);%n=x(3)
% params.var_weight=[x(4), x(5), x(6)];%w1=x(4) w2=x(5) w3=x(6)
%% writing the parameters file to params
WriteParamsFile(params,'params.dat');
%% Reading the parameters file from source
params=ReadParamsFile('params_read.dat');
%% params values
nb_realz=params.nb_realz;
params.nb_cells_x=Xsize_Ti;
params.nb_cells_y=Ysize_Ti;
params.nb_cells_z=Zsize_Ti;
nvalue=fix(x(3)*100);
params.t=x(1);%t=x(1)
params.f=x(2);%f=x(2)
params.n=nvalue.*ones(1,nbvar_Ti);%n=x(3)
%params.var_weight=[x(4), x(5), x(6)];%w1=x(4) w2=x(5) w3=x(6)
%% writing the parameters file to params
WriteParamsFile(params,'params.dat');
%% preallocating
MSE_realz=zeros(nsim,nbvar_Ti);
time_T=zeros(nsim,1);
realz=zeros(Xsize_Ti,Ysize_Ti,Zsize_Ti,nb_realz,nbvar_Ti);
ErrorT=zeros(1,nb_realz,nbvar_Ti);
if Zsize_Ti<=gapsize
    ngap=zeros(gap_num,2);
    Ti_crop=zeros(gapsize,gapsize,Zsize_Ti,gap_num,nbvar_Ti); %Ti_crop: gap information
    realz_crop=zeros(gapsize,gapsize,Zsize_Ti,gap_num,nb_realz,nbvar_Ti); %realz_crop: gap info for realizations
else
    ngap=zeros(gap_num,3);
    Ti_crop=zeros(gapsize,gapsize,gapsize,gap_num,nbvar_Ti);
    realz_crop=zeros(gapsize,gapsize,gapsize,gap_num,nb_realz,nbvar_Ti);
end
%% 
%number of simulations
nsim=1;
%number of the gaps to be calculated
gap_num=10;
%%
%preallocating
MSE_realz=zeros(nsim,nbvar_Ti);
time_T=zeros(nsim,1);
ngap=zeros(gap_num,2);
realz=zeros(Xsize_Ti,Ysize_Ti,Zsize_Ti,nb_realz,nbvar_Ti);
ErrorT=zeros(gap_num,nb_realz,nbvar_Ti);
if Zsize_Ti<gapsize
    Ti_crop=zeros(gapsize,gapsize,Zsize_Ti,gap_num,nbvar_Ti); %Ti_crop: gap information
    realz_crop=zeros(gapsize,gapsize,Zsize_Ti,gap_num,nb_realz,nbvar_Ti); %realz_crop: gap info for realizations
else
    Ti_crop=zeros(gapsize,gapsize,gapsize,gap_num,nbvar_Ti);
    realz_crop=zeros(gapsize,gapsize,gapsize,gap_num,nb_realz,nbvar_Ti);
end
%%
%Calculating Error function for nsim realizations
for nsimc=1: nsim
    Ti_new=Ti;
    if Zsize_Ti<=gapsize
        % all posible gaps for 2D case
        [XX,YY]=meshgrid(1:Xsize_Ti-gapsize+1,1:Ysize_Ti-gapsize+1);
        gap_pos=[XX(:),YY(:)]; %gap_pos: coordination of gap bottom left corner
        % selecting n gap without overlap
        ii=1;
        shuffle=randperm(top_gaps);
        kk=0;
        while and((ii<=gap_num),(kk<top_gaps))
            kk=kk+1;
            random_gap=sorted_std_index(shuffle(kk));
            npos=[gap_pos(random_gap,1),gap_pos(random_gap,2)];
            isOk=true;
            for j=1:ii-1
                ndif=max(abs(npos(1)-ngap(j,1)),abs(npos(2)-ngap(j,2)));
                if(ndif<=2*gapsize)
                    isOk=false;
                    break;
                end
                end
            if isOk
                ngap(ii,:)=npos;
                ii=ii+1;
            end
        end
        % save the gap info in Ti_crop and delete it in Ti_new
        for h=1:nbvar_Ti
            for i=1:gap_num
                Ti_crop(:,:,:,i,h)=Ti(ngap(i,1):(ngap(i,1)+gapsize-1),ngap(i,2):(ngap(i,2)+gapsize-1),:,h);
                Ti_new(ngap(i,1):(ngap(i,1)+gapsize-1),ngap(i,2):(ngap(i,2)+gapsize-1),:,h)=-999999999;
            end
        end
    else
        % all posible gaps for 3D case
        [XX,YY,ZZ]=meshgrid(1:Xsize_Ti-gapsize+1,1:Ysize_Ti-gapsize+1,1:Zsize_Ti-gapsize+1);
        gap_pos=[XX(:),YY(:),ZZ(:)];
        % selecting n gap without overlap
        ii=1;
        shuffle=randperm(top_gaps);
        kk=0;
        while and((ii<=gap_num),(kk<top_gaps))
            kk=kk+1;
            random_gap=sorted_std_index(shuffle(kk));
            npos=[gap_pos(random_gap,1),gap_pos(random_gap,2),gap_pos(random_gap,3)];
            isOk=true;
            for j=1:ii-1
                ndif=max(abs(npos(1)-ngap(j,1)),abs(npos(2)-ngap(j,2)),abs(npos(3)-ngap(j,3)));
                if(ndif<=2*gapsize)
                    isOk=false;
                    break;
                end
            end
            if isOk
                ngap(ii,:)=npos;
                ii=ii+1;
            end
        end
        for h=1:nbvar_Ti
            for i=1:gap_num
                Ti_crop(:,:,:,i,h)=Ti(ngap(i,1):(ngap(i,1)+gapsize-1),ngap(i,2):(ngap(i,2)+gapsize-1),ngap(i,3):(ngap(i,3)+gapsize-1),h);
                Ti_new(ngap(i,1):(ngap(i,1)+gapsize-1),ngap(i,2):(ngap(i,2)+gapsize-1),ngap(i,3):(ngap(i,3)+gapsize-1),h)=-999999999;
            end
        end
    end 
   %%
    %writing the new Ti file
    WrGrid(Ti_new,'GRID_mask.sgems',namevar_Ti);
    %simulation
    EfunTime=tic;
    !ds_x64.bat;
    sim_time=toc(EfunTime);

    for l=1:nb_realz
        str=['rec_realz_',num2str(l-1),'.sgems'];
        realz(:,:,:,:,l)=LoadGrid(str);
    end
    time_T(nsimc,:)=sim_time;
for i=1:nbvar_Ti
    figure;ViewGrid(realz(:,:,:,i))
end
figure;ViewGrid(Ti)
figure;ViewGrid(Ti_new)
%%
    %Error calculation for each simulation               
    if Zsize_Ti<=gapsize
        for l=1:nb_realz
            for i=1:gap_num
                realz_crop(:,:,:,i,l)=realz(ngap(i,1):(ngap(i,1)+gapsize-1),ngap(i,2):(ngap(i,2)+gapsize-1),:,l);
            end
        end
    else
        for l=1:nb_realz
            for i=1:gap_num
                realz_crop(:,:,:,i,l)=realz(ngap(i,1):(ngap(i,1)+gapsize-1),ngap(i,2):(ngap(i,2)+gapsize-1),ngap(i,3):(ngap(i,3)+gapsize-1),l);
            end
        end
    end
     for h=1:nbvar_Ti
        for l=1:nb_realz
            for i=1:gap_num
                if TiType==0
                    E1=realz_crop(:,:,:,i,l,h);
                    E2=Ti_crop(:,:,:,i,h);
                    E1==E2;
                    Egap=ans;
                    Error(i,l,h)=sqrt(mean((Egap(:)).^2));
                else
                    E1=realz_crop(:,:,:,i,l,h);
                    E2=Ti_crop(:,:,:,i,h);
                    Error(i,l,h)=sqrt(mean((E1(:)-E2(:)).^2));
                end
            end
        end
        ErrorT=Error(:,:,h);
        MSE_realz(nsimc,h)=mean(mean(ErrorT,2));
     end
end
%% Averaging errors for all simulations and calculating the total error
for h=1:nbvar_Ti
    SUM_E(1,h)=mean(MSE_realz(:,h));
    A=Ti(:,:,:,h);
    if TiType(h)==0
    else
    Norm_E(1,h)=SUM_E(1,h)/iqr(A(:));
    end
end
    if Funnum==1
        MSE_globalmain=sum(Norm_E);
        Timemain=mean(time_T);
        MSE_realz_global=MSE_realz;
        Time_realz_global=time_T;
    end
    if Funnum==2
        MSE_globalminus=sum(Norm_E);
        Timeminus=mean(time_T);
    end
    if Funnum==3
        MSE_globalplus=sum(Norm_E);
        Timeplus=mean(time_T);
    end
if stage==0
    MSE_T=sum(Norm_E);
else
    MSE_T=mean(time_T);
end