function MSE_T = E_fun(x)
%loading global variables
Global_var
%% Reading the parameters file from source
params=sgems_get_par('snesim_std');
%% params values
params.ti_file='Grid_mask.sgems';
params.XML.parameters.Nb_Realizations.value=nb_realz;   %nb_realz: number of realizations
params.XML.parameters.Nb_Facies.value=3;
params.XML.parameters.Marginal_Cdf.value=[0.3270,0.3413,0.3317];
%% Parameters to OPtimize
Max_Cond=round(x(1)*100,0);
Cmin=round(x(2)*100,0);
Nb_Multigrids_ADVANCED=round(x(3)*10,0);
params.XML.parameters.Nb_Multigrids_ADVANCED.value=Nb_Multigrids_ADVANCED;%number of multiple grids=x(1)
params.XML.parameters.Cmin.value=Cmin;%minimum replicates number=x(2)
params.XML.parameters.Max_Cond.value=Max_Cond;%maximum number of conditional data=x(3)
%% Grid used for simulation
%grid size
params.dim.nx=Xsize_Ti;
params.dim.ny=Ysize_Ti;
params.dim.nz=Zsize_Ti;
% grid cell size
params.dim.dx=1;
params.dim.dy=1;
params.dim.dz=1;
% grid origin
params.dim.x0=1;
params.dim.y0=1;
params.dim.z0=1;
%% preallocating
MSE_realz=zeros(nsim,nbvar_Ti);
time_T=zeros(nsim,1);
realz=zeros(Xsize_Ti,Ysize_Ti,Zsize_Ti,nb_realz,nbvar_Ti);
ErrorT=zeros(gap_num,nbvar_Ti);
Error=zeros(gap_num,nb_realz,nbvar_Ti);
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
                if(ndif<=gapsize)
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
                Ti_new(ngap(i,1):(ngap(i,1)+gapsize-1),ngap(i,2):(ngap(i,2)+gapsize-1),:,h)=-9966699;
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
                ndif=max(max(abs(npos(1)-ngap(j,1)),abs(npos(2)-ngap(j,2))),abs(npos(3)-ngap(j,3)));
                if(ndif<=gapsize)
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
                Ti_new(ngap(i,1):(ngap(i,1)+gapsize-1),ngap(i,2):(ngap(i,2)+gapsize-1),ngap(i,3):(ngap(i,3)+gapsize-1),h)=-9966699;
            end
        end
    end  

    % Conditional data
    [YY,XX,ZZ]=meshgrid(1:Ysize_Ti,1:Xsize_Ti,1:Zsize_Ti);
    Con_data=[XX(:) YY(:) ZZ(:) Ti_new(:)];
    params.d_obs=Con_data;
    %writing the new Ti file
    sgems_write_grid(1:Xsize_Ti,1:Ysize_Ti,1:Zsize_Ti,Ti,'Grid_mask.sgems','Ti','facies');
    %simulation
    EfunTime=tic;
    params=sgems_grid(params);
    sim_time=toc(EfunTime);

    for l=1:nb_realz
        realz(:,:,:,l)=params.D(:,:,:,l);
    end
    time_T(nsimc,:)=sim_time;
    
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
                E1=realz_crop(:,:,:,i,l,h);
                E2=Ti_crop(:,:,:,i,h);
                E1~=E2;
                Egap=ans;
                Error(i,l,h)=sqrt(mean((Egap(:)).^2));
            end
        end
        ErrorT(:,h)=mean(Error(:,:,h),2);
        MSE_realz(nsimc,h)=mean(ErrorT(:,h));
     end
end
%% Averaging errors for all simulations and calculating the total error
% SUM_E=zeros(1,nbvar_Ti);
SUM_E=mean(MSE_realz(:));
    if Funnum==1
        MSE_globalmain=sum(SUM_E);
        Timemain=mean(time_T);
        MSE_realz_global=MSE_realz;
        Time_realz_global=time_T;
    end
    if Funnum==2
        MSE_globalminus=sum(SUM_E);
        Timeminus=mean(time_T);
    end
    if Funnum==3
        MSE_globalplus=sum(SUM_E);
        Timeplus=mean(time_T);
    end
if stage==0
    MSE_T=sum(SUM_E);
else
    MSE_T=mean(time_T);
end