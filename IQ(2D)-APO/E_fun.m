    function MSE_T = E_fun(x)
%loading global variables
Global_var
gapconddis=1; % the additional pixel around the gaps
tilesize = round(x(1)*100,0); % Optimal size for the patch size is between 1/7 and 1/8 of the TI
% (To get better realizations tilesize should be minimum)
tilesize_vary = 0; % Randomize patch size (max +-10% of patch size)
overlap = round(x(2)*tilesize,0); % Optimal size for the overlap region is between 1/3 and 1/4 of the tilesize
if overlap<4
    overlap=4;
end
nbreplicates = round(x(3)*10,0); % Standard Number of replicas = 5 to 10 to avoid "Verbatim Copy"
cond = 1; % Define whether Conditional or Unconditional simulation (0 or 1)
c = 1000; % No of conditioning points (Don’t need to be defined if you load your Data file in the code)
w = 0.99; % Conditioning weight (Range [0 1])
w_v = 1; % Each variable's weight (For Multivariable IQ analysis)
% This is univariate IQ version
nbrealz = nb_realz; % Number of realizations to be produced
temp_split = 1; % temp_split = 1 Perform template splitting; temp_split = 0 Do not perform template splitting
do_cut = 1; % do_cut = 1 Perform optimum cutting; do_cut = 0 Do not perform optimum cutting

show_TI = 0; % Showing the TI
show_realizations = 0; % Showing all the simulated realizations
show_etype_stdev_maps = 0; % Showing the E-type and Standard deviation maps for all variables
show_histogram = 0; % Showing the histograms
show_variogram = 0; % Showing the variograms
show_data = 0; % Showing the hard data
show_connectivity = 0; % Showing the connectivity in both X and Y directions for multi-facies TI & Realizations
%% Define the size of final realizations
% realization size is similar to TI:
X1=Ti;
[p q Nvar] = size (X1);
%% Calculate m, n from tilesize and overlap
% m: The number of tiles to be placed in the output image, in X dimension
% n: The number of tiles to be placed in the output image, in Y dimension
m = (p-overlap)/(tilesize-overlap);
m = ceil(m);
n = (q-overlap)/(tilesize-overlap);
n = ceil(n);

%% Consider the full TI
X = X1;
%% preallocating
MSE_realz=zeros(nsim,nbvar_Ti);
time_T=zeros(nsim,1);
realz=zeros(Xsize_Ti,Ysize_Ti,Zsize_Ti,nb_realz,nbvar_Ti);
ErrorT=zeros(gap_num,nbvar_Ti);
Error=zeros(gap_num,nb_realz,nbvar_Ti);
Datagap=[];
if Zsize_Ti<=gapsize
    ngap=zeros(gap_num,2);
    Ti_crop=zeros(gapsize,gapsize,Zsize_Ti,gap_num,nbvar_Ti); %Ti_crop: gap information
    realz_crop=zeros(gapsize,gapsize,Zsize_Ti,gap_num,nb_realz,nbvar_Ti); %realz_crop: gap info for realizations
else
    ngap=zeros(gap_num,3);
    Ti_crop=zeros(gapsize,gapsize,gapsize,gap_num,nbvar_Ti);
    realz_crop=zeros(gapsize,gapsize,gapsize,gap_num,nb_realz,nbvar_Ti);
end
%% Calculating Error function for nsim realizations
for nsimc=1: nsim
    Ti_new=Ti;
    for nsimc=1: nsim
    Ti_new=Ti;
    if Zsize_Ti<=gapsize
        % all posible gaps for 2D case
        [XX,YY]=meshgrid(1+gapdistance:Xsize_Ti-gapsize-gapdistance+1,1+gapdistance:Ysize_Ti-gapsize-gapdistance+1);
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
                biggap=Ti_new(ngap(i,1)-gapconddis:(ngap(i,1)+gapsize+gapconddis-1),ngap(i,2)-gapconddis:(ngap(i,2)+gapsize+gapconddis-1),:,h);
                [XXD,YYD]=meshgrid(ngap(i,1)-gapconddis:(ngap(i,1)+gapsize+gapconddis-1),ngap(i,2)-gapconddis:(ngap(i,2)+gapsize+gapconddis-1));
                Condgap=[XXD(:) YYD(:) biggap(:)];
                Datagap=[Datagap;Condgap];     
            end
        end
    end 
%% Produce data set from TI with randomized locations or read given data set
if cond > 0
    %     fprintf('Produce data set from a unconditional realization with randomized locations:\n');
    %     X_data = imagequilt_Unconditional_v12(X, m, n, tilesize, overlap, 10, w_v,do_cut);
    %     X_data = X_data(1:p,1:q);
    %     Data = Produce_data(c, X_data);
    fprintf('Produce data set from TI with randomized locations:\n');
    
%     Datagap=datasample(Datagap,1000,'Replace',false);
    Data = Produce_data(c, Ti); % Produce data set from TI with randomized locations
%     Data(ismember(Data,Datagap,'rows'),:)=[];
%     Data=[Data;Datagap];
%     Data(Data(:,3)==-999999999,:)=[];
    D{1,1} = Data;
end
end
    
%% Simulation
EfunTime=tic;
% Run the IQ function to get realizations
Y=zeros(p,q,nbrealz); % preallocate - 3D array to store all realizations in the 3rd dimension

if cond > 0
    fprintf('Required time for Conditioning Image Quilting:\n');
else
    fprintf('Required time for Unconditioning Image Quilting:\n');
end

I = 0;
% loop all realizations
for i=1:nbrealz
    fprintf(['Output realization #',num2str(i),': ']);
    if tilesize_vary > 0
        I = I + 1;
        if I <= tilesize_vary
            tilesize_new = tilesize + I;
        else
            I = I - tilesize_vary;
            tilesize_new = tilesize;
        end
    else
        tilesize_new = tilesize;
    end
    
    if cond > 0
        % generate realization at the size needed for the tiles arrangement, possibly too large
        Y1 = imagequilt_multivariate_v7(X, m, n, D, w, tilesize_new, overlap, nbreplicates, 1, temp_split);
    else
        Y1 = imagequilt_Unconditional_v12(X, m, n, tilesize_new, overlap, nbreplicates, w_v,do_cut);
    end
    Y(:,:,i) = Y1(1:p,1:q); % crop the realization to its true dimensions and store in 3D array Y
end
sim_time=toc(EfunTime);

    for l=1:nb_realz
        realz(:,:,:,l,:)=Y(:,:,l);
    end
    time_T(nsimc,:)=sim_time;
    
    %Error calculation for each simulation               
    if Zsize_Ti<=gapsize
        for h=1:nbvar_Ti
            for l=1:nb_realz
                for i=1:gap_num
                    realz_crop(:,:,:,i,l,h)=realz(ngap(i,1):(ngap(i,1)+gapsize-1),ngap(i,2):(ngap(i,2)+gapsize-1),:,l,h);
                end
            end
        end
    else
        for h=1:nbvar_Ti
            for l=1:nb_realz
                for i=1:gap_num
                    realz_crop(:,:,:,i,l,h)=realz(ngap(i,1):(ngap(i,1)+gapsize-1),ngap(i,2):(ngap(i,2)+gapsize-1),ngap(i,3):(ngap(i,3)+gapsize-1),l,h);
                end
            end
        end
    end
    for h=1:nbvar_Ti
        for l=1:nb_realz
            for i=1:gap_num
                if TiType(h)==0
                    E1=realz_crop(:,:,:,i,l,h);
                    E2=Ti_crop(:,:,:,i,h);
                    E1~=E2;
                    Egap=ans;
                    Error(i,l,h)=sqrt(mean((Egap(:)).^2));
                else
                    E1=realz_crop(:,:,:,i,l,h);
                    E2=Ti_crop(:,:,:,i,h);
                    Error(i,l,h)=sqrt(mean((E1(:)-E2(:)).^2));
                end
            end
        end
        ErrorT(:,h)=mean(Error(:,:,h),2);
        MSE_realz(nsimc,h)=mean(ErrorT(:,h));
     end
end
%% Averaging errors for all simulations and calculating the total error
for h=1:nbvar_Ti
    SUM_E(1,h)=mean(MSE_realz(:,h));
    A=Ti(:,:,:,h);
    Norm_E(1,h)=SUM_E(1,h)/iqr(A(:));
end
    if Funnum==1
        MSE_globalmain=sum(Norm_E);
        %MSE_globalmain=Norm_E(1); %In multivariate cases the error for a specific variate can be chosen!
        Timemain=mean(time_T);
        for h=1:nbvar_Ti
            A=Ti(:,:,:,h); %A=Ti(:,:,:,h); for multiple variate error
            MSE_realz_global=MSE_realz./iqr(A(:));
            %MSE_realz_global=MSE_realz(:,1)./iqr(A(:));%In multivariate cases the error for a specific variate can be chosen!
        end    
        Time_realz_global=time_T;
    end
    if Funnum==2
        MSE_globalminus=sum(Norm_E);
        %MSE_globalminus=Norm_E(1); %In multivariate cases the error for a specific variate can be chosen!
        Timeminus=mean(time_T);
    end
    if Funnum==3
        MSE_globalplus=sum(Norm_E);
        %MSE_globalplus=Norm_E(1); %In multivariate cases the error for a specific variate can be chosen!
        Timeplus=mean(time_T);
    end
if stage==0
    MSE_T=sum(Norm_E);
    %MSE_T=Norm_E(1); %In multivariate cases the error for a specific variate can be chosen!
else
    MSE_T=mean(time_T);
end