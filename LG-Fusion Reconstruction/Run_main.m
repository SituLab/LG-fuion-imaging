%% Code for LG-fusion reconstruction
% Author: Haiming Yuan (hmyuan@siom.ac.cn)
% corresponding to the Fig. 1,3 of manuscript.
% Non-commercial use, copying, and modification of this code is permitted only when you cite our paper.


%Reference:
%Lei Zhu , Fernando Soldevila , Claudio Moretti , Alexandra d’Arco , Antoine Boniface ,
% Xiaopeng Shao , Hilton B. de Aguiar,et al.2022."Large field-of-view non-invasive 
% imaging through scattering layers using fluctuating random illumination".
% Nature Communications,13(1):.https://doi.org/https://doi.org/10.1038/s41467-022-29166-y

%% Input
clear all, close all, clc;
imgnum=2;  %1:pi   2:smile

% Ask for human input to control the object to be reconstructed
%% Parameter
dat = strings(1,2);
dat(1)='datapi.mat';
dat(2)='datasmile.mat';
EstiRank=[21,30];%Rank estimation of pi and smile
imamax=500; %quantity
rexx=100; %speckle image scaling

%% Load enhanced speckle data
load(strcat('./data/',dat(imgnum)),'datapool');
datapoolsum=datapool;
T=imamax;
datapool=datapoolsum(1:T,:);

%% Rank estimation
% rankmax=35;
% repeat=5;
% Resi=zeros(repeat,rankmax);
% for i=1:rankmax
%     for j=1:repeat
%     [W,H,resi]=nnmf(datapool,i,'algorithm','als');
%     Resi(j,i)=resi;
%     disp(i);
%     disp(resi);
%     end
% end
% resimin=min(Resi);
% resimax=max(Resi);
% resimean=sum(Resi)/repeat;
% f = figure("Name","Residual");
% Ed=resimean-resimin;
% Eu=resimax-resimean;
% Ex=1:rankmax;
% e = errorbar(Ex,resimean,Ed,Eu);
% legend('mean and error of 5 realizations')
% xlabel('Rank','FontSize',10)
% ylabel('Residual','FontSize',10)
% e.LineWidth = 1;  
% e.CapSize = 10; 


%% NMF
EstimatedNumber=EstiRank(imgnum);
iter=0;
while iter<20
opt = statset('MaxIter',30,'Display','final');
[W0,H0] = nnmf(datapool,EstimatedNumber,'Replicates',10,...
    'options',opt,'algorithm','mult');

opt = statset('MaxIter',1000,'Display','iter','TolFun',1e-6);
[W,H,resi,iter] = nnnmf(datapool,EstimatedNumber,'W0',W0,'H0',H0,...
    'options',opt,'algorithm','als');
end



%% Fingerprint PSFs reshape
xpixel = rexx;
ypixel = rexx;
M = cell(EstimatedNumber,1);
figure;
set(gcf,'position',[0,0,600,600]);
ha = tight_subplot(floor(sqrt(EstimatedNumber)+1),floor(sqrt(EstimatedNumber)+1),EstimatedNumber,[.01 .01],[.01 .01],[.01 .01]);
for kk=1:EstimatedNumber
    M{kk} = reshape(H(kk,:),xpixel,ypixel);
    axes(ha(kk));
    imagesc(M{kk}), colormap hot;
    text(5,5,['PSF#',num2str(kk)], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize',10, 'Color', 'white');
    axis square;
    axis off;
end

%% Intensity fluctuations    
Wi=W'; 
hymax=max(max(Wi));
figure;
tiledlayout(floor(sqrt(EstimatedNumber)+1),floor(sqrt(EstimatedNumber)+1));
for kk=1:EstimatedNumber  
    nexttile
    plot(Wi(kk,:))
    title(['#',num2str(kk)])
    axis([0,500,0,hymax]);
end

%% Deconvolution for vector distance between two object points
mu = 20;
opts.rho_r   = 1;
opts.rho_o   = 1;
opts.beta    = [1.0 1.0 0];
opts.gamma  = 2;
opts.print   = false;
opts.alpha   = 0.1;
opts.tol   = 1e-5;
opts.method  = 'l2';
opts.max_itr  = 100;

recon_image = zeros(xpixel,ypixel,EstimatedNumber);
Maximum_intensity = zeros(EstimatedNumber,EstimatedNumber);
xx = zeros(EstimatedNumber,EstimatedNumber);
yy = zeros(EstimatedNumber,EstimatedNumber);

for k=1:EstimatedNumber
    PSF = M{k};
    for i = 1:EstimatedNumber
        out = deconvtv(M{i},PSF, mu, opts);
        recon_image(:,:,i) = out.f;
        [xx_temp,yy_temp] = find(recon_image(:,:,i) == max(max(recon_image(:,:,i))));      
        if size(xx_temp,1)>1 || size(yy_temp,1)>1
            xx(i,k) = NaN;
            yy(i,k) = NaN;
        else
            xx(i,k) = floor(mean(xx_temp(:)));
            yy(i,k) = floor(mean(yy_temp(:)));
        end
        Maximum_intensity(i,k) = max(max(recon_image(:,:,i)));
    end
end

%% Delete results that exceed MER
WW=sum(W',2);
Rrco=zeros(EstimatedNumber,EstimatedNumber,3);
Rrco(:,:,1)=Maximum_intensity;
Rrco(:,:,2)=xx-50;
Rrco(:,:,3)=yy-50;
for i=1:EstimatedNumber
    for j=1:EstimatedNumber
        if Rrco(i,j,1)<0.01
            Rrco(i,j,1)=0;
            Rrco(i,j,2)=0;
            Rrco(i,j,3)=0;
        else
            Rrco(i,j,1)=1;
        end     
    end
end
for i=1:EstimatedNumber
    for j=1:EstimatedNumber
        if Rrco(i,j,2)==0
            Rrco(i,j,2)=-Rrco(j,i,2);
            Rrco(i,j,1)=Rrco(j,i,1);
        end
        if Rrco(i,j,3)==0
            Rrco(i,j,3)=-Rrco(j,i,3);
            Rrco(i,j,1)=Rrco(j,i,1);
        end
    end
end
%% Delete invalid results from a reserved noise layer(see Supplementary 6)
changed = true;
while changed
    changed = false;
    if any(sum(abs(Rrco(:,:,1)), 2) == 0) || any(sum(abs(Rrco(:,:,1)), 1) == 0)       
        WW=WW(sum(abs(Rrco(:,:,1)), 2) > 0);
        Rrco = Rrco(sum(abs(Rrco(:,:,1)), 2) > 0, :, :);
        Rrco = Rrco(:, sum(abs(Rrco(:,:,1)), 1) > 0, :);
        changed = true;
    end
end
EstimatedNumber2=size(Rrco,1);
Rrco2=zeros(EstimatedNumber2,EstimatedNumber2,2);
deco2=Rrco(:,:,1);
Rrco2(:,:,1)=Rrco(:,:,2);
Rrco2(:,:,2)=Rrco(:,:,3);

%% Undirected Graph and Minimum Spanning Tree
deco3 = deco2 - diag(diag(deco2));
G = graph(deco3);
figure;
plot(G, 'Layout', 'force');
title('Undirected Graph');
% If two points are within the memory effect range link with a line segment
Tmin = minspantree(G);
deco4 = adjacency(Tmin);
G = graph(deco4);
figure;
plot(G, 'Layout', 'force');
title('Minimum Spanning Tree');
% The most concise relationship between object points

%% Breadth First Search
P=zeros(EstimatedNumber2,2);
P(1,:)=[50 50];
pl=zeros(EstimatedNumber2,1);
hl=zeros(EstimatedNumber2,1);
i=1;
while 1
    for j=1:EstimatedNumber2
        if pl(j)==0
            P(j,1)=P(i,1)+Rrco2(i,j,1);
            P(j,2)=P(i,2)+Rrco2(i,j,2);
            if deco2(i,j)==0
                pl(j)= pl(j) || 0 ;
            else
                pl(j)=1;
            end
        end
    end
    hl(i)=1;
    nl=abs(pl-hl);
    nonZeroIndices = find(nl); 
    if ~isempty(nonZeroIndices)
        i = nonZeroIndices(1); 
    else
        disp('finished');
        break
    end   
end

%% LG-fusion
imag=zeros(100,100);
for i=1:EstimatedNumber2
    imag(P(i,1),P(i,2))=WW(i);
end 
figure
imshow(mat2gray(imrotate(imag,180)));
colormap hot






 