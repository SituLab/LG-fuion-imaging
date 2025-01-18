%% Code for Reconstruction-capability function and Reconstruction-benefit function
% Author: Haiming Yuan (hmyuan@siom.ac.cn)
% corresponding to the Fig. 1,3 of manuscript.
% Non-commercial use, copying, and modification of this code is permitted only when you cite our paper.

%Reference: Zhu, L. et al. Repository for ’Large field-of-view non-invasive imaging through
%scattering layers using fluctuating random illumination’ (v1.0.0). Zenodo.
%https://doi.org/10.5281/zenodo.5850465 (2021).

%% Speckles reading and enhancement
% Spe = cell(17,2);
% num=1;
% for sdoc=1:0.2:3 
%     doc=num2str(sdoc);
%     dat='./data20240517/';
%     rexx=100;
%     dx=1000;
%     imamax=1500;
%     datapool=zeros(rexx*rexx,imamax);
%     allimage=fullfile(strcat(dat,doc));
%     searchname=dir(fullfile(allimage,'*.tif'));
%     name={searchname.name};
%     
%     for kk=1:imamax
%         imageinit=imread(strcat(strcat(dat,doc,'/'),name{1,kk}));
%         I=double(imresize((imageinit(1:dx,1:dx)),[rexx rexx],'nearest'));
%         I=I-min(min(I));
%         %    image=I;
%         %    I=mat2gray(I);%%%%%%%%%
%         s=fftshift(fft2(I));
%         absa=mat2gray(abs(s));
%         [a,b]=size(s);
%         a0=round(a/2);
%         b0=round(b/2);
%         d=10;
%         p=0.2;q=0.5;
%         for i=1:a
%             for j=1:b
%                 distance=sqrt((i-a0)^2+(j-b0)^2);
%                 if distance<=d
%                     h=0;
%                 else
%                     h=1;
%                 end
%                 s(i,j)=(p+q*h)*s(i,j);
%             end
%         end
%         image=(real(ifft2(ifftshift(s))));
%         image=mat2gray(image);
%         image=imadjust(image);
%         image=imrotate(image,0);
%         %     figure
%         %     imshow(mat2gray(image));
%         %     if kk<2
%         %         figure
%         %         imagesc(image),daspect([1 1 1]), colormap hot;
%         %     end
%         datapool(:,kk)=reshape(image,[rexx*rexx,1]);
%     end
%     datapool=datapool';
%     Spe{num,1}=datapool;
%     Spe{num,2}=sdoc/10;
%     num=num+1;
% end
% 
% for sdoc=4:9
%     
%     doc=num2str(sdoc);
%     dat='./data20240517/';
%     rexx=100;
%     dx=1000;
%     imamax=1500;
%     datapool=zeros(rexx*rexx,imamax);
%     allimage=fullfile(strcat(dat,doc));
%     searchname=dir(fullfile(allimage,'*.tif'));
%     name={searchname.name};
%     
%     for kk=1:imamax
%         imageinit=imread(strcat(strcat(dat,doc,'/'),name{1,kk}));
%         I=double(imresize((imageinit(1:dx,1:dx)),[rexx rexx],'nearest'));
%         I=I-min(min(I));
%         %    image=I;
%         %    I=mat2gray(I);
%         s=fftshift(fft2(I));
%         absa=mat2gray(abs(s));
%         [a,b]=size(s);
%         a0=round(a/2);
%         b0=round(b/2);
%         d=10;
%         p=0.2;q=0.5;
%         for i=1:a
%             for j=1:b
%                 distance=sqrt((i-a0)^2+(j-b0)^2);
%                 if distance<=d
%                     h=0;
%                 else
%                     h=1;
%                 end
%                 s(i,j)=(p+q*h)*s(i,j);
%             end
%         end
%         image=(real(ifft2(ifftshift(s))));
%         image=mat2gray(image);
%         image=imadjust(image);
%         image=imrotate(image,0);
%         %     figure
%         %     imshow(mat2gray(image));
%         %     if kk<2
%         %         figure
%         %         imagesc(image),daspect([1 1 1]), colormap hot;
%         %     end
%         datapool(:,kk)=reshape(image,[rexx*rexx,1]);
%     end
%     datapool=datapool';
%     Spe{num,1}=datapool;
%     Spe{num,2}=sdoc/10;
%     num=num+1;
% end
% save(strcat(dat,'/alldata.mat'));


%% Image reconstruction for different encoding sparsities and quantities.
% load('alldata.mat');
% N=[30 40 50 60 70 100 150 300 500];%quantities
% spa=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];%sparsity labels
% Spanum=size(spa,2);
% Spa=zeros(1,Spanum);%sparsity
% 
% Nnum=size(N,2);
% Spanum=size(Spa,2);
% rep=3;
% Reconf=zeros(Spanum,Nnum,rep);%Reconstruction capability
% Reconiter=zeros(Spanum,Nnum,rep);%Number of iterations
% Reconimg=cell(Spanum,Nnum,rep);%Reconstruction image 
% 
% for spakk=1:Spanum
%     datapoolsum=Spe{spa(spakk),1};
%     tic
%     for nkk=1:Nnum
%         T=N(nkk);
%         datapool=datapoolsum(1:T,:);
%         for repkk=1:rep
%             [Reconf(spakk,nkk,repkk),Reconiter(spakk,nkk,repkk),Reconimg{spakk,nkk,repkk}]=nmfdeconv(datapool);
%             fprintf('sparse:%3.2f N:%d ff:%4.1f \n',Spa(1,spakk),T,Reconf(spakk,nkk,repkk));
%         end
%     end
%     toc
% end
% save('Recon3.6.mat')


%% Results visualization
V=3.6;
load(['Recon',num2str(V),'.mat'])
for i=1:Spanum
    Spa(1,i)=1-Spe{spa(i),2};
end                                                                                                                                                           
% Reconstruction-capability function
ssize=14;
[X, Y] = meshgrid(N,Spa);
reconf=max(Reconf,[],3);
reconiter=max(Reconiter,[],3);
fg=figure 
surf(X, Y,reconf,'FaceColor','interp');
xlabel('Encoding quantity N','FontSize',ssize)
ylabel('Encoding sparsity \theta','FontSize',ssize)
zlabel(' R(N,\theta)','FontSize',ssize)
colormap jet
axis([30,500,0.2,0.9,0,20]);
set(gca, 'XTick',[30,60,100,150,300,500])
ax=gca;
ax.XAxisLocation='bottom';
ax.YAxisLocation='right';
cbar=colorbar(ax,'Position',...
    [0.158035714285714 0.625476190476197 0.0616071428571408 0.230476190476192],...
    'Color',[1 1 1]);
cbar.TickLabels={};
cbar.Label.String={'max','','','','','min'};
cbar.Label.Rotation=0;
cbar.Label.Position=[1.35018315018315 1.52578130210797 0];
annotation(fg,'textbox',...
    [0.460821428571428 0.795238096905608 0.142857139131853 0.0797619030943939],...
    'Color',[1 1 1],...
    'String',{'R(\theta,N)'},...
    'LineWidth',1,...
    'LineStyle','none',...
    'FontSize',14,...
    'EdgeColor',[1 1 1]);
view(270, 90)

% Reconstruction-benefit function
frecon=reconf./X;
fg2=figure; 
surf(X, Y,frecon,'FaceColor','interp');
xlabel('Encoding quantity N','FontSize',ssize)
ylabel('Encoding sparsity \theta','FontSize',ssize)
zlabel('f(\theta,N)/N','FontSize',ssize)
colormap jet
colorbar
axis([30,500,0.2,0.9,0,0.3]);
set(gca, 'XTick',[30,60,100,150,300,500])
ax=gca;
ax.XAxisLocation='bottom';
ax.YAxisLocation='right';
cbar=colorbar(ax,'Position',...
    [0.158035714285714 0.625476190476197 0.0616071428571408 0.230476190476192],...
    'Color',[1 1 1]);
cbar.TickLabels={};
cbar.Label.String={'max','','','','','min'};
cbar.Label.Rotation=0;
cbar.Label.Position=[1.47632850241546 0.0127149165045357 0];
annotation(fg2,'textbox',...
    [0.460821428571428 0.795238096905608 0.142857139131853 0.0797619030943939],...
    'Color',[1 1 1],...
    'String',{'f(\theta,N)'},...
    'LineWidth',1,...
    'LineStyle','none',...
    'FontSize',14,...
    'EdgeColor',[1 1 1]);
view(270, 90)

% Reconstruction image display
Num=60;
sparse=0.78;
Sd=find(abs(Spa-sparse)<0.001);
Nd=find(N == Num);
figure
for jjj=1:3
    GGG=Reconimg{Sd,Nd,jjj};
    subplot(1,3,jjj)
    imshow(GGG,[]);colormap hot;
end
Reconf(Sd,Nd,:)

figure
subplot('Position', [0, 0, 1, 1]);
set(gcf,'position',[0,0,300,300]);
jjj=3;
GGG=Reconimg{Sd,Nd,jjj};
dx=4;
dy=0;
imshow(GGG(50+(-16:10)+dx,50+(-16:10)+dy),[]);colormap hot;
fprintf('%3.2f %d %d %d %d\n,',sparse,Num,jjj,dx,dy);


