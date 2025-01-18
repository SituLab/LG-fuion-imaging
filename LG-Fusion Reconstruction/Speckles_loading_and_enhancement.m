%% Speckles loading and enhancement from raw data

doc='1';
dat='./data/';
rexx=100;
dx=1000;
imamax=500;
datapool=zeros(rexx*rexx,imamax);
allimage=fullfile(strcat(dat,doc));
searchname=dir(fullfile(allimage,'*.tif'));
name={searchname.name};

for kk=1:imamax
    imageinit=imread(strcat(strcat(dat,doc,'/'),name{1,kk}));
    I=double(imresize((imageinit(1:dx,1:dx)),[rexx rexx],'nearest'));
    I=I-min(min(I));
    s=fftshift(fft2(I));
    absa=mat2gray(abs(s));
    [a,b]=size(s);
    a0=round(a/2);
    b0=round(b/2);
    d=10;
    p=0.2;q=0.5;
    for i=1:a
        for j=1:b
            distance=sqrt((i-a0)^2+(j-b0)^2);
            if distance<=d
                h=0;
            else
                h=1;
            end
            s(i,j)=(p+q*h)*s(i,j);
        end
    end
    image=(real(ifft2(ifftshift(s))));
    image=imrotate(image,0);
    %image=imadjust(image);%binary object
    if kk<2
        figure
        imagesc(image),daspect([1 1 1]), colormap hot;
    end
    datapool(:,kk)=reshape(image,[rexx*rexx,1]);
end
datapool=datapool';
save(strcat(dat,doc,'/data.mat'),'datapool');