function [ff,iter,GGG] = nmfdeconv(datapool)  
% UseParallel NMF and deconvlution


%% UseParallel NMF
EstimatedNumber =19; 
rexx=100;

opt = statset('MaxIter',30,'Display','off');
[W0,H0] = nnmf(datapool,EstimatedNumber,'Replicates',5,...
    'options',opt,'algorithm','mult');

opt = statset('MaxIter',1000,'Display','off','TolFun',1e-6,'UseParallel',1);
[W,H,resi,iter] = nnnmf(datapool,EstimatedNumber,'W0',W0,'H0',H0,...
    'options',opt,'algorithm','als');

%% Fingerprint PSF reshape
global xpixel ypixel
xpixel = rexx; % the size of the image in the X direction
ypixel = rexx; % the size of the image in the Y direction
M = cell(EstimatedNumber,1);% create the empty matrix to store the fingerprint

for kk=1:EstimatedNumber
    M{kk} = reshape(H(kk,:),xpixel,ypixel);
end

%% Deconvlution
mu = 20;
opts.rho_r   = 2;
opts.rho_o   = 50;
opts.beta    = [1.0 1.0 0];
opts.gamma  = 2;
opts.print   = false;
opts.alpha   = 0.7;
opts.tol   = 1e-5;
opts.method  = 'l2';
opts.max_itr  = 100;
O = cell(EstimatedNumber,1);
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
    O{k}=sum(recon_image,3);
end
%% Merge images
First_pattern = 2;
[Reached_Pattern,Global_image,ff1] = mergeimages(Maximum_intensity,O,EstimatedNumber,xx,yy,xpixel,ypixel,First_pattern);
[Reached_Pattern,Global_image,ff2] = mergeimages(Maximum_intensity,O,EstimatedNumber,xx,yy,xpixel,ypixel,Reached_Pattern(end));
ff=max(ff1,ff2);

Global_image=mat2gray(Global_image);
GGG=Global_image(250+(-50:49),250+(-50:49));

end

