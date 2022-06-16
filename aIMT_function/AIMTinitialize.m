function [aorta_coor, lab_aorta]=AIMTinitialize(xorig, xroi,pars, dbf)

dim=pars.dim;
sigma=pars.sigma;
theta=pars.theta;

%% Filter bank
xim=imresize(xorig,pars.resize);
BWroi=imresize(xroi,pars.resize);
[imx,imy]=meshgrid(-dim:dim,-dim:dim);
xfilt=zeros(size(xim,1),size(xim,2),length(sigma),length(theta));
xf=(xim-mean(xim(find(BWroi)))).*BWroi;
for ct=1:length(sigma)
    %K=1/sigma(ct)^2.*(imy.^2/sigma(ct)^2-1).*exp(-0.5*imy.^2/(sigma(ct))^2-0.5*imx.^2/(dim/3)^2);
    %K=exp(-0.5*imx.^2/sigma(ct)^2)-exp(-0.5*(imx-sigma(ct)).^2/sigma(1)^2)-exp(-0.5*(imx+sigma(ct)).^2/sigma(1)^2);
    
    K1=exp(-0.5*imy.^2/sigma(ct)^2);
    K2=exp(-0.5*(imy+1.5*sigma(ct)).^2/(0.25*sigma(ct))^2);
    K3=exp(-0.5*(imy-1.5*sigma(ct)).^2/(0.25*sigma(ct))^2);
    K=(-K1+K2+K3).*exp(-0.5*imx.^2/(dim/3)^2);
    n=find(K>0);
    K(n)=K(n)/sum(K(n));
    n=find(K<0);
    K(n)=K(n)/abs(sum(K(n)));
    %K=0.5-K;
    %K=K-sum(K(:))/(2*dim+1)^2;
    for ctt=1:length(theta)
        xfilt(:,:,ct,ctt)=conv2(xf,imrotate(K,theta(ctt),'crop'),'same');
    end;
end;

%% ------------------------------------------------------------------------------
m1=max(xfilt,[],4);
[m,xmax]=max(m1,[],3);
mnorm=m-min(m(:));
mnorm=mnorm/max(mnorm(:));
th=graythresh(mnorm(BWroi));
mask=mnorm>1.25*th;
mask=mask.*BWroi;
stats=regionprops(bwlabel(mask),'Eccentricity','Centroid','Area');
CenEcc=zeros(length(stats),2);
for ct=1:length(stats)
    tmp=stats(ct,1).Centroid;
    %CenEcc(ct,1)=abs(size(xim,1)/2-tmp(2));
    CenEcc(ct,2)=stats(ct,1).Eccentricity;
end

%okCen=find(CenEcc(:,1)<median(CenEcc(:,1)));
okEcc=find(CenEcc(:,2)>0.9);
okCenEcc=okEcc;%intersect(okCen,okEcc);

%% Region identification
%scarto quelli che non stanno in mezzo o che sono poco eccentrici. Se ne
%restano più di uno tengo quello di area massima.
if length(okCenEcc)>1
    areamax=stats(okCenEcc(1),1).Area;
    cttmp=1;
    for ct=2:length(okCenEcc)
        if stats(okCenEcc(ct),1).Area>areamax
            areamax=stats(okCenEcc(ct),1).Area;
            cttmp=ct;
        end
    end
    
    okCenEcc=okCenEcc(cttmp);
end

if ~isempty(okCenEcc)
    lab_aorta=(imerode(bwlabel(mask)==okCenEcc,ones(3)));
end
lab_aorta=imresize(lab_aorta.*imerode(BWroi,ones(3)),1/pars.resize);
lab_aorta=imdilate(imopen(lab_aorta>0,ones(3)),ones(3));
aorta_coor=bwboundaries(lab_aorta);

