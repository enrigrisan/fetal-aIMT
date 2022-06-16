function [tratto_aortax,distrid]=AIMT_findStrip2(xim,xSkel,ySkel,distSkel,tratto_aortax0,nf,dbf)

step=10;

[gx,gy]=gradient(xim);
ximbordisup=zeros(5,length(xSkel));
ximbordiinf=zeros(5,length(xSkel));

minGradyneg=zeros(size(xSkel));
maxGradypos=zeros(size(xSkel));
sp=find(distSkel>mean(distSkel),1,'first');
lp=find(distSkel>mean(distSkel),1,'last');
for ct=sp:lp%4:length(xSkel)-2 %parto dal quarto perchè prima le distanze non sono corrette
    rxSkel=round(xSkel(ct));
    npix=round(distSkel(ct)/2)+step;
    npixGrLev=round(distSkel(ct)/2)+5;
    ximcol=xim(round(ySkel(ct))-npixGrLev:round(ySkel(ct))+npixGrLev,rxSkel);
    %     figure(ctframe);
    %     hold on
    %gradiente a cavallo del bordo; livelli di grigio oltre il bordo
    %         plot(rxSkel,round(ySkel(ct))-npix:round(ySkel(ct))+npix,'b.')
    %         gradcol=abs(gy(round(ySkel(ct))-npix:round(ySkel(ct))+npix,rxSkel));
    %         graycol=xim(round(ySkel(ct))-npixGrLev:round(ySkel(ct))+npixGrLev,rxSkel);
    ximbordisup(:,ct)=ximcol(1:5);
    ximbordiinf(:,ct)=ximcol(end-4:end);
    gradyneg=gy(1:round(ySkel(ct)),rxSkel);
    try
        gradyneg=gradyneg(end-npix:end-npix+step);
        
        gradypos=gy(round(ySkel(ct)):end,rxSkel);
        gradypos=gradypos(npix-step+1:npix);
        minGradyneg(ct)=min(gradyneg);
        maxGradypos(ct)=max(gradypos);
    catch
        keyboard
    end;
    %     if min(gradyneg)<-0.1
    %         testgradbordis(1,ct)=1;
    %     end
    %     if max(gradypos)>0.1
    %         testgradbordii(1,ct)=1;
    %     end
end
meanGyneg=mean(minGradyneg(4:end-2));
meanGypos=mean(maxGradypos(4:end-2));
testximbordis=zeros(size(xSkel));
testximbordii=zeros(size(xSkel));
testGyneg=zeros(size(xSkel)); %il gradiente negativo è quello sopra
testGypos=zeros(size(xSkel));


for ct=1:length(xSkel)
    %controllo che la parte esterna superiore sia chiara
    ctok=0;
    for nrighe=1:size(ximbordisup,1)
        mediariga=mean(ximbordisup(nrighe,:));
        if ximbordisup(nrighe,ct)>=mediariga
            ctok=ctok+1;
        end
    end
    if ctok>=size(ximbordisup,1)-1
        testximbordis(1,ct)=1;
    end
    %controllo che la parte esterna inferiore sia chiara
    ctok=0;
    for nrighe=1:size(ximbordisup,1)
        mediariga=mean(ximbordiinf(nrighe,:));
        if ximbordiinf(nrighe,ct)>=mediariga
            ctok=ctok+1;
        end
    end
    if ctok>=size(ximbordisup,1)-1
        testximbordii(1,ct)=1;
    end
    
    %controllo che la parte ext sup abbia gradiente molto negativo
    if minGradyneg(ct)<meanGyneg;
        testGyneg(ct)=1;
    end
    
    %controllo che la parte ext inf abbia gradiente molto positivo
    if maxGradypos(ct)>meanGypos;
        testGypos(ct)=1;
    end
    
end

prova=testximbordis+testximbordii+testGyneg+testGypos;
prova4=prova==4;
prova4conv=conv(double(prova4),ones(1,7)/7,'same');
bw4=bwlabel(prova4conv);
if max(bw4)>0
    lungh=zeros(max(bw4),1);
    for ct=1:max(bw4)
        lungh(ct)=length(find(bw4==ct));
    end
    [maxv,maxp]=max(lungh);
    tratto_aortax=xSkel(bw4==maxp);
    tratto_aortax=round(tratto_aortax(1)):round(tratto_aortax(end));
    %     if ~isempty(tratto_aortax0)
    %         if isempty(intersect(tratto_aortax,tratto_aortax0))
    %             tratto_aortax=tratto_aortax0;
    %         end
    %     end
    if dbf
        d1=(xSkel-tratto_aortax(1)).^2;
        n1=find(d1==min(d1));
        
        d2=(xSkel-tratto_aortax(end)).^2;
        n2=find(d2==min(d2));
        figure(nf); hold on; plot(xSkel(n1:n2),ySkel(n1:n2),'c.')
    end
    distrid=distSkel(bw4==maxp);
    %         meandiam(ctframe,1)=mean(distrid);
    %         stdiam(ctframe,1)=std(distrid);
else
    tratto_aortax=[];
    distrid=[];
end
