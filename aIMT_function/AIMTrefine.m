function [diamd,diamu,seld,selu]=AIMTrefine(xim,xp,yp,dir,mdiam,nker,dbf);

if(dbf)
    imagesc(xim)
    colormap(gray)
    hold on
    plot(xp,yp,'y');
end;

rllim=0.25;
rulim=2;
drlim=rulim-rllim;
thlev=0.3;

mrad=0.5*mdiam;
lseg=round(drlim*mrad);

for ct=1:length(xp),
    min=0;
    for ctp=-nker:nker,
        xstart=xp(ct)+ctp*cos(dir(ct));
        ystart=yp(ct)+ctp*sin(dir(ct));
        
        lx=xstart+mrad*cos(dir(ct)+pi/2)*[rllim,rulim];
        ly=ystart+mrad*sin(dir(ct)+pi/2)*[rllim,rulim];
        
        p(ctp+nker+1,:)=improfile(xim,lx ,ly, lseg);
        
        %lx=xstart+mrad*cos(dir(ct)+pi/2)*[-rllim,rllim];
        %ly=ystart+mrad*sin(dir(ct)+pi/2)*[-rllim,rllim];
        %pin(ctp+nker+1,:)=improfile(xim,lx ,ly,round(mrad));
    end;
    
    pf=mean(p);
    try
        [diamu(ct),selu(ct)]=AIMTprofile(pf,thlev,rllim*mrad,mrad,0.02, 0.2, dbf);
    catch
        keyboard
    end;
    diamu(ct)=rllim*mrad+diamu(ct);
    
    %if(selu(ct)>0)
    %    if(dbf)
    %        plot(xp(ct)+diamu(ct)*cos(dir(ct)+pi/2),yp(ct)+diamu(ct)*sin(dir(ct)+pi/2),'or');
    %    end;
    %end;
end;
md=median(diamu);
vd=abs(diamu-md)/md;
diamu(find(vd>0.2))=md;
if(dbf)
    ns=find(selu==1);
    plot(xp+diamu.*cos(dir+pi/2),yp+diamu.*sin(dir+pi/2),'r');
    plot(xp(ns)+diamu(ns).*cos(dir(ns)+pi/2),yp(ns)+diamu(ns).*sin(dir(ns)+pi/2),'or');
end;


for ct=1:length(xp),
    for ctp=-nker:nker,
        xstart=xp(ct)+ctp*cos(dir(ct));
        ystart=yp(ct)+ctp*sin(dir(ct));
        
        lx=xstart+mrad*cos(dir(ct)-pi/2)*[rllim,rulim];
        ly=ystart+mrad*sin(dir(ct)-pi/2)*[rllim,rulim];
        
        p(ctp+nker+1,:)=improfile(xim,lx ,ly, lseg);
    end;
    
    pf=mean(p);
    [diamd(ct),seld(ct)]=AIMTprofile(pf,thlev,rllim*mrad,mrad,0.02, 0.2, dbf);
    diamd(ct)=rllim*mrad+diamd(ct);
    
    %if(seld(ct)>0)
    %    if(dbf)
    %        plot(xp(ct)+diamd(ct)*cos(dir(ct)-pi/2),yp(ct)+diamd(ct)*sin(dir(ct)-pi/2),'or');
    %    end;
    %end;
end;

md=median(diamd);
vd=abs(diamd-md)/md;
diamd(find(vd>0.2))=md;
if(dbf)
    ns=find(seld==1);
    plot(xp+diamd.*cos(dir-pi/2),yp+diamd.*sin(dir-pi/2),'r');
    plot(xp(ns)+diamd(ns).*cos(dir(ns)-pi/2),yp(ns)+diamd(ns).*sin(dir(ns)-pi/2),'or');
end;
