function [aIMT_vec]=aIMT_thickness2(xim,nframe,xp,yp,diam,ncol,tol_aIMT,dbf)
% xim = immagine gia passata al BWroi e senza freccia/cursore
% col_ok = colonne in cui ha senso cercare aIMT
% AC = vettore coordinate bordo aorta
% tol_aIMT = numero di pixel massimi per uno spessore accettabile. Dipende
% se l'immagine è zoomata o no

%N=30;
%NN=60;
tol=0.2; %quanto l'ext dev'essere più chiaro dell'int per considerare il punto
passo = 1;

fdiam=0.5;
sdiam=1-fdiam;
% tol_aIMT=10; %numero di pixel massimi per uno spessore accettabile

aIMT_vec=zeros(length(ncol),4);
dx=diff(xp);
dx=[dx(1),dx];
dy=diff(yp);
dy=[dy(1),dy];
nd=sqrt(dx.^2+dy.^2);
dx=dx./nd;
dy=dy./nd;

if dbf
    figure(nframe)
    imagesc(xim); colormap gray
    hold on
    plot(xp,yp,'y')
end
for ctcol=1:passo:length(ncol)
    
    nsel=ncol(ctcol);
    %tmpcol=abs(xp-col_ok(ctcol));
    %nsel=find(tmpcol==min(tmpcol));
    
    N=fdiam*diam(nsel);
    NN=2*N;
    
    rowUp=xp(nsel)-dy(nsel)*0.5*diam(nsel);
    colUp=yp(nsel)+dx(nsel)*0.5*diam(nsel);
    rowLow=xp(nsel)+dy(nsel)*0.5*diam(nsel);
    colLow=yp(nsel)-dx(nsel)*0.5*diam(nsel);
    
    xlim=fix([rowUp+dy(nsel)*N,rowUp-dy(nsel)*N]);
    ylim=fix([colUp-dx(nsel)*N,colUp+dx(nsel)*N]);
    [xup,yup,NVPup]=improfile(xim,xlim,ylim,2*NN);
    if dbf
        figure(nframe)
        %plot(xlim,ylim,'-og')
    end
    
    xlim=fix([rowLow+dy(nsel)*N,rowLow-dy(nsel)*N]);
    ylim=fix([colLow-dx(nsel)*N,colLow+dx(nsel)*N]);
    [xlow,ylow,NVPlow]=improfile(xim,xlim,ylim,2*NN);
    
    if dbf
        figure(nframe)
        %plot(xlim,ylim,'-og')
    end
    
    %NVPup=xim(rowUp-fix(N/2):rowUp+fix(N/2),col_ok(ctcol));
    %NVPlow=xim(rowLow-fix(N/2):rowLow+fix(N/2),col_ok(ctcol));
    %controllo se il punto è affidabile verificando che ext>>int
    imax=max(xim(:,fix(xp(nsel))));
    %[lim,peak]=aIMTfindpeak(NVPup,tol,imax,[],dbf);
    %[lim,peak]=aIMTfindmix(NVPup,3,mod(ctcol-1,10)==0);
    [lim,peak]=aIMTfindmix(NVPup,3,0);
    
    if (~isempty(peak))
        aIM_Up1=[xup(lim),yup(lim)];
        aIM_Up2=[xup(peak),yup(peak)];
        aIMT1=sqrt(sum((aIM_Up2-[rowUp,colUp]).^2));
        aIMT2=sqrt(sum((aIM_Up2-aIM_Up1).^2));
        if aIMT1<=tol_aIMT
            aIMT_vec(ctcol,1)=aIMT1;
            aIMT_vec(ctcol,3)=aIMT2;
            if dbf
                figure(nframe)
                hold on
                %plot(aIM_Up2(1),aIM_Up2(2),'r.')
                %plot(xup(lim),yup(lim),'.y')
                plot([aIM_Up2(1),xup(lim)],[aIM_Up2(2),yup(lim)],'o-w')
            end
        end
        
    end
    
    imax=max(xim(:,fix(xp(nsel))));
    %[lim,peak]=aIMTfindpeak(flipud(NVPlow),tol,imax,[],dbf);
    %[lim,peak]=aIMTfindmix(flipud(NVPlow),3,mod(ctcol-1,10)==0);
    [lim,peak]=aIMTfindmix(flipud(NVPlow),3,0);
    
    if ~isempty(peak)
        peak=length(NVPlow)-peak+1;
        lim=length(NVPlow)-lim;
        
        aIM_Low1=[xlow(lim),ylow(lim)];
        aIM_Low2=[xlow(peak),ylow(peak)];
        aIMT1=sqrt(sum((aIM_Low2-[rowLow,colLow]).^2));
        aIMT2=sqrt(sum((aIM_Low2-aIM_Low1).^2));
        if aIMT1<=tol_aIMT
            aIMT_vec(ctcol,2)=aIMT1;
            aIMT_vec(ctcol,4)=aIMT2;
            if dbf
                figure(nframe)
                hold on
                %plot(aIM_Low2(1),aIM_Low2(2),'r.')
                %plot(xlow(lim),ylow(lim),'.y')
                plot([aIM_Low2(1),xlow(lim)],[aIM_Low2(2),ylow(lim)],'o-w')
            end
        end
    end
end
