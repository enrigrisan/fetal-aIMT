function [aIMT_vec]=aIMT_thickness(xim,nframe,xp,yp,diam,col_ok,tol_aIMT,dbf)
% xim = immagine gia passata al BWroi e senza freccia/cursore
% col_ok = colonne in cui ha senso cercare aIMT
% AC = vettore coordinate bordo aorta
% tol_aIMT = numero di pixel massimi per uno spessore accettabile. Dipende
% se l'immagine è zoomata o no

%N=30;
%NN=60;
tol=0.2; %quanto l'ext dev'essere più chiaro dell'int per considerare il punto
passo = 1;
% tol_aIMT=10; %numero di pixel massimi per uno spessore accettabile

aIMT_vec=zeros(length(col_ok),2);
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
for ctcol=1:passo:length(col_ok)
    
    tmpcol=abs(xp-col_ok(ctcol));
    nsel=find(tmpcol==min(tmpcol));
    
    N=0.5*diam(nsel);
    NN=2*N;
    
    rowUp=xp(nsel)-dy(nsel)*0.5*diam(nsel);
    colUp=yp(nsel)+dx(nsel)*0.5*diam(nsel);
    rowLow=xp(nsel)+dy(nsel)*0.5*diam(nsel);
    colLow=yp(nsel)-dx(nsel)*0.5*diam(nsel);
    
    xlim=fix([rowUp+dy(nsel)*N/2,rowUp-dy(nsel)*N/2]);
    ylim=fix([colUp-dx(nsel)*N/2,colUp+dx(nsel)*N/2]);
    [xup,yup,NVPup]=improfile(xim,xlim,ylim,2*NN);
    if dbf
        figure(nframe)
        %plot(xlim,ylim,'-og')
    end
    
    xlim=fix([rowLow+dy(nsel)*N/2,rowLow-dy(nsel)*N/2]);
    ylim=fix([colLow-dx(nsel)*N/2,colLow+dx(nsel)*N/2]);
    [xlow,ylow,NVPlow]=improfile(xim,xlim,ylim,2*NN);
    
    if dbf
        figure(nframe)
        %plot(xlim,ylim,'-og')
    end
    
    %NVPup=xim(rowUp-fix(N/2):rowUp+fix(N/2),col_ok(ctcol));
    %NVPlow=xim(rowLow-fix(N/2):rowLow+fix(N/2),col_ok(ctcol));
    %controllo se il punto è affidabile verificando che ext>>int
    extUp=mean(NVPup(end-fix(N/2)+1:end));
    intUp=mean(NVPup(1:fix(N/2)));
    extLow=mean(NVPlow(1:fix(N/2)));
    intLow=mean(NVPlow(end-fix(N/2)+1:end));
    if intUp<tol*extUp %punto superiore affidabile
        peaksn=fpeak([1:length(NVPup)],NVPup,10,[1,length(NVPup),0,max(xim(:,fix(xp(nsel))))]);
        if (~isempty(peaksn) && peaksn(1,1)==1)
            peaksn=peaksn(2:end,:);
        end
        nlim=find(NVPup<mean(NVPup));
        peaksn(peaksn(:,1)<nlim(end),:)=[];
        if (~isempty(peaksn) && any(peaksn(:,2)<0.05)) % tolgo i minimi sotto una certa soglia
            to_del=find(peaksn(:,2)<0.05);
            for k=1:length(to_del)
                peaksn(to_del(k),:)=[0 0];
            end
            peaksntmp1=nonzeros(peaksn(:,1));
            peaksntmp2=nonzeros(peaksn(:,2));
            clear peaksn
            peaksn=[peaksntmp1';peaksntmp2']';
            clear peaksntmp1 peaksntmp2
        end
        if (~isempty(peaksn))
            %             figure; plot(xim(rowUp-fix(N/2)-NN:rowUp+fix(N/2)+NN,col_ok(ctcol))); title(num2str(ctcol))
            %             hold on; plot(peaksn(:,1),peaksn(:,2),'ko')
            %             hold on; plot(peaksn(1,1),peaksn(1,2),'ro')
            %aIM_Up2=peaksn(1,1)+(rowUp-fix(N/2)-NN-1);   %indice del bordo esterno dell'aIMT superiore
            %aIMT_Up=-peaksn(1,1)+fix(N/2)+NN+1;
            aIM_Up2=[xup(peaksn(1,1)),yup(peaksn(1,1))];
            aIMT_Up=sqrt(sum((aIM_Up2-[rowUp,colUp]).^2));
            nlim=find(NVPup>NVPup(peaksn(1,1)));
            %aIMT_Up=length(xup)-peaksn(1,1);
            if aIMT_Up<=tol_aIMT
                aIMT_vec(ctcol,:)=aIMT_Up;
                if dbf
                    figure(nframe)
                    hold on
                    plot(aIM_Up2(1),aIM_Up2(2),'r.')
                    plot(xup(nlim(1)),yup(nlim(1)),'y.')
                end
            end
            
        end
    end
    if intLow<tol*extLow %punto inferiore affidabile
        peakdx=fpeak([1:length(NVPlow)],NVPlow,10,[1,length(NVPlow),0,max(xim(:,col_ok(ctcol)))]);
        if (~isempty(peakdx) && peakdx(end,1)==length([1:N+2*NN]))
            peakdx=peakdx(1:end-1,:);
        end
        
        nlim=find(NVPlow>mean(NVPlow));
        peakdx(peakdx(:,1)>nlim(end),:)=[];
        if (~isempty(peakdx) && any(peakdx(:,2)<0.05)) % tolgo i minimi sotto una certa soglia
            to_del=find(peakdx(:,2)<0.05);
            for k=1:length(to_del)
                peakdx(to_del(k),:)=[0 0];
            end
            peakdxtmp1=nonzeros(peakdx(:,1));
            peakdxtmp2=nonzeros(peakdx(:,2));
            clear peakdx
            peakdx=[peakdxtmp1';peakdxtmp2']';
            clear peakdxtmp1 peakdxtmp2
        end
        
        if ~isempty(peakdx)
            % pongo condizioni sullo spessore, per scartare quelli
            % implausibili
            %aIM_Low2=peakdx(1,1)+(rowLow-fix(N/2)-NN-1);   %indice del bordo esterno dell'aIMT inferiore
            %aIMT_Low=peakdx(1,1)+(-fix(N/2)-NN-1);
            aIM_Low2=[xlow(peakdx(end,1)),ylow(peakdx(end,1))];
            aIMT_Low=sqrt(sum((aIM_Low2-[rowLow,colLow]).^2));
            nlim=find(NVPlow>NVPlow(peakdx(end,1)));
            if aIMT_Low<=tol_aIMT
                aIMT_vec(ctcol,2)=aIMT_Low;
                if dbf
                    figure(nframe)
                    hold on
                    plot(aIM_Low2(1),aIM_Low2(2),'r.')
                    plot(xlow(nlim(end)),ylow(nlim(end)),'y.')
                end
            end
        end
    end
    
end