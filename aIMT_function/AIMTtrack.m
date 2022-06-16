function [x,y,diam,dir]=AIMTtrack(xim,xstart,ystart,dstart,maxfail,maxiter,maxd,step,dbf)

x(1)=xstart;
y(1)=ystart;
dir(1,:)=dstart;

sed=round(max(1,maxd*0.05));

if(dbf)
    figure;
    imagesc(xim);
    colormap(gray);
    hold on;
end;

failcount=0;
ct=1;
while(ct<maxiter & failcount<maxfail)
    
    if(ct>1)
        mdiam=mean(diam(1:end-1));
        mdir=mean(dir);
        mind=max(1, size(dir,1)-5);
        vd=mean(dir(mind:end,:));
        vd=vd/norm(vd);
    else
        %mdiam=diam(1);
        vd=dir(ct,:);
        vd=vd/norm(vd);
    end;
    
    xorig=x(ct);
    yorig=y(ct);
    vn=[vd(2),-vd(1)];
    
    lx=x(ct)+maxd*[vn(2),-vn(2)];
    ly=y(ct)+maxd*[-vn(1),vn(1)];
    lseg=sqrt(diff(lx).^2+diff(ly).^2);
    
    p=improfile(xim,lx ,ly ,lseg);
    p(isnan(p))=0;
    
    c=fix(lseg/2);
    [center, U, obj_fcn] = fcm_my(p(1:c), 2);
    th1=max(mean(center),0.1);
    
    [center, U, obj_fcn] = fcm_my(p(c+1:end), 2);
    th2=max(mean(center),0.1);
    
    try
        fl=[imopen(p(1:c)>th1,ones(sed,1));imopen(p(c+1:end)>th2,ones(sed,1))];
    catch
        keyboard
    end;
    nl=find(fl(1:c)==1);
    nr=c+find(fl(c+1:end)==1);
    
    if(isempty(nl) | isempty(nr))
        failcount=failcount+1;
        dir(ct+1,:)=vd;
        x(ct+1)=x(ct)+step*vd(2);
        y(ct+1)=y(ct)+step*vd(1);
        if(ct==1)
            diam(ct)=-1;
        end;
        diam(ct+1)=diam(ct);
    else
        nl=nl(end);
        nr=nr(1);
        n_new=0.5*(nl+nr);
        db=(c-n_new);
        diam(ct)=nr-nl;
        
        if(ct>1 && (abs(mdiam-diam(ct))/mdiam>0.2 | db>mdiam*0.2))
            failcount=failcount+1;
            diam(ct)=mdiam;
            dir(ct+1,:)=vd;
            x(ct+1)=x(ct)+step*vd(2);
            y(ct+1)=y(ct)+step*vd(1);
            diam(ct+1)=mdiam;
        else
            failcount=0;
            x(ct)=x(ct)-db*vn(2);
            y(ct)=y(ct)-db*vn(1);
            
            if(ct==1)
                Dn=vd;
            else
                D=[y(ct)-y(ct-1),x(ct)-x(ct-1)];
                Dn=D/sqrt(sum(D.^2));
                
            end;
            
            dir(ct+1,:)=Dn;
            x(ct+1)=x(ct)+step*Dn(2);
            y(ct+1)=y(ct)+step*Dn(1);
            diam(ct+1)=diam(ct);
        end;
    end;
    
    if(dbf)
        plot(xorig,yorig,'oy');
        if(failcount)
            plot(x(ct),y(ct),'*m');
            plot(x(ct)-0.5*diam(ct)*vn(2),y(ct)-0.5*diam(ct)*vn(1),'*m');
            plot(x(ct)+0.5*diam(ct)*vn(2),y(ct)+0.5*diam(ct)*vn(1),'*m');
        else
            plot(x(ct),y(ct),'*r');
            plot(x(ct)-0.5*diam(ct)*vn(2),y(ct)-0.5*diam(ct)*vn(1),'*r');
            plot(x(ct)+0.5*diam(ct)*vn(2),y(ct)+0.5*diam(ct)*vn(1),'*r');
        end;
    end;
    
    ct=ct+1;
end;