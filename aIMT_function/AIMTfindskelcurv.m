function [xp,yp,dist, ddist] = AIMTfindskelcurv(xbw, smoothxy, mindist, remperc, dbf)

%% init



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the skeleton and prun it

% thinning
xbwsk = bwmorph(xbw,'thin',inf);

xbwdist=(bwdist(1-xbw,'Quasi-Euclidean'));
xbwdist4 = (xbwdist>0).*(xbwdist<mindist);

% spur those ends distant up to 4 pixel from border
xbwends = xbwsk-bwmorph(xbwsk,'spur',1);
xbwdiff = xbwends.*xbwdist4;
ec = any(any(xbwdiff));
ct = 0;
while ec
    ct = ct+1; % max 5 times
    xbwsk = xbwsk.*(1-xbwdiff);
    xbwends = xbwsk.*(1-bwmorph(xbwsk,'spur',1));
    xbwdiff = xbwends.*xbwdist4;
    ec = any(any(xbwdiff)) && (ct<4);
end

xbwsk0 = xbwsk;
xbwprun = xbwsk;

redo=1;
while redo

    % find ends
    endslab = [];
    xend = []; yend = [];

    xbwsksp = bwmorph(xbwsk,'spur',1);
    xbwends = xbwsk-xbwsksp;

    endslab = bwlabel(xbwends);
    yend=[];
    xend=[];
    for ct=1:max(max(endslab))
        %% EG 2007-09-14
        %% Sometimes a two-pixel endpoint results
        %% We save both points
        [ytmp, xtmp]= find(endslab==ct);
        yend=[yend;ytmp];
        xend=[xend;xtmp];
        %[yend(ct), xend(ct)] = find(endslab==ct);
    end

    % find nodes
    xnode = [];  % node's coordinates
    ynode = [];  % node's coordinates
    nadj = [];   % node's number of adiacency
    idnode = []; % primary key of the node
    N = 0;       % number of nodes

    xbwskd = double(xbwsk);
    [rig col] = find(xbwsk);
    for ct = 1:length(rig)
        i = rig(ct);
        j = col(ct);
        numadj = xbwskd(i-1,j) + xbwskd(i,j-1) + xbwskd(i+1,j) + xbwskd(i,j+1) ;
        d1 = double( xbwsk(i-1,j-1) & ( ~( xbwsk(i-1,  j) | xbwsk(  i,j-1))));
        d2 = double( xbwsk(i+1,j-1) & ( ~( xbwsk(  i,j-1) | xbwsk(i+1,  j))));
        d3 = double( xbwsk(i+1,j+1) & ( ~( xbwsk(i+1,  j) | xbwsk(  i,j+1))));
        d4 = double( xbwsk(i-1,j+1) & ( ~( xbwsk(i-1,  j) | xbwsk(  i,j+1))));
        numadj = numadj  + d1 + d2 + d3 + d4;
        if (numadj >= 3)
            N = N+1;
            xnode(N) = j;
            ynode(N) = i;
            nadj(N) = numadj;
            idnode(N) = N;
        end
    end

    if N==0
        break;
    end

    % define structure of arches
    % xs, ya: coordinates of arch
    % l: lenght
    % idnode: ID of related node
    % archnode: archnode(i) is the ID of the node related to arch i
    arch = struct('xa',[],'yx',[],'l',[],'idnode',[]);
    archnode = [];

    % find arches: from every end up to a node
    for ct=1:length(xend)
        ec=0;
        xcurr = xend(ct);
        ycurr = yend(ct);
        l=1;
        arch(ct).xa(l) = xcurr;
        arch(ct).ya(l) = ycurr;
        xbwcurr = xbwsk;
        while ~ec
            % find next point...
            l=l+1;
            xbwcurr(ycurr,xcurr) = 0;
            [yup, xup] = find(xbwcurr(ycurr-1:ycurr+1,xcurr-1:xcurr+1));
            % if next is a node
            if length(yup)>1
                [yup, xup] = find(xbwcurr(ycurr-1:ycurr+1,xcurr-1:xcurr+1).*([0,1,0;1,0,1;0,1,0]));
            end
            xcurr = xcurr-2+xup;
            ycurr = ycurr-2+yup;
            arch(ct).xa(l) = xcurr;
            arch(ct).ya(l) = ycurr;
            % ...until a node is found
            indN = find(all(eq([xcurr; ycurr]*ones(1,length(xnode)),[xnode;ynode])));
            if ~isempty(indN)
                ec = 1;
                % save info about length and node
                arch(ct).l = l;
                arch(ct).idnode = indN;
                archnode(ct) = indN;
            end
        end
    end

    % count number of arch for every node
    narchnode = zeros(1,N);
    for ct=1:numel(arch)
        narchnode(arch(ct).idnode) = narchnode(arch(ct).idnode)+1;
    end

    % prun the skeleton
    % find nodes with number of arch>0
    indprun = find(narchnode);
    for ct=1:length(indprun)
        % find related arches
        indarch = find(archnode==indprun(ct));
        % sort them by length
        [larch indpra] = sort([arch(indarch).l]);
        ct2=1;
        if numel(indpra)~=0
            while (nadj(indprun(ct))>2 && narchnode(indprun(ct))>0) || (larch(ct2)<5 && nadj(indprun(ct))>1)
            %while (nadj(indprun(ct))>1 && narchnode(indprun(ct))>0) || (larch(ct2)<5 && nadj(indprun(ct))>1)
                xprun = arch(indarch(indpra(ct2))).xa;
                yprun = arch(indarch(indpra(ct2))).ya;
                for ct3=1:length(xprun)
                    xbwprun(yprun(ct3),xprun(ct3)) = 0;
                end
                % update
                archnode(ct) = 0;
                nadj(indprun(ct)) = nadj(indprun(ct))-1;
                narchnode(indprun(ct)) = narchnode(indprun(ct))-1;
                ct2 = ct2+1;
                if ct2>length(indpra), break; end
            end
        end
    end

    % redraw nodes
    for ct=1:length(xnode)
        xbwprun(ynode(ct),xnode(ct))=1;
    end
    % adjust m-adiacency
    xbwprun = bwmorph(xbwprun,'thin',1);

    % delete arch
    arch(find(archnode)==0)=[];

    %if dbf, figure();imagesc(double(xbw)+double(xbwsk)+2*double(xbwprun)); end

    % update image
    xbwsk = xbwprun;

    % update end condition
    if isempty(find(nadj>2, 1))
        redo=0;
    end
end

%xbwsk = bwmorph(xbwsk, 'spur', 30);

if dbf, figure();imagesc(double(xbw)+double(xbwsk)+2*double(xbwprun)); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find spline on skeleton

% find ends
xend = []; yend = [];
xbwsksp = bwmorph(xbwsk,'spur',1);
xbwends = xbwsk-xbwsksp;
endslab = bwlabel(xbwends);
if(max(max(endslab))==0)
    [rc,cc]=find(xbwsksp);
    [r,c]=find(xbw-imerode(xbw,ones(3)));
    nr=length(r);
    dmat=(r*ones(1,nr)-ones(nr,1)*r').^2+(c*ones(1,nr)-ones(nr,1)*c').^2;
    [rm,cm]=find(dmat==max(dmat(:)));
    yend=r(rm);
    xend=c(rm);
    for ct=1:length(yend),
        [cx,cy,c]=improfile(xbwsk,[cc(1),xend(ct)],[rc(1),yend(ct)]);
        for ct2=1:length(cx),
            xbwsksp(round(cx(ct2)),round(cy(ct2)))=1;
        end;
    end;
    xbwsk=bwmorph(imdilate(xbwsksp,ones(3)),'thin','Inf');
    xend = []; yend = [];
    xbwsksp = bwmorph(xbwsk,'spur',1);
    xbwends = xbwsk-xbwsksp;
    endslab = bwlabel(xbwends);
end

for ct=1:max(max(endslab))
    [yend(ct), xend(ct)] = find(endslab==ct);
end



% start from the first end
try
    xcurr = xend(1);
    ycurr = yend(1);
catch
    keyboard
end;
xskel = xcurr; yskel = ycurr;
l=1; ec=0;
xbwcurr = xbwsksp;
while ~ec
    % find next point...
    l=l+1;
    xbwcurr(ycurr,xcurr) = 0;
    [yup, xup] = find(xbwcurr(ycurr-1:ycurr+1,xcurr-1:xcurr+1));
    xcurr = xcurr-2+xup;
    ycurr = ycurr-2+yup;
    dist = sqrt((xcurr-xskel(end)*ones(length(xcurr),1)).^2+(ycurr-yskel(end)*ones(length(xcurr),1)).^2);
    ind = find(dist==min(dist));
    xcurr = xcurr(ind);
    ycurr = ycurr(ind);
    xskel = [xskel xcurr];
    yskel = [yskel ycurr];
    % ...until the other end is found
    try
        if isempty(xcurr) | any(all(eq([xcurr; ycurr]*ones(size(xend)),[xend; yend])));
            ec = 1;
        end
    catch
        keyboard
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sample the skeleton
l = length(xskel);
switch l
    case {1,2,3}
        xsp = xskel;
        ysp = yskel;

        % recalculating points
        xre1up = xsp(1); yre1up = ysp(1);
        xre2up = xsp(end); yre2up = ysp(end);
        xre1dw = xsp(end); yre1dw = ysp(end);
        xre2dw = xsp(1); yre2dw = ysp(1);

    case {4,5,6,7,8}
        xsp = xskel;
        ysp = yskel;

        % recalculating points
        xre1up = xsp(2); yre1up = ysp(2);
        xre2up = xsp(end-1); yre2up = ysp(end-1);
        xre1dw = xsp(end-1); yre1dw = ysp(end-1);
        xre2dw = xsp(2); yre2dw = ysp(2);

    case {9,10,11,12}
        xsp = xskel(1:1:length(xskel));
        ysp = yskel(1:1:length(yskel));

        % recalculating points
        xre1up = xsp(2); yre1up = ysp(2);
        xre2up = xsp(4); yre2up = ysp(4);
        xre1dw = xsp(end-1); yre1dw = ysp(end-1);
        xre2dw = xsp(end-3); yre2dw = ysp(end-3);

    case {13,14,15,16,17,18,19,20}
        xsp = xskel(1:2:length(xskel));
        ysp = yskel(1:2:length(yskel));

        % recalculating points
        xre1up = xsp(2); yre1up = ysp(2);
        xre2up = xsp(3); yre2up = ysp(3);
        xre1dw = xsp(end-1); yre1dw = ysp(end-1);
        xre2dw = xsp(end-2); yre2dw = ysp(end-2);

%     otherwise
%         xsp = xskel(1:3:length(xskel));
%         ysp = yskel(1:3:length(yskel));
% 
%         % recalculating points
%         xre1up = xsp(3); yre1up = ysp(3);
%         xre2up = xsp(4); yre2up = ysp(4);
%         xre1dw = xsp(end-2); yre1dw = ysp(end-2);
%         xre2dw = xsp(end-3); yre2dw = ysp(end-3);
        
	otherwise
        xsp = xskel(1:3:length(xskel));
        ysp = yskel(1:3:length(yskel));

        onset=4;
        offset=4;
        % recalculating points
        xre1up = xsp(onset); yre1up = ysp(onset);
        xre2up = xsp(onset+1+offset); yre2up = ysp(onset+1+offset);
        xre1dw = xsp(end-onset-1); yre1dw = ysp(end-onset-1);
        xre2dw = xsp(end-onset-offset); yre2dw = ysp(end-onset-offset);
end

if dbf, axis equal; hold on, plot(xsp,ysp,'+g'); end

%% final points recalculation (up)
step = 2;
xre = xsp(1);
yre = ysp(1);

updateup = [xre1up-xre2up yre1up-yre2up];
normaup = norm(updateup);
if normaup==0, normaup = 1; end
updateup = updateup/normaup;

while ((xre>1 && yre>1) && (yre<size(xbw,1) && xre<size(xbw,2)) && ...
        all(xbw(yre-1:yre+1,xre)==1) && all(xbw(yre,xre-1:xre+1)==1))
    xsp(1) = xre;
    ysp(1) = yre;
    xre = fix(xre+step*updateup(1));
    yre = fix(yre+step*updateup(2));
    step=step;
end

% final points recalculation (down)
step = 2;
xre = xsp(end);
yre = ysp(end);

updatedw = [xre1dw-xre2dw yre1dw-yre2dw];
normadw = norm(updatedw);
if normadw==0, normadw = 1; end
updatedw = updatedw/normadw;

while ((xre>1 && yre>1) && (yre<size(xbw,1) && xre<size(xbw,2)) && ...
        all(xbw(yre-1:yre+1,xre)==1) && all(xbw(yre,xre-1:xre+1)==1))
    xsp(end) = xre;
    ysp(end) = yre;
    xre = fix(xre+step*updatedw(1));
    yre = fix(yre+step*updatedw(2));
    step = step;
end

%% Remove starting and ending parts
lsp=sum(sqrt(diff(xsp).^2+diff(ysp).^2));
dsp=cumsum(sqrt(diff(xsp).^2+diff(ysp).^2));
cutl=remperc*lsp;

n1=find(dsp>cutl);
dsp(1:n1(1))=[];
xsp(1:n1(1))=[];
ysp(1:n1(1))=[];

n2=find(dsp>(lsp-cutl));
dsp(n2(1):end)=[];
xsp(n2(1):end)=[];
ysp(n2(1):end)=[];

if dbf, plot(xsp(1),ysp(1),'+y'); plot(xsp(end),ysp(end),'+y'); end

%% Curvilinear coordinate
param(1,1) = 0;
for ct=2:length(xsp),
    dist = sqrt((xsp(ct)-xsp(ct-1))^2+(ysp(ct)-ysp(ct-1))^2);
    param(ct,1)=param(ct-1,1)+dist;
end;

%% Spline interpolation
try
    ppx = csaps(param, xsp, smoothxy);
    ppy = csaps(param, ysp, smoothxy);
    xp = fnval(ppx, ppx.breaks(1):step:ppx.breaks(end));
    yp = fnval(ppy, ppy.breaks(1):step:ppy.breaks(end));
catch
    keyboard;
end;

if dbf, hold on, plot(xp,yp,'r'); end

%% find axis' curvature

k = CROMOcurv(ppx,ppy,2,0);

%% find contour function

derivx=fnder(ppx);
derivy=fnder(ppy);
dx=fnval(derivx,ppx.breaks(1):step:ppx.breaks(end));
dy=fnval(derivy,ppy.breaks(1):step:ppy.breaks(end));
theta=atan2(dy,dx);
distup=[];
distdown=[];

% Steps along the axis
for ct=1:1:length(xp)

    %Initialize variables of the up cycle
    stepup = 0.1;
    colretup = fix(xp(ct));
    rigretup = fix(yp(ct));
    xrup=xp(ct);
    yrup=yp(ct);
    %[rigretup colretup] = CROMOcart2ind(fix(xrup),fix(yrup),size(xbw,1));
    distanceup = 0;

    % Steps going up
    while ((rigretup>0 && colretup>0) && (rigretup<=size(xbw,1) && colretup<=size(xbw,2))) && (xbw(rigretup,colretup)==1)
        xretup = xrup;
        yretup = yrup;
        distanceup = sqrt((xretup-xp(ct))^2+(yretup-yp(ct))^2);
        xrup = xp(ct)+stepup*cos(theta(ct)+pi/2);
        yrup = yp(ct)+stepup*sin(theta(ct)+pi/2);
        colretup = fix(xrup);
        rigretup = fix(yrup);
        %[rigretup colretup] = CROMOcart2ind(fix(xrup),fix(yrup),size(xbw,1));
        stepup = stepup+0.1;
    end
    distup = [distup,distanceup];

    % Initialize variables of the down cycle
    stepdown = 0.1;
    xrdown=xp(ct);
    yrdown=yp(ct);
    colretdown = fix(xp(ct));
    rigretdown = fix(yp(ct));
    %[rigretdown colretdown] = CROMOcart2ind(fix(xrdown),fix(yrdown),size(xbw,1));
    distancedown = 0;

    % Steps going down
    while ((rigretdown>0 && colretdown>0) && (rigretdown<=size(xbw,1) && colretdown<=size(xbw,2))) && (xbw(rigretdown,colretdown)==1)
        xretdown=xrdown;
        yretdown=yrdown;
        distancedown=sqrt((xretdown-xp(ct))^2+(yretdown-yp(ct))^2);
        xrdown=xp(ct)+stepdown*cos(theta(ct)-pi/2);
        yrdown=yp(ct)+stepdown*sin(theta(ct)-pi/2);
        %[rigretdown colretdown]=CROMOcart2ind(fix(xrdown),fix(yrdown),size(xbw,1));
        colretdown = fix(xrdown);
        rigretdown = fix(yrdown);
        stepdown=stepdown+0.1;
    end
    distdown=[distdown,distancedown];

end

%Find contour function and differential contour function
dist = distup+distdown;
ddist = (abs(distup-distdown));
ddistq = ddist.^2;

%% Find straits

smcc = [dist, dist]; %smooth(dist,5);

[ccmin indccmin] = min(smcc);

dersmcc = smcc;


%% debug info

if (dbf)
    iptsetpref('Imshowborder','tight')
    figure,
    
    subplot(3,1,1)
    plot(ppx.breaks(1):step:ppx.breaks(end), distup,'-g',ppx.breaks(1):step:ppx.breaks(end), -distdown,'-b')
    title('Chromosome border')
    
    subplot(3,1,2)
    plot(ppx.breaks(1):step:ppx.breaks(end), ddist,'-k')
    hold on;
    %plot(ppx.breaks(1):step:ppx.breaks(end), ddistq,'-b')
    %hold off
    title('difference between the two contour')

    subplot(3,1,3)
    plot(ppx.breaks(1):step:ppx.breaks(end), dist,'-c')
    title('Contour function')

    %figure()
    %plot(dersmcc,'-b');
    %hold on
    %plot(dist,'-k');
    %title('derivate of the contour function')
end