function res=AIMTthick(xim, xseg, vert, smooth, parthick, nframe, dbf);

[xp,yp,dist, ddist] = AIMTfindskelcurv(xseg,smooth,20,0.15,dbf);
try
    [tratto_aortax,distrid]=AIMT_findStrip2(xim,xp,yp,dist,[],1,dbf);
catch
    keyboard
end;

if ~isempty(tratto_aortax)
    %AIMT_vec=AIMT_thickness(xim, tratto_aortax, vert, parthick,dbf);
    aIMT_vec=aIMT_thickness2(xim,nframe,xp,yp,dist,tratto_aortax,20,dbf);
else
    aIMT_vec=-1;
end;

res=struct('nframe',nframe,'ACborder',vert,'aIMT_col',tratto_aortax,'xydist',[xp;yp;dist],'AIMT',aIMT_vec);