function callbackfn_rawovplot(src,edat,d,varargin)
% ** function callbackfn_rawovplot(src,edat,d,varargin)
% is part of the spike detection gui (threshdetgui.m).
% It is the callback of ButtonDownFcn of a 'dummy' patch (a rectangle) in
% subplot sp.rawOv.axH. Extracts an excerpt of data d centered on the x
% coordinate of a mouse click in the subplot and plots it in subplot
% sp.rawExc.axH. If input argument 'x' is explicitly specified the
% coordinate of the mouse click is ignored

global wp ap sp ds

x=[];
pvpmod(varargin);
etslconst;

% -----  preliminaries
% make sure excerpt is leq total # points
wp.excLen_pts=min(wp.rawNRow,wp.excLen_pts);
% ----- do it
% catch inadvertent clicks on empty subplot
if ~isempty(d)
  if isempty(x)
    % get x coordinate of last mouse click and store it
    wp.curPtX=get(sp.rawOv.axH,'currentpoint');
    wp.curPtX=wp.curPtX(1,1);
  else
    wp.curPtX=x;
  end
  % transform x axis coordinate to indices into d
  ix=round(wp.curPtX/diff(get(sp.rawOv.axH,'xlim'))*wp.rawNRow + round(wp.excLen_pts*[-.5 .5]));
  % curb and/or shift indices to avoid negative or too large values
  if diff(ix)+1>wp.rawNRow
    % in this case, chosen excerpt length is larger than total raw data
    % trace
    ix=[1 wp.rawNRow];
  else
    if ix(1)<1
      ix=ix-ix(1)+1;
    end
    if ix(2)>wp.rawNRow
      ix=ix+wp.rawNRow-ix(2);
    end
  end
  % store ix as field in wp
  wp.excXLim_pts=ix;
  subplot(sp.rawExc.axH);
  % 'plot' excerpt (instead of plotting excerpt data in real experimental
  % time plot vs index and change axis labels to experimental time in s
  % further down)
  for g=1:min(2,size(d,2))
    set(sp.rawExc.excH(g),'xdata',1:diff(ix)+1,'ydata',d(ix(1):ix(2),g));
  end
  mima=[min(min(d(ix(1):ix(2),:))) max(max(d(ix(1):ix(2),:)))];
  tmpxl=[1 diff(ix)+1];
  if length(find(isfinite(wp.excYLim)))==2
    tmpyl=wp.excYLim;
  else
    tmpyl=[mima(1)-(mima(2)-mima(1))/20  mima(2)+(mima(2)-mima(1))/20];
  end
  set(sp.rawExc.axH,'ylim',tmpyl,'xlim',tmpxl);
  % adjust patch so that it fills out entire area
  set(sp.rawExc.patchH,'xdata',tmpxl([1 1 2 2])','ydata',tmpyl([1 2 2 1])');
  % adjust threshold line
  set(sp.rawExc.threshH,'ydata',ap.thresh*[1 1],'xdata',tmpxl);

  % delete previous markers for detected events
  if any(ishandle(sp.rawExc.evMarkH))
    delete(sp.rawExc.evMarkH);
    sp.rawExc.evMarkH=nan;
  end
  % plot new markers
  if ~isempty(wp.evtTsl)
    for ui=1:length(wp.evtTsl)
      % all in the present excerpt, offset corrected for local axis
      tsl=wp.evtTsl{ui}(wp.evtTsl{ui}>=ix(1) & wp.evtTsl{ui}<=ix(2))-ix(1)+1;
      % in case no ts is in current view make up a fake tsl for the plot
      % function so we will in any case get back a handle
      if isempty(tsl)
        tsl=nan;
      end
      tmpy=repmat(tmpyl(1)+diff(tmpyl)*sp.rawExc.EvMarkRelHeight, [length(tsl) 1]);
      pCol=wp.tsMarkerCol{mod(ui-1,length(wp.tsMarkerCol))+1};
      sp.rawExc.evMarkH(ui)=plot(tsl,tmpy,[pCol '.']);
    end
  end

  % delete previous lines for bursts & silent periods
  if any(ishandle(sp.rawExc.burstLh))
    delete(sp.rawExc.burstLh);
    sp.rawExc.burstLh=nan;
  end
  if any(ishandle(sp.rawExc.silentPerLh))
    delete(sp.rawExc.silentPerLh);
    sp.rawExc.silentPerLh=nan;
  end
  % draw lines marking detected bursts
  if ~isempty(wp.buEtsl)
    % all burst with a foot in the present excerpt, offset-corrected for local axis
    etsl=wp.buEtsl;
    stopTsl=sum(etsl,2);
    bIx=~(etsl(:,etslc.tsCol)>ix(2) | stopTsl<ix(1));
    etsl=etsl(bIx,:);
    etsl(:,etslc.tsCol)=etsl(:,etslc.tsCol)-ix(1)+1;
    tmpy=repmat(tmpyl(1)+diff(tmpyl)*sp.rawExc.BurstLineRelHeight, [2 size(etsl,1)]);
    sp.rawExc.burstLh=line(cumsum((etsl(:,[etslc.tsCol etslc.durCol]))',1),tmpy,'color','r','linewidth',2);
  end
  % draw lines marking silent periods
  if ~isempty(wp.buSilentEtsl)
    % all silent periods with a foot in the present excerpt, offset-corrected for local axis
    etsl=wp.buSilentEtsl;
    stopTsl=sum(etsl,2);
    bIx=~(etsl(:,etslc.tsCol)>ix(2) | stopTsl<ix(1));
    etsl=etsl(bIx,:);
    etsl(:,etslc.tsCol)=etsl(:,etslc.tsCol)-ix(1)+1;
    tmpy=repmat(tmpyl(1)+diff(tmpyl)*sp.rawExc.BurstLineRelHeight, [2 size(etsl,1)]);
    sp.rawExc.silentPerLh=line(cumsum((etsl(:,[etslc.tsCol etslc.durCol]))',1),tmpy,'color','b','linewidth',2);
  end

  % if a burst start had been marked manually (wp.buMaStatus==2) delete it,
  % then replot it should it be in current excerpt
  if wp.buMaStatus==2 
    if any(ishandle(sp.rawExc.burstStartMarkH)) 
      delete(sp.rawExc.burstStartMarkH);
      sp.rawExc.burstStartMarkH=nan;
    end
    % x coordinate of marker in pts
    buMaX=cont2discrete(wp.buMaX(1),wp.si/1e3,'intv',0);
    if buMaX>=ix(1) && buMaX<=ix(2)
      sp.rawExc.burstStartMarkH=plot(buMaX-ix(1)+1,sp.rawExc.BurstLineRelHeight,'r^');
    end
  end
  
  % put out information on threshold position expressed in terms of
  % base line noise
  sp.rawExc.textH=setThreshInfoText(wp.baselinePar,ap.thresh,sp.rawExc.textH);
  
  % set order of things: make sure 
  % - patch is in background 
  % - raw data come next (conditioned data overplotting unconditioned)
  % - markers in foreground (=upper positions as children)
  set(sp.rawExc.axH,'children',...
    [setdiff(get(sp.rawExc.axH,'children'),[sp.rawExc.excH; sp.rawExc.patchH]);...
    flipud(sp.rawExc.excH); sp.rawExc.patchH]);
  % x axis labels
  t=discrete2cont(get(gca,'xtick')+ix(1),wp.si/1e6,'intv',0);
  set(gca,'xticklabel',num2str(t','%5.3f'));
  xlabel('time (s)');
end