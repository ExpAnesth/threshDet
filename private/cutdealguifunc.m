function cutdealGUIfunc(src,eventdata,job,varargin)
% ** function cutdealGUIfunc(src,eventdata,job,varargin)
% Collection of callback routines for cutdealGUI.m
%                         >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% job              cell array of char    jobs to accomplish
% sph              struct                handles to subplots in main figure
%                                        window of cutdealgui

% global and persistent data are only partly the same as in threshdetguifunc
global evt wp ap evtCutout sph
persistent wop

pvpmod(varargin);

% §§§ to do:

done=0;
while ~done
  partJob=job{1};
  switch partJob
    case 'init'
      % default layout of cutouts in parallel plot
      if ~all(isfield(wop,{'nRow', 'nCol'}))
        wop.nRow=15;
        wop.nCol=8;
      end
      % number of cutouts to be displayed in parallel plot
      wop.nDispCutout=wop.nRow*wop.nCol;
      % index to cutouts (see further below)
      wop.cutoutIx=[];
      % tags for cutouts
      wop.OKTag=1;    % OK
      wop.killTag=2;  % kick
      % array holding these tags
      wop.cutoutTagArr=[];
      % colors for differently tagged cutouts in same order as above so the
      % tag values can be used for color indexing
      wop.tagCol={[0 0 0],[1 .5 .5]};
      % handles to cutouts
      wop.parPh=nan*zeros(wop.nDispCutout,1);
      % length of cutouts in pts and number of cutouts
      [wop.lenCutout,wop.nCutout]=size(evtCutout{1});
      if numel(evt.tsl{1}) ~= wop.nCutout
        error('internal: number of cutouts and time stamps do not match');
      end
      % placeholder for cutouts needed for initializing parallel plot
      wop.mockData=repmat(nan,wop.lenCutout,1);
      % percentiles needed for estimating good y offset
      wop.prc=prctile(evtCutout{1}(:),[5 50 95]);
      % time axis for cutouts in average plot
      wop.avPlotT=discrete2cont(1:wop.lenCutout,wp.si/1000,'intv',0)+ap.winEvtCutout(1);
      % initialize graphics
      subplot(sph.parall.axH), cla
      set(sph.parall.axH,'color',wp.stdAxCol,'box','on','nextplot','add',...
        'xtick',[]);
      subplot(sph.overlay.axH), cla
      set(sph.overlay.axH,'color',wp.stdAxCol);
      sph.parall.lineH=line(repmat(nan,2,wop.nCol),repmat(nan,2,wop.nCol),...
        'linestyle',':','color','w');
      % next jobs: computing offsets, initializing parallel plot & plotting
      % excerpts 
      job(3:end+2)=job;
      job(1:3)={'computeOffset','initParPlot','plotSame'};
      
    case 'setPlotLayout'
      % (re-)define layout of parallel plot
      if any(wop.cutoutTagArr~=wop.OKTag)
        a=questdlg('Resizing the plot will discard all tags - continue?', ...
          'A second, please','Yes', 'No', 'Yes');
      else
        a='Yes';
      end
      if strcmp(a,'Yes')
        h=findobj('tag', 'setPlotLayout', 'type', 'uicontrol');
        [tmpNum,convOK]=str2num(get(h,'string'));
        if convOK && isequal(size(tmpNum),[1 2]) && all(tmpNum>0)
          wop.nRow=tmpNum(1);
          wop.nCol=tmpNum(2);
          % number of cutouts to be displayed in parallel plot changes
          wop.nDispCutout=wop.nRow*wop.nCol;
          job(3:end+2)=job;
          job(1:3)={'computeAxLim','initParPlot','plotSame'};
        else
          warndlg('number of rows and columns of cutout plot not OK')
          job(1)=[];
        end
      else
        job(1)=[];
      end
      
    case 'computeOffset'
      % x offset (gap+trace length) between cutouts (if they are displayed
      % in columns) in pts
      wop.xOffs=round(wop.lenCutout*1.1);
      % vertical offset between cutouts and y limits of plot: negative -
      % plot traces from top to bottom
      wop.yOffs=-2*mean(abs(diff(wop.prc)));
      job(1)={'computeAxLim'};

    case 'computeAxLim'
      % x and y limits of plot
      wop.xLim=[1 wop.nCol*wop.lenCutout+(wop.nCol-1)*(wop.xOffs-wop.lenCutout)];
      wop.yLim=[wop.prc(1)+(wop.nRow-.5)*wop.yOffs wop.prc(3)-.5*wop.yOffs]; 
      set(sph.parall.axH,'xlim',wop.xLim,'ylim',wop.yLim);
      set(sph.parall.lineH,'ydata',wop.yLim);
      job(1)=[];
      
    case 'initParPlot'
      % initialize/recreate parallel plot:
      % (re-) generate array of handles to excerpts by plotting mock data
      % (as many as cutouts displayed)
      tmp=ishandle(wop.parPh);
      if any(tmp)
        delete(wop.parPh(tmp));
      end
      wop.parPh=nan*zeros(wop.nDispCutout,1);
      subplot(sph.parall.axH)
      % broken lines corresponding to t=0 of cutouts
      sph.parall.lineH=line(...
        [1; 1]*((0:wop.nCol-1)*wop.xOffs+cont2discrete(-ap.winEvtCutout(1),wp.si/1000)),...
        repmat(wop.yLim',1,wop.nCol),...
        'linestyle',':','color','w');
      % plot mock cutouts
      for ct=1:wop.nDispCutout
        % column-dependent x offset
        tmp=floor((ct-1)/wop.nRow)*wop.xOffs;
        wop.parPh(ct)=plot(tmp+1:tmp+wop.lenCutout,wop.mockData,'-');
      end
      % ** set callbacks & color
      set(wop.parPh,'buttondownfcn',{@toggleCutoutStatus},'color',wop.tagCol{wop.OKTag});
      wop.th=title(' ','color','b');
      % *****************************************************************
      % wop.cutoutIx is the index to all cutouts, initialized as [1 2 3..
      % wop.nCutout]. The first wop.nDispCutout entries will be plotted.
      % When the user scrolls through the cutouts wop.cutoutIx will be
      % circularly shifted (in jobs 'plotPrev' and 'plotNext') . If the
      % total number of cutouts is less than the number that can be
      % displayed according to the plot layout wop.cutoutIx will be
      % extended to a length matching that of wop.nDispCutout and contain
      % redundant entries, e.g. [1 2 3.. wop.nCutout 1 2 3..]
      % *****************************************************************
      if wop.nDispCutout<=wop.nCutout
        wop.cutoutIx=1:wop.nCutout;
      else
        tmp=repmat(1:wop.nCutout,1,ceil(wop.nDispCutout/wop.nCutout));
        wop.cutoutIx=tmp(1:wop.nDispCutout);
      end
      % *****************************************************************
      % wop.cutoutTagArr contains the tags for the cutouts (for tag codes see
      % job 'init'). In contrast to wop.cutoutIx it will not be circularly
      % shifted and contains the tags of cutouts in their original
      % (chronological) order irrespective of the display. 
      % *****************************************************************
      wop.cutoutTagArr=repmat(wop.OKTag,wop.nCutout,1);
      job(1)=[];
      
    case 'incYOffs'
      wop.yOffs=wop.yOffs*1.5;
      job(2:end+1)=job;
      job(1:2)={'computeAxLim','plotSame'};
    
    case 'decYOffs'
      wop.yOffs=wop.yOffs/1.5;
      job(2:end+1)=job;
      job(1:2)={'computeAxLim','plotSame'};

    case {'plotPrev','plotNext','plotSame'}
      if strcmp(job{1},'plotPrev')
        wop.cutoutIx=circshift(wop.cutoutIx,[0 wop.nDispCutout]);
      elseif strcmp(job{1},'plotSame')
        % do nothing
      elseif strcmp(job{1},'plotNext')
        wop.cutoutIx=circshift(wop.cutoutIx,[0 -wop.nDispCutout]);
      else
        error('internal: undef parallel plot mode');
      end
      subplot(sph.parall.axH)
      for ct=1:wop.nDispCutout
        % set y data, index into evtCutout as userdata and appropriate
        % color
        set(wop.parPh(ct),'ydata',...
          evtCutout{1}(:,wop.cutoutIx(ct))+mod(ct-1,wop.nRow)*wop.yOffs,...
          'userdata',wop.cutoutIx(ct),...
          'color',wop.tagCol{wop.cutoutTagArr(wop.cutoutIx(ct))});
      end
      if any(diff(wop.cutoutIx(1:wop.nDispCutout))<0)
        if wop.nDispCutout>wop.nCutout
          set(wop.th,'string',['all of ' int2str(wop.nCutout) ' cutouts (with redundancies)']);
        else
          set(wop.th,'string',['cutouts # ' int2str(wop.cutoutIx(1)) ' - ' int2str(wop.cutoutIx(ct)) ' of ' int2str(wop.nCutout) ' (circularly shifted)']);
        end
      else
        set(wop.th,'string',['cutouts # ' int2str(wop.cutoutIx(1)) ' - ' int2str(wop.cutoutIx(ct)) ' of ' int2str(wop.nCutout)]);
      end
      job(1)=[];

    case 'rmEvt'
      % do an irreversible deletion:
      tmpIx=wop.cutoutTagArr==wop.killTag;
      evtCutout{1}(:,tmpIx)=[];
      evt.tsl{1}(tmpIx)=[];
      if isfield(evt,'amp') && ~isempty(evt.amp)
        evt.amp{1}(tmpIx)=[];
      end
      wp.evtTsl{1}(tmpIx)=[];
      job(4:end+3)=job;
      job(1:4)={'init','computeOffset','initParPlot','plotSame'};
      
    case 'averageEvt'
      subplot(sph.overlay.axH), cla, 
      th=ultext('un momento, per favore...',0.02,'fontsize',12,'color','r');      
      % mean of all traces...
      m=mean(evtCutout{1},2);
      ylim=[min(m) max(m)];
      ylim=ylim+diff(ylim)*[-.2 .2];
      % ...but plot only limited number of cutouts
      tmpv1=min(wp.maxNPlotCutout,wop.nCutout);
      tmpix=unique(ceil((1:tmpv1)/tmpv1*wop.nCutout));
      delete(th);
      plot(wop.avPlotT,evtCutout{1}(:,tmpix),'k-');
      hold on
      ph=plot(wop.avPlotT,m,'-');
      set(ph,'color',[.2 1 .2],'linewidth',2);
      axis tight
      set(gca,'ylim',ylim);
      xlabel('time (ms)');
      job(1)=[];

    case 'alignEvt'
      % §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
      % this job has been taken 1:1 from threshdetgui and must still be
      % adjusted
      % §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
      % check for contents of tsl and cutouts
      if isempty(evt.tsl) || isempty(evtCutout)
        warndlg('events cannot be aligned because none were detected');
        job(1)=[];
      elseif numel(evt.tsl)>1
        warndlg('aligning events does not make sense because they seem to have been separated already')
        job(1)=[];
      else
        if ~wp.alignFlag
          cla(sp.cutout.axH);
          % interval (in pts) within which to determine steepest slope:
          % t0 +/- half of pretrig interval
          intv=cont2discrete(-[1/2 3/2]*ap.winEvtCutout(1),wp.si/1000,'intv',1);
          intv=[max(1,intv(1)) min(size(evtCutout{1},1),intv(2))];
          if ap.thresh<=0
            alignMeth='negSlope';
          else
            alignMeth='posSlope';
          end
          % interpolation: target freq is 50 kHz
          ipFac=ceil(wp.si/20);
          evtCutout{1}=halign(evtCutout{1},intv,alignMeth,ipFac);
          wp.alignFlag=true;
          nTs=length(evt.tsl{1});
          pCol=wp.tsMarkerCol{1};
          % plot specified number of cutouts
          tmpv1=min(wp.maxNPlotCutout,nTs);
          tmpix=unique(ceil((1:tmpv1)/tmpv1*nTs));
          time=discrete2cont(1:size(evtCutout{1},1),wp.si/1000,'intv',0)+ap.winEvtCutout(1);
          subplot(sp.cutout.axH); hold on
          plot(time,evtCutout{1}(:,tmpix),pCol);
          axis tight
          set(sp.cutout.axH,'color',wp.stdAxCol);
          xlabel('time (ms)');
        else
          warndlg('events have been aligned already - will only align freshly detected events')
        end
        job(1)=[];
      end

    case 'done'
      % cleanup
      % - remove tslTags
      tmph=findobj('name', 'Processing of Cutouts', 'type', 'figure');
      if ~isempty(tmph)
        delete(tmph);
      end
      job(1)=[];
      
    otherwise
      % error(['illegal job: ' partJob]);
      disp(['job ' partJob ' not yet implemented']);
      job(1)=[];
  end

  done=isempty(job);
end

% callback function (buttondownfcn) of individual cutouts
  function toggleCutoutStatus(src,evt)
    % the figure may hold individual traces more than once, so we need to
    % find them all and set the properties of their handles accordingly
    tmpIx=get(src,'userdata');
    tmpIx2=find(wop.cutoutIx(1:wop.nDispCutout)==tmpIx);
    if wop.cutoutTagArr(tmpIx)~=wop.killTag
      wop.cutoutTagArr(tmpIx)=wop.killTag;
      set(wop.parPh(tmpIx2),'color',wop.tagCol{wop.killTag});
    else
      wop.cutoutTagArr(tmpIx)=wop.OKTag;
      set(wop.parPh(tmpIx2),'color',wop.tagCol{wop.OKTag});
    end
  end

end


% figure(findobj('name', 'Processing of Cutouts', 'type', 'figure'));
