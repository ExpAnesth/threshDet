function callbackfn_rawexcplot(src,edat)
% ** function callbackfn_rawexcplot(src,edat)
% Callback of ButtonDownFcn of dummy patch (rectangle) in subplot
% sp.rawExc.axH:
% - regular click (SelectionType==´'normal') sets threshold and plots it as
%   a line (or cancels manual burst marking)
% - shift-click (SelectionType=='extend') deletes silent period clicked on 
%   (or cancels manual burst marking)
% - ctrl-click (SelectionType=='alt') sets burst begin or end in manual 
%   burst marking mode 

% historical note: the code relying on lists of transitions (wp.transIA
% and .transitAI) had been developed after the code working with bu and sp
% etsls; therefore, in places, we have a mix of both kinds of code


global wp ap evt bu sp
etslconst

% catch inadvertent clicks on empty subplot
if ~isempty(sp.rawExc.excH)
  subplot(sp.rawExc.axH);
  yl=get(sp.rawExc.axH,'ylim');
  xl=get(sp.rawExc.axH,'xlim');
  % coordinate of mouse pointer
  pointerCo=get(sp.rawExc.axH,'currentpoint');
  % determine kinda click 
  sType=get(gcf,'SelectionType');
  % if it is any click following a ctrl-click remove marker
  if wp.buMaStatus==2 && any(ishandle(sp.rawExc.burstStartMarkH))
    delete(sp.rawExc.burstStartMarkH);
    sp.rawExc.burstStartMarkH=nan;
  end
  % ** if it is a normal click or a shift-click following a ctrl-click
  % cancel manual burst detection (reset wp.buMaStatus to 1) and do NOT
  % execute the normal action associated with the clicks
  if wp.buMaStatus==2 && ismember(sType,{'normal','extend'})
    wp.buMaStatus=1;
    sType='none';
  end
  switch sType
    case 'normal'
      % -------------------------------------------------------------------
      %                      set threshold 
      % -------------------------------------------------------------------
      % get y coordinate of mouse click
      y=pointerCo(1,2);
      % accept only if mouse click inside current y limits
      if y>=yl(1) || y<=yl(2)
        % any click on this window is interpreted as a desire to get rid
        % of the previously determined time stamps, etc., so kill them &
        % associated plots here ***
        cutouts=[];
        evt.tsl=[];
        bu.etsl=zeros(0,etslc.nCol);
        bu.silentEtsl=zeros(0,etslc.nCol);
        wp.evtTsl=[];
        wp.buEtsl=zeros(0,etslc.nCol);
        wp.buSilentEtsl=zeros(0,etslc.nCol);
        wp.manualDeletedSilentPerIx=[];
        wp.transIA=zeros(0,1);
        wp.transAI=zeros(0,1);
        if any(ishandle(sp.rawOv.burstLh))
          delete(sp.rawOv.burstLh);
          sp.rawOv.burstLh=nan;
        end
        if any(ishandle(sp.rawOv.silentPerLh))
          delete(sp.rawOv.silentPerLh);
          sp.rawOv.silentPerLh=nan;
        end
        if any(ishandle(sp.rawOv.evMarkH))
          delete(sp.rawOv.evMarkH);
          sp.rawOv.evMarkH=nan;
        end
        if any(ishandle(sp.rawExc.burstLh))
          delete(sp.rawExc.burstLh);
          sp.rawExc.burstLh=nan;
        end
        if any(ishandle(sp.rawExc.silentPerLh))
          delete(sp.rawExc.silentPerLh);
          sp.rawExc.silentPerLh=nan;
        end
        if any(ishandle(sp.rawExc.evMarkH))
          delete(sp.rawExc.evMarkH);
          sp.rawExc.evMarkH=nan;
        end
        if any(ishandle(wp.legH))
          delete(wp.legH);
          wp.legH=nan;
        end
        cla(sp.iei.axH);
        cla(sp.cutout.axH);

        % now set thresh and accomplish according graphics jobs
        ap.thresh=y;
        set(sp.rawExc.threshH,'ydata',y*[1 1]);

        % put out information on threshold position expressed in terms of
        % base line noise
        sp.rawExc.textH=setThreshInfoText(wp.baselinePar,y,sp.rawExc.textH);
        
        % make sure line is always on top
        tmpChi=get(sp.rawExc.axH,'children');
        tmpChiIx=tmpChi==sp.rawExc.threshH;
        set(sp.rawExc.axH,'children',[tmpChi(tmpChiIx); tmpChi(~tmpChiIx)]);
      end
      
    case 'extend'
      % -------------------------------------------------------------------
      %                delete selected silent period 
      % -------------------------------------------------------------------      
      if isempty(bu.etsl) || isempty(bu.silentEtsl)
        warndlg('shift-click in the excerpt window is an attempt to delete a silent period but time stamp lists are currently empty');
      else
        % get x coordinate of mouse click and convert to scale in etsl (ms)
        x=(pointerCo(1,1)+wp.excXLim_pts(1))/1000*wp.si;
        % find silent period 
        spIx=find(bu.silentEtsl(:,etslc.tsCol)<=x & ...
          sum(bu.silentEtsl(:,[etslc.tsCol etslc.durCol]),2)>x);
        % do nothing if silent period is missed
        if ~isempty(spIx)
          if spIx==1 && bu.silentEtsl(1,etslc.tsCol)<bu.etsl(1,etslc.tsCol)
            warndlg('cannot delete silent period at beginning of recording');
          elseif spIx==size(bu.silentEtsl,1) && bu.silentEtsl(end,etslc.tsCol)>bu.etsl(end,etslc.tsCol)
            warndlg('cannot delete silent period at end of recording');
          else
            % § storing manually deleted sp does not make sense if any
            % burst had been marked manually
            disp(['manually deleting silent period # ' int2str(spIx)]);
            % save this information in ap - it may be useful 
            wp.manualDeletedSilentPerIx=cat(2,wp.manualDeletedSilentPerIx,spIx);
            % find bursts flanking sp in question
            if bu.silentEtsl(1,etslc.tsCol)<bu.etsl(1,etslc.tsCol)
              bIx=[spIx-1 spIx];
            else
              bIx=[spIx spIx+1];
            end
            % ditto for transitions, but remove them in the same step
            wp.transIA(find(wp.transIA>x,1))=[];
            wp.transAI(find(wp.transAI<x,1,'last'))=[];
            % 'glue'
            bu.etsl(bIx(1),etslc.durCol)=sum(bu.etsl(bIx(2),[etslc.tsCol etslc.durCol]),2)-bu.etsl(bIx(1),etslc.tsCol);
            bu.etsl(bIx(2),:)=[];
            bu.silentEtsl(spIx,:)=[];
            % ** don't forget to delete the wp. versions here (will be
            % recreated instantly in job afterBurst)
            wp.buEtsl=zeros(0,etslc.nCol);
            wp.buSilentEtsl=zeros(0,etslc.nCol);
            % also, we must delete current lines (no need to check for
            % their existence; has been done above)
            delete(sp.rawOv.burstLh);
            sp.rawOv.burstLh=nan;
            delete(sp.rawOv.silentPerLh);
            sp.rawOv.silentPerLh=nan;
            delete(sp.rawExc.burstLh);
            sp.rawExc.burstLh=nan;
            delete(sp.rawExc.silentPerLh);
            sp.rawExc.silentPerLh=nan;
            % *** (recursive, in a sense) call to threshdetguifunc
            threshdetguifunc(sp.rawExc.axH,[],{'afterBurst'});
          end
        end
      end

    case 'alt'
      % -------------------------------------------------------------------
      %                     add burst manually 
      % -------------------------------------------------------------------
      % get x coordinate of mouse click and convert to scale in etsl (ms)
      x=(pointerCo(1,1)+wp.excXLim_pts(1))/1000*wp.si;
      % clicks representing burst beginning and end: 
      % - if wp.transIA and wp.transAI have no entries, they may be
      %   anywhere
      % - if bu.etsl and bu.silentEtsl have no entries, they may be
      %   anywhere
      % - if bu.etsl has entries, must both be within a single silent
      %   period OR before first burst OR after last burst
      % - if bu.etsl has no entries but bu.silentEtsl does, they must be
      %   within that silent period
      % find silent period within which click ocurred (if any)
      spIx=find(bu.silentEtsl(:,etslc.tsCol)<=x & ...
        sum(bu.silentEtsl(:,[etslc.tsCol etslc.durCol]),2)>x);
      % find burst within which click ocurred (if any)
      buIx=find(bu.etsl(:,etslc.tsCol)<=x & ...
        sum(bu.etsl(:,[etslc.tsCol etslc.durCol]),2)>x);
      if ~isempty(buIx)
        % if click ocurred within burst, ignore it, and set status to 1
        % (=click will have no consequences, no matter whether it is the
        % first or second one)
        warning('manually marking bursts is only allowed in silent periods'); 
        curBuMaStatus=1;
      elseif isempty(spIx) && ~isempty(bu.silentEtsl) && isempty(bu.etsl)
        % this is one special case mentioned above, separately dealt with
        % here for oversight: if only one silent period had been detected
        % (flanked by bursts not recognized as such because they straddle
        % the borders of the recording) any burst to be marked must be in
        % it
        warning('manually marking bursts is only allowed in silent periods'); 
        curBuMaStatus=1;
      elseif isempty(wp.transAI) && ~isempty(wp.transIA) && x>=max([0 wp.transIA(1)])
        %             ----
        %            |
        %        ----
        warning('manually marking bursts is only allowed left of a single i->a transition'); 
        curBuMaStatus=1;
      elseif ~isempty(wp.transAI) && isempty(wp.transIA) && x<=max([0 wp.transAI(1)])
        %        ----
        %            |
        %             ----
        warning('manually marking bursts is only allowed right of a single a->i transition'); 
        curBuMaStatus=1;
      else
        % if this is the first click within a pair (beginning of burst)...
        if wp.buMaStatus==1
          % § inactivate all buttons?
          % transfer current coordinate
          wp.buMaX(1)=x;
          % add marker
          sp.rawExc.burstStartMarkH=plot(pointerCo(1,1),sp.rawExc.BurstLineRelHeight,'r^');
          % identify rightmost allowed x coordinate given present one and
          % place it in second position of wp.buMaX:
          if isempty(spIx)
            % if click was NOT in identifiable silent period, see whether
            % there is a i->a transition (potential or real burst start)
            % right of current x coordinate
            nextBuIx=find(wp.transIA>x,1);
            if isempty(nextBuIx)
              % if none exists, set to inf (end of recording)
              wp.buMaX(2)=inf;
            else
              % set to that transition
              wp.buMaX(2)=wp.transIA(nextBuIx(1));
            end
          else
            % if click was in identifiable silent period, use its end
            wp.buMaX(2)=sum(bu.silentEtsl(spIx,[etslc.tsCol etslc.durCol]),2);
          end
          disp('ctrl-click to mark burst end, or any other click to cancel this burst');          
          % finally, toggle status
          curBuMaStatus=2;
        elseif wp.buMaStatus==2
          % toggle status regardless of whether click was within legal
          % x limits or not
          curBuMaStatus=1;
          % § re-activate all buttons?

          % check whether click was within allowed interval:
          % - left of limit determined above
          % - right of first click (burst begin)
          if x<wp.buMaX(2) && x>wp.buMaX(1)
            wp.buMaX(2)=x;
            % impute manually marked burst first in tsl of transitions...
            wp.transIA=sort([wp.transIA; wp.buMaX(1)]);
            wp.transAI=sort([wp.transAI; wp.buMaX(2)]);
            % ...then in etsl and silentEtsl 
            % ** note that we cannot use ds.fileInfo.recTime to know the
            % eor (end of recording) input arg into transit2etsl for
            % multi-sweep data (see similar comment in job 'plotEvt' in
            % threshdetguifunc.m)
            [bu.etsl,bu.silentEtsl]=transit2etsl(wp.transIA,wp.transAI,...
              wp.rawNRow/1e3*wp.si);
%               'minActive',ap.minEventWidth,'minInactive',ap.maxBurstGapWidth,...
%             'elimOrder','inactive');
            % just to be on the safe side...
            wp.buMaX=[nan nan];
            % ** don't forget to delete the wp. versions of all etsls (will
            % be recreated instantly in job afterBurst)
            wp.buEtsl=zeros(0,etslc.nCol);
            wp.buSilentEtsl=zeros(0,etslc.nCol);
            % also, we must delete current marker lines for burst & silent
            % periods, if any
            if any(ishandle(sp.rawOv.burstLh))
              delete(sp.rawOv.burstLh);
              sp.rawOv.burstLh=nan;
            end
            if any(ishandle(sp.rawOv.silentPerLh))
              delete(sp.rawOv.silentPerLh);
              sp.rawOv.silentPerLh=nan;
            end
            if any(ishandle(sp.rawExc.burstLh))
              delete(sp.rawExc.burstLh);
              sp.rawExc.burstLh=nan;
            end
            if any(ishandle(sp.rawExc.silentPerLh))
              delete(sp.rawExc.silentPerLh);
              sp.rawExc.silentPerLh=nan;
            end
            % *** (recursive, in a sense) call to threshdetguifunc
            threshdetguifunc(sp.rawExc.axH,[],{'afterBurst'});
          else
            warning('illegal end point of burst was marked');
          end
        end
      end
      wp.buMaStatus=curBuMaStatus;
    
    case 'none'
      % 'intentially left blank'
  end
end

