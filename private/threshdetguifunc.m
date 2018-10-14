function threshdetguifunc(~,~,job,varargin)
% ** function threshdetguifunc(~,~,job,varargin)
% Collection of callback routines for threshdetGUI.m
%                         >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% job              cell array of char    jobs to accomplish
% sp               struct                handles to subplots of main
%                                         figure

% -------------------------------------------------------------------------
% Version 5.10, October 2018
% (C) Harald Hentschke (University Hospital of Tuebingen)
% -------------------------------------------------------------------------

% to do:
% - at one point, do a serious cleanup: rename variables, delegate the
% various jobs to functions, rethink the concept of using global variables
% - 'events' (read spikes) were originally believed to be sufficiently
% represented by tsl only; however, in the meantime we have sweep indices
% (evt.sweepIx) and also amplitudes (evt.amp) which may characterize them.
% Instead of having these parameters as separate fields of evt it would be
% much cleaner to have an etsl for events, too, and use this to hold all
% relevant information on the events. The downside, however, is that in
% the case of extracellularly recorded spikes, for which threshdetui was
% originally conceived, we'd be wasting memory. Maybe having additional
% fields of evt isn't such a bad idea after all?
% - multiple-threshold detection (+alignment) to capture 'emerging' units
% of different polarity
% - plotBurst (equivalent of plotEvt for bursts)
% - minor things --> see the '§' characters within code for sections to be
% improved
% - systematic 'intermediate cleanup' func to be called when thresh
% changes, new data loaded, etc.
% - option to remove individual bursts or events via mouse click (buttons
% close to subplots which set different callback functions for zoom
% rectangle)

% We need some global data:
%   evt.tsl=cell array of time stamp lists of single 'events' 
%   bu.etsl=extended time stamp list of phases of activity ('bursts')
%   evtCutout, buCutout = cutouts of events/bursts (=fixed/variable length) 
% Strucures:
%   ds='data set' listing properties of current file
%   ap='analysis parameters' describing details of current analyis
%   wp='working parameters' (like colors)
%   sp=subplot handles
global evt bu ds ap wp sp evtCutout buCutout

% We need persistent data:
%   d=raw data
persistent d  

pvpmod(varargin);
etslconst;

jobsToDo=true;
while jobsToDo
  if ~isempty(wp) && wp.isBatchMode
    % check if current batch job is to be cancelled and set job accordingly
    if get(wp.mainGuiHandles.cancelBatchJobBttn,'value')~=0
      disp('cancellation of batch job requested...');
      job={'killBatchJob'};
    end
  end
  
  partJob=job{1};
  switch partJob
    case 'init'
      % *******************************************************************
      % in this ini job, all fields of wp, ap and ds are set up. This is a
      % necessity because the callbacks set here are handles to functions
      % expecting fully defined variables. Besides, it is helpful to have
      % an overview of all fields of above-mentioned key variables. Last,
      % but not least, default values for some parameters are set.
      % *******************************************************************
      % ----- set up major variables for raw data and time stamp list
      d=[];
      evtCutout=[];
      buCutout=[];
      evt.tsl=[];
      evt.amp=[];
      evt.tRise=[];
      bu.etsl=zeros(0,etslc.nCol);
      bu.silentEtsl=zeros(0,etslc.nCol);
      % --------------------------
      % ----- set up ds (data set)
      % --------------------------
      % uigetfile returns 0 as output arguments if the cancel button is
      % hit, so use these as initializing values
      ds.dataFn=0;
      ds.dataPath='e:\_data\otc_ctx\';
      % will contain a number of file-specific parameters after upload
      ds.fileInfo=[];
      
      % -------------------------------------
      % ----- set up ap (analysis parameters)
      % -------------------------------------      
      % name of major results file
      ap.resFn='lastFile_res';
      % name of file holding preconditioned raw data
      ap.rawFn='lastFile_raw';
      % name of file holding evtCutout
      ap.evtCutoutFn='lastFile_evtCutout';
      % name of file holding buCutout
      ap.buCutoutFn='lastFile_buCutout';
      % name of ascii file holding time stamp list of events (in case
      % events will be separated into units a specifier will be appended)
      ap.asciiEvtFn='lastFile_evtTsl';

      % ~~~~~~~ preconditioning section
      % artifact elimination
      ap.afElimination=[];
      % lowpass filter freq
      ap.loCFreq=nan;
      % highpass filter freq
      ap.hiCFreq=nan;
      % notch center freq
      ap.notchCFreq=nan;
      % settings for differentiator filter
      ap.differfi=[];
      % downsampling factor
      ap.sampFac=nan;
      % custom command
      ap.customCommand='why';
      % the question whether cutouts shall be produced from partially
      % conditioned raw data
      ap.doUncondCutout=false;
      
      % ~~~~~~~ event detection section
      % detection threshold (absolute, that is, in magnitude units)
      ap.thresh=nan;
      % detection threshold (relative, that is, relative to noise level)
      ap.relativeThresh=nan;
      % event detection method (see tcd)
      ap.evtDetMode='cross';
      % fraction of threshold; determines burst ends 
      ap.threshFract=1;
      % dead time
      ap.evDeadT=nan;
      % cutout window (ms)
      ap.winEvtCutout=[-1 2];
      % extension of cutout window (ms)
      ap.winBuCutout=[-10 10];

      % ~~~~~~~ burst detection section
      % minimal width of an event to be accepted as part of a burst
      ap.minEventWidth=0;
      % gap length (ms) separating series of events into bursts 
      ap.maxBurstGapWidth=0;
      % input parameter 'threshLag' into tbt (ms)
      ap.threshLag=0;

      % --------------------------------------
      % ----- set up wp ('working' parameters)
      % --------------------------------------
      % ~~~~~~~ display options & matlab version section
      % which version of Matlab?
      wp.mver=ver;
      % note that in standalone deployed code function ver may produce
      % several entries with .Name equal to 'Matlab', so we have to opt for
      % one
      tmpIx=find(strcmpi('matlab',{wp.mver.Name}),1);
      wp.mver=str2double(wp.mver(tmpIx).Version);
      % which version of threshdetgui?
      wp.ver=5.10;
      % length (ms) of raw data excerpt displayed in second subplot from top
      wp.excLen=1000;
      % same in pts
      wp.excLen_pts=nan;
      % x limits of excerpt in points
      wp.excXLim_pts=[nan nan];
      % y limits of excerpt (mV) - if nan limits will enclose all data
      wp.excYLim=[nan nan]; 
      % max number of cutouts to plot
      wp.maxNPlotCutout=100;
      % standard background color of subplots in main figure window
      wp.stdAxCol=[.7 .7 .6];
      % color of line representing threshold
      wp.threshcol='y';
      % color of unconditioned data trace in overview and excerpt plots (in
      % case it is displayed along conditioned data)
      wp.uncondCol=[.8 .8 .8];
      % colors of markers representing spx of different units (will be
      % cycled through in case of more units than colors
      wp.tsMarkerCol={'b','r','g','c','m','y','k'};
      % handle to legend for display of firing rates
      wp.legH=nan;
      % flag for preconditioning of raw data
      wp.precondFlag=false;
      % inter-event-interval histogram (ieih): bin width
      wp.ieiBinWid=1;
      % ieih: upper bin limit (bin centers for inter-spike-interval
      % histogram will be computed later)
      wp.ieiBinULim=100;
      % stop band attenuation of highpass filter (will change depending on
      % setting of corner frequency by user)
      wp.hiRs=5;
      % ~~~~~~~ saving options section
      wp.saveEvt=0;
      wp.saveBu=0;
      wp.saveEvtCutout=0;
      wp.saveBuCutout=0;
      wp.savePrecondData=0;
      wp.resFnString='';
      wp.notesString='';
      % ----- the following working pars are not accessible in the options
      % dialog:
      % handles to main window and its children
      wp.mainFigureHandle=findobj('Tag','threshdetgui','type','figure');
      wp.mainGuiHandles=guihandles(wp.mainFigureHandle);
      % directory for parameter files
      if isdeployed
        % in deployed mode, set to a fail-safe location
        wp.pfDir='c:';
      else
        wp.pfDir=strrep(mfilename('fullpath'),['private' filesep mfilename],'parameterFiles');
      end
      % flag determining whether batch mode is on
      wp.isBatchMode=false;
      % list of files to be processed in batch mode
      wp.batchFileList={};
      % field names of handles to objects which are supposed to be hidden
      % while batch mode runs
      wp.batchHideOHFNames=fieldnames(wp.mainGuiHandles);
      tmpIx=false(size(wp.batchHideOHFNames));
      for g=1:numel(wp.batchHideOHFNames)
        tmpIx(g)=~isempty(strfind(wp.batchHideOHFNames{g},'Bttn')) &&...
          ~strcmp(wp.batchHideOHFNames{g},'cancelBatchJobBttn');
      end
      wp.batchHideOHFNames=wp.batchHideOHFNames(tmpIx);
      % a containers for event numbers and rates of batch processed files
      % for output on command line
      wp.batchEvNumber=[];
      wp.batchEvRate=[];
      % a list of physical units and the factor by which the data will be
      % multiplied 
      wp.unitList={...
        'µV', .001;
        'mV', 1;
        'V', 1000;
        'pA', 1;
        'nA',1000};
      wp.scaleFac=1;
      % number of rows, columns and slices and name of current channel
      wp.rawNRow=0;
      wp.rawNCol=0;
      wp.rawNSli=0;
      % for gap-free (=continuously recorded) data the following two
      % parameters are nan; for multi-sweep (episodic) data they will be
      % set to the number of points per episode and number of episodes,
      % respectively
      wp.rawNEpisodePt=nan;
      wp.rawNEpisode=nan;      
      wp.dataChanName={''};
      wp.si=nan;
      % base line parameters
      wp.baselinePar=[nan nan];
      % x coordinate of last mouse click in excerpt window
      wp.curPtX=0;
      % current status of manual burst markup
      wp.buMaStatus=1;
      % burst markup x coordinates
      wp.buMaX=[nan nan];
      % time stamp list of events in ticks (for excerpt plot)
      wp.evtTsl=[];
      % time stamp list of bursts in ticks (for excerpt plot)
      wp.buEtsl=zeros(0,etslc.nCol);
      % time stamp list of silent periods in ticks (for excerpt plot)
      wp.buSilentEtsl=zeros(0,etslc.nCol);
      % a list of manually, sequentially eliminated silent periods
      wp.manualDeletedSilentPerIx=[];
      % list of 'raw' transitions
      wp.transIA=zeros(0,1);
      wp.transAI=zeros(0,1);

      % ----- initialize subplots
      % -- raw, overview 
      % plot is only partially initialized here because 
      % function pllplot does funny things to the axis and its children
      subplot(sp.rawOv.axH), cla, hold on
      set(sp.rawOv.axH,'color',wp.stdAxCol);
      tmpyl=get(gca,'ylim');
      tmpxl=get(gca,'xlim');
      % handle to patch object whose ButtonDownFcn callback generates an
      % excerpt plot (callbackfn_rawovplot)
      sp.rawOv.patchH=patch(tmpxl([1 1 2 2])',tmpyl([1 2 2 1])',wp.stdAxCol);
      % setting the callback here does not make sense because variables sp
      % and d are not yet set up (fully or at all)
      % set(sp.rawOv.patchH,'ButtonDownFcn',{@callbackfn_rawovplot,d});
      % lines indicating bursts...
      sp.rawOv.burstLh=nan;
      % ...and silent periods
      sp.rawOv.silentPerLh=nan;
      % relative height of these lines in the plot
      sp.rawOv.BurstLineRelHeight=.02;
      % markers indicating events
      sp.rawOv.evMarkH=nan;
      % relative height of these markers in the plot 
      sp.rawOv.EvMarkRelHeight=.98;

      % -- raw, excerpt 
      subplot(sp.rawExc.axH), hold on,
      set(sp.rawExc.axH,'color',wp.stdAxCol);
      % excerpt plots
      sp.rawExc.excH=plot(ones(2)*nan,'k');
      % plot the dummy patch
      tmpyl=get(gca,'ylim');
      tmpxl=get(gca,'xlim');
      sp.rawExc.patchH=patch(tmpxl([1 1 2 2])',tmpyl([1 2 2 1])',wp.stdAxCol);
      % plot invisible threshold line 
      sp.rawExc.threshH=line(tmpxl,ap.thresh*[1 1],'color',wp.threshcol,...
        'linestyle','--','linewidth',1.5);
      % lines indicating bursts
      sp.rawExc.burstLh=nan;
      % same for silent periods
      sp.rawExc.silentPerLh=nan;
      % relative height of these lines in the plot
      sp.rawExc.BurstLineRelHeight=.02;
      % markers indicating events
      sp.rawExc.evMarkH=nan;
      % relative height of these markers in the plot 
      sp.rawExc.EvMarkRelHeight=.98;
      % marker indicating manually marked burst begin
      sp.rawExc.burstStartMarkH=nan;
      % info text (currently: threshold expressed as amplitude of base line
      % noise
      sp.rawExc.textH='';
      % set the callback (*** important: all variables/fields to be accessed in
      % callback function MUST EXIST AT THIS POINT ***)
      set(sp.rawExc.patchH,'ButtonDownFcn',{@callbackfn_rawexcplot});
      % -- cutouts
      set(sp.cutout.axH,'color',wp.stdAxCol);
      % -- iei
      set(sp.iei.axH,'color',wp.stdAxCol);
      job(1)=[];
      
    case 'openOptionsDialog'
      threshdet_optgui;
      job(1)={'writeOptions2Gui'};
      
    case 'writeOptions2Gui'
      % set the various uicontrols' strings and values to those of
      % corresponding fields of ap and wp, all along checking for errors
      % get handles to all uicontrols via uihandles.m ...
      handles=guihandles(findobj('tag','options'));
      uicFn=fieldnames(handles);
      apFn=fieldnames(ap);
      wpFn=fieldnames(wp);
      structName={'ap','wp'};
      % ..and set their 'string' properties to the values of the
      % matching fields of ap or wp
      for g=1:length(uicFn)
        structIx=[~isempty(strmatch(uicFn{g},apFn,'exact')),...
          ~isempty(strmatch(uicFn{g},wpFn,'exact'))];
        if length(find(structIx))==1
          eval(['cType=get(handles.' uicFn{g} ',''style'');']);
          switch cType
            case 'edit'
              % depending on the type of the field...
              switch uicFn{g}
                case {'customCommand','resFnString','notesString','afElimination','differfi','evtDetMode'}
                  eval(['set(handles.' uicFn{g} ',''string'',' structName{structIx} '.' uicFn{g} ');']);
                case {'thresh'}
                  % threshold must be written back with comparatively high
                  % precision
                  eval(['set(handles.' uicFn{g} ',''string'',num2str(' structName{structIx}  '.' uicFn{g} ',''% 2.5f''));']);
                case {'excYLim'}
                  % same for y limits of excerpt plot
                  eval(['set(handles.' uicFn{g} ',''string'',num2str(' structName{structIx}  '.' uicFn{g} ',''% 4.3f''));']);
                otherwise
                  eval(['set(handles.' uicFn{g} ',''string'',num2str(' structName{structIx}  '.' uicFn{g} ',''% 6.2f''));']);
              end
            case 'checkbox'
              eval(['set(handles.' uicFn{g} ',''value'',' structName{structIx}  '.' uicFn{g} ');']);
            case 'text'
              % do nothing because text uicontrols do not hold any
              % information
            otherwise
              errordlg('internal: encountered uicontrol other than edit, checkbox or text');
          end
        elseif length(find(structIx))>1
          errordlg(['internal: uicontrol tagged ''' uicFn{g} ''' has more than one matching fields in ap and wp']);
        end
      end
      job(1)=[];

    case 'readOptionsFromGui'
      % the inverse of job 'writeOptions2Gui': retrieve the various 
      % uicontrols' strings and values and set corresponding fields of ap 
      % and wp. All checks for major pitfalls are done in 'writeOptions2Gui' 
      % so they are omitted here.
      handles=guihandles(findobj('tag','options'));
      uicFn=fieldnames(handles);
      apFn=fieldnames(ap);
      wpFn=fieldnames(wp);
      structName={'ap','wp'};
      for g=1:length(uicFn)
        structIx=[~isempty(strmatch(uicFn{g},apFn,'exact')),...
          ~isempty(strmatch(uicFn{g},wpFn,'exact'))];
        if length(find(structIx))==1
          eval(['cType=get(handles.' uicFn{g} ',''style'');']);
          switch cType
            case 'edit'
              % depending on the type of the field...
              switch uicFn{g}
                case {'customCommand','resFnString','notesString','afElimination','differfi','evtDetMode'}
                  % §§for unknown reasons differfi returns a cell, so catch
                  % that
                  tmpString=get(handles.(uicFn{g}),'string');
                  if iscell(tmpString)
                    tmpString=tmpString{1};
                  end
                  eval([structName{structIx} '.' uicFn{g} '=tmpString;']);
                  % eval([structName{structIx} '.' uicFn{g} '=get(handles.' uicFn{g} ',''string'');']);
                otherwise
                  eval(['[tmpNum,conversionOK]=str2num(get(handles.' uicFn{g} ',''string''));']);
                  if conversionOK
                    eval([structName{structIx} '.' uicFn{g} '=tmpNum;']);
                  else
                    warndlg(['could not read value for ' structName{structIx} '.' uicFn{g} ' - only numeric values are allowed (typographic error?)'])
                  end
                  % § to do: if threshold or dead time were changed delete
                  % the results of the previous event detection (if it
                  % exists) and plot the threshold line accordingly. Do the
                  % same for burst detection. Probably it would be best to 
                  % dedicate a nested function to this
              end
            case 'checkbox'
              eval([structName{structIx}  '.' uicFn{g} '=get(handles.' uicFn{g} ',''value'');']);
            otherwise
          end
        end
      end
      job(1)={'digestOptions'};

    case 'writeOptions2File'
      % *** it is extremely important not to dump wp/ap to file because
      % each of these has fields which cannot be set in the options dialog
      % but are determined e.g. after upload of fresh raw data. Instead,
      % the whole figure including the uicontrols will be saved (and loaded
      % by 'readOptionsFromFile'). This is certainly less than elegant, but 
      % it is relatively fail-safe
      [tmpDataFn,tmpDataPath] = uiputfile([wp.pfDir filesep '*.fig']);
      if ischar(tmpDataFn) && ischar(tmpDataPath) 
        saveas(findobj('tag','options'),[tmpDataPath tmpDataFn ],'fig');
      end
      job(1)=[];
      
    case 'readOptionsFromFile'
      % look for files in default dir
      [tmpOptFn,tmpOptPath] = uigetfile([wp.pfDir filesep '*.fig'],'pick options file');
      if ischar(tmpOptFn) && ischar(tmpOptPath)
        close(findobj('tag','options'));
        openfig([tmpOptPath tmpOptFn]);
        % set wp.pfDir
        wp.pfDir=tmpOptPath;
      end
      job(1)=[];

    case 'digestOptions'
      % this part job must be run whenever i) options were read from gui
      % and ii) data were loaded
      disp('** processing & checking options..');
      % ----- checks of parameters:
      % --- artifact elimination: check for syntax, all other checks will
      % be performed in function elim_artifact
      if ~isempty(ap.afElimination)
        try
          tmp=eval(ap.afElimination);
          if ~iscell(tmp)
            warndlg('input in field ''artifact elimination'' must evaluate to a cell array - not eliminating artifacts');
            ap.afElimination={};
          end
        catch
          warndlg('faulty syntax in field ''artifact eliminiation'': must evaluate to a cell array - not eliminating artifacts');
          ap.afElimination={};
        end
      end
      % --- differentiator filter: check for syntax, all other checks will
      % be performed in function differfi
      if ~isempty(ap.differfi)
        try
          tmp=eval(ap.differfi);
          if ~iscell(tmp)
            warndlg('input in field ''differentiator filter'' must evaluate to a cell array - not filtering');
            ap.differfi={};
          end
        catch
          warndlg('faulty syntax in field ''differentiator filter'': must evaluate to a cell array - not filtering');
          ap.differfi={};
        end
      end
      % --- downsampling factor
      if ~isfinite(ap.sampFac)
        ap.sampFac=1;
      else
        if ap.sampFac<=0
          warndlg('downsampling factor is zero or negative - data will not be downsampled');
          ap.sampFac=1;
        end
      end
      % --- filter parameters: can only be checked if sampling interval is
      % known (=data loaded)
      if ~isempty(d)
        % Highpass filtering at low cutoff frequencies and with the default
        % value of Rs (stopband attenuation, see buttord) in hifi (30
        % dB/octave) produces ugly wobbles or even fails completely.
        % Therefore, implement a filter with Rs depending on the cutoff
        % frequency in a linear way: start with a minimal Rs of 5 dB/octave
        % at 0 Hz and increase linearly to 40 dB/octave at 30 Hz:
        wp.hiRs=5+min(40,5+(40-5)/30);
      end
      % --- threshold section      
      if ap.threshFract<=0 || ap.threshFract>1
        warndlg('threshold fraction must be between 0 and 1');
      end
      if ap.relativeThresh<2 && ap.relativeThresh>-2
        warndlg('a relative threshold value in [-2 2] usually does not make sense');
      end
      % --- event detection
      if length(ap.winEvtCutout)~=2
        warndlg('event cutout window is an interval, so it must contain two values')
      elseif diff(ap.winEvtCutout)<=0
        warndlg('event cutout window must contain two values, the left one being lower than the right one');
      end
      % --- burst detection
      if length(ap.winBuCutout)~=2
        warndlg('burst ''cutout win extension'' must contain two values')
      elseif diff(ap.winBuCutout)<=0
        warndlg('burst ''cutout win extension'' must contain two values, the left one being lower than the right one');
      end
      if ap.threshLag<0 || ~isfinite(ap.threshLag)
        warndlg('threshLag must be zero or positive finite (set to zero to ignore this criterion)');
      end
      % --- display options:
      % excerpt length in pts may need to be recomputed
      wp.excLen_pts=cont2discrete(wp.excLen,wp.si/1000,'intv',1);
      if any(isfinite(wp.excYLim))
        if length(wp.excYLim)~=2
          warndlg('''y limits excerpt'' must contain two values')
        elseif diff(wp.excYLim)<=0
          warndlg('''y limits excerpt'' must contain two values in ascending order')
        end
      end
      % bin centers for inter-spike-interval histogram
      wp.ieiBin=wp.ieiBinWid/2:wp.ieiBinWid:wp.ieiBinULim;
      if isempty(wp.ieiBin)
        warndlg('check settings for Inter-Event Interval bins (IEI)');
      end
      if ~isempty(wp.resFnString)
        wp.resFnString=deblank(wp.resFnString);
      end
      job(1)=[];      

    case 'closeOptionsDialog'
      close(findobj('tag','options'));
      job(1)=[];      
      
    case 'readData'
      % ******************************************************************
      % delete/reset results/plots/variables from last file/channel
      d=[];
      evtCutout=[];
      buCutout=[];
      wp.si=nan;
      wp.excXLim_pts=[nan nan];
      % empty tsl
      evt.tsl=[];
      evt.amp=[];
      evt.tRise=[];
      bu.etsl=zeros(0,etslc.nCol);
      bu.silentEtsl=zeros(0,etslc.nCol);
      wp.evtTsl=[];
      wp.buEtsl=zeros(0,etslc.nCol);
      wp.buSilentEtsl=zeros(0,etslc.nCol);
      wp.manualDeletedSilentPerIx=[];
      wp.transIA=zeros(0,1);
      wp.transAI=zeros(0,1);
      wp.baselinePar=[nan nan];
      % delete all lines & markers & set handles to nan
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
      if any(ishandle(sp.rawExc.burstStartMarkH))
        delete(sp.rawExc.burstStartMarkH);
        sp.rawExc.burstStartMarkH=nan;
      end
      if any(ishandle(sp.rawExc.textH))
        delete(sp.rawExc.textH);
        sp.rawExc.textH=nan;
      end
      if any(ishandle(wp.legH))
        delete(wp.legH);
        wp.legH=nan;
      end
      % wipe plots
      cla(sp.rawOv.axH)
      % for some mysterious reason clearing the overview plot makes it
      % invisible, so revert that
      set(sp.rawOv.axH,'visible','on');
      set(sp.rawExc.excH,'xdata',nan*[1 1]');
      set(sp.rawExc.excH,'ydata',nan*[1 1]');
      cla(sp.cutout.axH);
      cla(sp.iei.axH);
      % (re-)set preconditioning flag
      wp.precondFlag=false;
      % reset color of unconditioned data trace to black
      set(sp.rawExc.excH(1),'color','k');
      % drawnow;
      % ******************************************************************
      % thresh is deliberately left at old value, as is the channel name
      % ******************************************************************
      
      % if we are in batch mode pick file from list
      if wp.isBatchMode
        tmpDataPath=ds.dataPath;
        % extract & kick current file from list
        tmpDataFn=wp.batchFileList{1};
        wp.batchFileList(1)=[];
      else
        % regular mode: ui interface for file selection
        [tmpDataFn,tmpDataPath] = uigetfile({'*.abf';'*.h5';'*.mat'},'pick data file',ds.dataPath,'MultiSelect','on');
        if iscell(tmpDataFn)
          % multiple files were selected for batch detection: set flag 
          wp.isBatchMode=true;
          % file list (sorted, because some OS don't do that)
          wp.batchFileList=sort(tmpDataFn);
          % prevent further execution of current partJob by setting tmpDataFn to []
          tmpDataFnA=[];
          % next two tasks: a placeholder (will be deleted further below in present partJob)
          % and the batch job-creating partJob
          job(2:end+1)=job;
          job(1:2)={'nada','createBatchJob'};
        end
      end
      
      if ischar(tmpDataFn) && ischar(tmpDataPath)
        ds.dataFn=tmpDataFn;
        ds.dataPath=tmpDataPath;
        % retrieve information about file
        [~,~,dataFileExt]=fileparts(ds.dataFn);
        switch dataFileExt
          case '.abf'
            [~,~,ds.fileInfo]=abfload([ds.dataPath ds.dataFn],'info');
          case '.mat'
            load([ds.dataPath ds.dataFn],'fi');
            ds.fileInfo=fi;
          case '.h5'
            [~,~,ds.fileInfo]=read_mcdhdf5([ds.dataPath ds.dataFn],'info');
          otherwise
            % this case should not be reached because of the filter in
            % uigetfile
            errordlg('illegal data file type')
            return
        end
        % keep previously chosen channel (if any) as default selection for
        % subsequent files
        if ~isempty(wp.dataChanName{1})
          chix=strmatch(wp.dataChanName{1},ds.fileInfo.recChNames,'exact');
          % in case no common channels exist default to 1
          if isempty(chix)
            chix=1;
          end
        else
          chix=1;
        end
        % if not in batch mode, open dialog for channel selection
        if ~wp.isBatchMode
          chix=picklistitem(strvcat(ds.fileInfo.recChNames),'defaultVal',chix);
        end
        % *** once complete file/channel info is available set some wp vars
        wp.dataChanName=ds.fileInfo.recChNames(chix);
        % determine whether information about physical units is available
        % and if so, possibly convert data
        if isfield(ds.fileInfo,'recChUnits')
          ix=strmatch(ds.fileInfo.recChUnits{chix},wp.unitList(:,1));
          % if unit not found don't scale data and warn
          if isempty(ix)
            warndlg(['physical units on chosen channel unknown (' ds.fileInfo.recChUnits{chix} ') - data will not be scaled']);
          else
            wp.scaleFac=wp.unitList{ix,2};
          end
        else
          warndlg(['information on physical units on chosen channel unavailable - data will not be scaled']);
        end
        % copy sampling interval to wp because due to downsampling it may
        % change
        wp.si=ds.fileInfo.si;
        % length of raw data excerpt in points
        wp.excLen_pts=cont2discrete(wp.excLen,wp.si/1000,'intv',1);
        % load data
        switch dataFileExt
          case '.abf'
            d=abfload([ds.dataPath ds.dataFn],'channels',wp.dataChanName);
          case '.mat'
            d=matDload([ds.dataPath ds.dataFn],'channels',wp.dataChanName);
          case '.h5'
            d=read_mcdhdf5([ds.dataPath ds.dataFn],'channels',wp.dataChanName);
        end
        % first thing to do: scale data
        if wp.scaleFac~=1
          d=d*wp.scaleFac;
        end
        % at this point d can only have one column (=data from one channel)
        % but more than one 'slices' (=multi-sweep data)
        [wp.rawNRow,wp.rawNCol,wp.rawNSli]=size(d);
        % ** if there is more than one 'slice' (=d is truly 3D) we're
        % dealing with multi-sweep, fixed-length data
        if wp.rawNSli>=2
          warning('data file contains discontinuous, multi-sweep data - the sweeps will be concatenated (association of events with sweep number will be retained)');
          wp.rawNEpisodePt=wp.rawNRow;
          wp.rawNEpisode=wp.rawNSli;
          d=d(:);
          % repeat
          [wp.rawNRow,wp.rawNCol,wp.rawNSli]=size(d);
        else
          wp.rawNEpisodePt=nan;
          wp.rawNEpisode=nan;
        end
        % next two jobs: plotting data & checking options
        job(2:end+1)=job;
        job(1:2)={'plotOv','digestOptions'};
        clear tmp*
      else
        job(1)=[];
      end
      
    case 'createBatchJob'      
      disp('starting batch job...');
      % make batch job control button surface
      set(wp.mainGuiHandles.cancelBatchJobBttn,...
        'visible','on',...
        'enable','on');
      % disable all other buttons
      for bIx=1:numel(wp.batchHideOHFNames)
        set(wp.mainGuiHandles.(wp.batchHideOHFNames{bIx}),'enable','off');
      end
      % change main window's background color
      set(wp.mainFigureHandle,'color',[.6 .6 .6]);
      % temporary containers for event number and rate of batch of files
      % for quick display at end of job
      wp.batchEvRate=nan(numel(wp.batchFileList),1);
      wp.batchEvNumber=wp.batchEvRate;
      % for each file to be processed, this is the sequence of jobs that
      % must be run
      tmpJob={'readData','preconditionRawData','detEvt','saveResults'};
      % replicate this set of jobs an according number of times, in the
      % process overwriting variable job
      job=repmat(tmpJob,1,numel(wp.batchFileList));
      
    case 'killBatchJob'
      % last file has been processed or usere requested cancellation of
      % batch job - set parameters accordingly
      wp.isBatchMode=false;
      wp.batchFileList={};
      % reset main window's background color
      set(wp.mainFigureHandle,'color',[.9 .9 .9]);
      % disable cancel button
      set(wp.mainGuiHandles.cancelBatchJobBttn,...
        'value',0,...
        'visible','off',...
        'enable','off');
      % re-enable all other buttons
      for bIx=1:numel(wp.batchHideOHFNames)
        set(wp.mainGuiHandles.(wp.batchHideOHFNames{bIx}),'enable','on');
      end
      disp('event numbers:');
      % ** note flipud
      disp(flipud(wp.batchEvNumber));
      wp.batchEvNumber=[];
      disp('event rates (Hz):');
      disp(flipud(wp.batchEvRate));
      wp.batchEvRate=[];

      disp('*** batch job done');
      job=cell(1,0);
      
    case 'plotOv'
      figure(findobj('Tag','threshdetgui','type','figure'));
      subplot(sp.rawOv.axH), cla, hold on
      if ~isempty(d)
        % unconditioned data...
        [~,~,~,ph]=pllplot(d(:,1),'si',wp.si,'noscb',1);
        tmpNShift=1;
        % data in first column (unconditioned trace) possibly in other color
        if wp.precondFlag
          set(ph,'color',wp.uncondCol);
          % conditioned data
          pllplot(d(:,2),'si',wp.si,'noscb',1);
          tmpNShift=2;          
        end
        niceyax;
        % dummy patch has to be re-created because pllplot does funny
        % things to the axis
        tmpyl=get(gca,'ylim');
        tmpxl=get(gca,'xlim');
        sp.rawOv.patchH=patch(tmpxl([1 1 2 2])',tmpyl([1 2 2 1])',wp.stdAxCol);
        % set the callback
        set(sp.rawOv.patchH,'ButtonDownFcn',{@callbackfn_rawovplot,d});
        % information on file
        th=ultext([ds.dataFn ', ' wp.dataChanName{1}],0.005,'color','b','fontsize',10,'fontweight','bold','interpreter','none');
        % put patch to background without changing order of other children
        chi=get(sp.rawOv.axH,'children');
        tmpIx=chi==sp.rawOv.patchH;
        set(sp.rawOv.axH,'children',[chi(~tmpIx);sp.rawOv.patchH]);
        % plot excerpt starting at beginning of trace
        callbackfn_rawovplot(sp.rawExc.axH,[],d,'x',0);
      end
      job(1)=[];
      
    case 'preconditionRawData'
      % the logics behind preconditioning:
      % - only if wp.precondFlag is false (data have been untouched) will the
      % data be preconditioned
      % - column 1 will hold raw data
      % - column 2 will hold the fully preconditioned data
      % - column 3 will hold partly preconditioned data (artifact-cleaned,
      % lowpass-filtered and downsampled, but nothing else further down the
      % line; may be needed for generation of cutouts and determination of
      % event amplitude)
      % - if d is originally episodic data it will be temporarily reshaped
      % into its original form, i.e. different sweeps will reside in
      % different 'slices' because filtering, resampling etc. may result in
      % funny artifacts if it is done on concatenated sweeps. After
      % preconditioning d will be reshaped into the default concatenated
      % form.
      if ~isempty(d)
        if ~wp.precondFlag
          disp('** preconditioning data...');
          wp.precondFlag=true;
          % set color of unconditioned data to its background hue
          set(sp.rawExc.excH(1),'color',wp.uncondCol);
          % 0. copy
          tmpD=d;
          % if we're dealing with episodic data these have to be reshaped
          % temporarily for (most of) the preconditioning to work properly
          if ~isnan(wp.rawNEpisode) && wp.rawNEpisode>1
            % reshape into 3D
            d=reshape(d,[wp.rawNEpisodePt 1 wp.rawNEpisode]);
          end
          % i. removal of line hum (all columns at once via elim_hum)
          if all(isfinite(ap.notchCFreq))
            d(:,1,:)=elim_hum(d(:,1,:),wp.si,ap.notchCFreq);
          end
          % ii. artifact elimination
          if ~isempty(ap.afElimination)
            tmpInputArg=eval(ap.afElimination);
            d(:,1,:)=elim_artefact(d(:,1,:),wp.si,tmpInputArg{:});
          end
          % iii.a lowpass filter
          if isfinite(ap.loCFreq)
            d(:,1,:)=lofi(d(:,1,:),wp.si,ap.loCFreq);
          end
          % iii.b resampling (don't forget tmpD)
          if isfinite(ap.sampFac) && ap.sampFac~=1.0
            disp('** resampling...');
            d=resample(d,wp.si,wp.si*ap.sampFac);
            tmpD=resample(tmpD,wp.si,wp.si*ap.sampFac);
            wp.si=wp.si*ap.sampFac;
            wp.rawNEpisodePt=size(d,1);
          end
          % *** replicate d and do all further preconditioning on the
          % second column
          d=cat(2,d,d);
          % iv. highpass filter
          if isfinite(ap.hiCFreq)
            d(:,2,:)=hifi(d(:,2,:),wp.si,ap.hiCFreq,'rs',wp.hiRs);
          end
          % v. differentiator filter
          if ~isempty(ap.differfi)
            tmpInputArg=eval(ap.differfi);
            d(:,2,:)=differfi(d(:,2,:),wp.si,tmpInputArg{:});
          end
          % vi. custom command
          try
            eval(ap.customCommand);
          catch MExc
            errordlg(['custom command failed. Here''s the error message: ' MExc.message]);
          end
          % if this is episodic data reshape back 
          if ~isnan(wp.rawNEpisode) && wp.rawNEpisode>1
            d=reshape(permute(d,[1 3 2]),[wp.rawNEpisodePt*wp.rawNEpisode 2]);
          end
          % vi. re-concatenate, order: downsampled raw, full precond,
          % partly precond
          d=cat(2,tmpD,d(:,[2 1]));
          % don't forget to recompute length of raw data excerpt in points
          wp.excLen_pts=cont2discrete(wp.excLen,wp.si/1000,'intv',1);
          wp.rawNRow=size(d,1);
          % vii. compute base line characteristics, poor man's version:
          % divide data into short segments of 200 ms and compute noise of
          % each; in histogram of these pick maximum 
          excPts=cont2discrete(200,wp.si/1000,'intv',1);
          % note tmpD being re-used
          tmpD=d(1:end-rem(wp.rawNRow,excPts),2);
          tmpD=reshape(tmpD,excPts,size(tmpD,1)/excPts);
          [~,v]=detbaseline(tmpD,'meth','median');
          % the simple assumption here is that intervals with the lowest
          % noise level correspond to periods in which there is little
          % neuronal activity, so these are the baseline ones we are
          % looking for
          % - identify segments in bin with maximal count
          [n,~,binIx]=histcounts(v,'BinMethod','fd');
          [~,ix]=max(n);
          % - select these segments
          tmpD=tmpD(:,binIx==ix);
          % from these, compute median, iqr/2
          [m,v]=detbaseline(tmpD(:),'meth','median');
          wp.baselinePar=[m v];
          % viii. set relative threshold
          if ~isempty(ap.relativeThresh) && isfinite(ap.relativeThresh)
            ap.thresh=wp.baselinePar(1)+wp.baselinePar(2)*ap.relativeThresh;
          end
          clear tmpD
        else
          h=warndlg('current data have already been preconditioned');
          uiwait(h);
        end
      else
        h=warndlg('no data loaded');
        uiwait(h);
      end
      % next job: re-plotting data
      job(1)={'plotOv'};
      
    case 'detEvt'
      % tsl must be emptied because it may exist from former session
      evt.tsl=[];
      evt.amp=[];
      evt.tRise=[];
      wp.evtTsl=[];
      evtCutout=[];
      if any(ishandle(sp.rawOv.evMarkH))
        delete(sp.rawOv.evMarkH);
        sp.rawOv.evMarkH=nan;
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
      if ~isnan(ap.thresh) && ~isempty(d)
        % § this line is needed in case the threshold had been set in the
        % options menu (as opposed to clicking on the subplot). It should
        % go once this is checked for and the threshold is set
        % automatically, see line ~310)
        set(sp.rawExc.threshH,'ydata',ap.thresh*[1 1]);
        % wipe iei & cutouts
        subplot(sp.cutout.axH);
        th=ultext('se patiencier, svp...',0.02,'fontsize',12,'color','r');
        % drawnow;
        % column of d in which to detect events
        if wp.precondFlag
          detEvCol=2;
        else
          detEvCol=1;
        end
        % time stamps in ticks
        evt.tsl=tcd(d(:,detEvCol),'idx',ap.thresh,'detMode',ap.evtDetMode);
        delete(th);
        % kill all spx following predecessor by less than dead time
        if isfinite(ap.evDeadT) && ~isempty(evt.tsl{1})
          evt.tsl={tsldeadt(evt.tsl{1},cont2discrete(ap.evDeadT,wp.si/1e3,'intv',0))};
        end
        nTs=length(evt.tsl{1});
        tmpwinEvtCutout=cont2discrete(ap.winEvtCutout,wp.si/1e3,'intv',1);
        % column of d from which to extract cutouts
        if wp.precondFlag
          tmpCol=3-double(~ap.doUncondCutout);
        else
          tmpCol=1;
        end
        if nTs>0
          [evtCutout,cti]=tsl2exc(d(:,tmpCol),'idx',evt.tsl,'win',tmpwinEvtCutout);
          if any(~cti{1})
            evt.tsl{1}=evt.tsl{1}(cti{1});
            warning('given current interval for cutouts at least one event could not be retained because it is too close to either border');
          end
          % ** important to convert to cell array
          evtCutout={evtCutout};
        else
          warning('no events detected (threshold is beyond range of data amplitudes)');
        end
        % copy tsl
        wp.evtTsl=evt.tsl;
        % after all is done convert tsl to ms
        evt.tsl{1}=discrete2cont(evt.tsl{1},wp.si/1e3,'intv',0);
      else
        warndlg('no data loaded or threshold not set');
      end
      % next job: plotting markers & histogram
      job(1)={'plotEvt'};
      
    case 'treatEvt'
      % check for contents of tsl and cutouts
      if isempty(evt.tsl) || isempty(evtCutout)
        warndlg('events cannot be treated because none were detected');
        job(1)=[];
      elseif numel(evt.tsl)>1
        % if there is more than one population of events (=units) offer the
        % option to kill one of them (useful to get rid of artifacts)
        answer=inputdlg({'Enter index of event population to purge'},...
          'Event pop organizer',1,{'1'});
        if ~isempty(answer)
          killIx=str2double(answer{1});
          if killIx>0 && killIx<=numel(evt.tsl)
            evt.tsl(killIx)=[];
            if ~isempty(evt.amp)
              evt.amp(killIx)=[];            
              evt.tRise(killIx)=[];            
            end
            evtCutout(killIx)=[];
            wp.evtTsl(killIx)=[];
            job(1)={'plotEvt'};
          else
            warndlg('illegal index entered')
            job(1)=[];
          end
        else
          job(1)=[];
        end
      else
        % call the cutout-treating gui
        tmph=cutdealgui;
        uiwait(tmph);
        job(1)={'plotEvt'};
      end
      
    case {'sepEvt'}
      % check for contents of tsl and cutouts
      if isempty(evt.tsl) || isempty(evtCutout)
        warndlg('events cannot be separated because none were detected');
        job(1)=[];
      elseif numel(evt.tsl)>1
        warndlg('events seem to have been separated already')
        job(1)=[];
      else
        % separation of events (observations) via principal components
        [~,nd]=PCexplore(evtCutout{1}','nPC',4,'omitMd','obs','normalize',0);
        % GUI for separating events - use 'userdata' of cutouts plot for
        % temporary storage of results
        PCexploreGUI(evtCutout{1}',nd,sp.cutout.axH);
        % *** block all further execution of code until unit separation is
        % done ***
        uiwait(findobj('tag','PCexplore'));
        pop=get(sp.cutout.axH,'userdata');
        set(sp.cutout.axH,'userdata',[]);
        nPop=numel(pop);
        if nPop>0
          %  now really split them up
          tmptsl=evt.tsl{1};
          if ~isempty(evt.amp)
            tmpAmp=evt.amp{1};
            tmpTRise=evt.tRise;            
          end
          tmpCutouts=evtCutout{1};
          allTsIx=1:numel(tmptsl);
          % first group of marked events will be the first population,
          % the second will be the second, and so on; the remainder will
          % be assigned to the last pop
          for ui=1:nPop+1
            if ui==nPop+1
              helperIx=setdiff(allTsIx,cat(1,pop.ix));
              evt.tsl{ui}=tmptsl(helperIx);
              evtCutout{ui}=tmpCutouts(:,helperIx);
              if ~isempty(evt.amp)
                evt.amp{ui}=tmpAmp(helperIx);
                evt.tRise{ui}=tmpTRise(helperIx);
              end
            else
              evt.tsl{ui}=tmptsl(pop(ui).ix);
              evtCutout{ui}=tmpCutouts(:,pop(ui).ix);
              if ~isempty(evt.amp)
                evt.amp{ui}=tmpAmp(pop(ui).ix);
                evt.tRise{ui}=tmpTRise(pop(ui).ix);
              end
            end
          end
          clear tmp* nUnit allTsIx
          % make sure all pops contain at least one spike
          emptyIx=[];
          for ui=1:length(evt.tsl)
            if isempty(evt.tsl{ui})
              emptyIx=[emptyIx ui];
            end
          end
          if ~isempty(emptyIx)
            evt.tsl(emptyIx)=[];
            evtCutout(emptyIx)=[];
            if ~isempty(evt.amp)
              evt.amp(emptyIx)=[];
              evt.tRise(emptyIx)=[];
            end
            warndlg(['unit separation resulted in ' int2str(numel(emptyIx)) ' empty time stamp lists (which were deleted)']);
          end
        end
        % convert time stamps of detected events to ticks for markers on
        % excerpt plot (time is ms)
        for ui=1:length(evt.tsl)
          wp.evtTsl{ui}=cont2discrete(evt.tsl{ui},wp.si/1e3,'intv',0);
        end
        % restore original graphics defaults as defined in threshdetgui.m 
        % § could be defined by wp 
        labelscale('fontSz',8,'scaleFac',1.0,'lineW',.25,'markSz',6); 
        % next job: plotting markers & histogram
        job(1)={'plotEvt'};
      end
  
    case 'detEvtAmp'
      % this is a very specialized job, computing amplitudes of PSCs, that
      % is, events with a steep rise and a slower decay, even if they occur
      % in bursts, that is, if they are stacked upon each other. 
      doDetEvtAmp=true;
      % perform a number of checks
      if isempty(evt.tsl) || isempty(evtCutout)
        warndlg('you must perform event detection prior to determination of event amplitudes');
        doDetEvtAmp=false;
      end
      if numel(evt.tsl)>1
        warndlg('event amplitudes cannot be computed after separation into populations')
        doDetEvtAmp=false;
      end
      if ~wp.precondFlag && doDetEvtAmp
        warndlg('data must be preconditioned for computation of event amplitude');
        doDetEvtAmp=false;
      end
      if ~ap.doUncondCutout
        warndlg('cutouts must be produced from *partially* preconditioned trace');
        doDetEvtAmp=false;
      end        
      if (ap.winEvtCutout(1)>=0 || ap.winEvtCutout(2)<=0) && doDetEvtAmp
        warndlg('computation of event amplitude requires a cutout window embracing the time stamp');
        doDetEvtAmp=false;
      end
      if doDetEvtAmp
        % pop up window right away so user knows it's in the make
        tmpFh=findobj('Tag','evAmpFigure','type','figure');
        if isempty(tmpFh)
          tmpFh=figure('Units','normalized', ...
            'Name','EvAmplitude Plot', ...
            'NumberTitle','off', ...
            'Position',[0.005 0.33 0.99 0.35], ...
            'Tag','evAmpFigure'...
            );
        else
          figure(tmpFh);
          cla
        end
        % start by creating cutouts from fully preconditioned column of d
        % in which events were detected (using wp.evtTsl because it's in
        % points)
        detEvCol=2;
        tmpwinEvtCutout=cont2discrete(ap.winEvtCutout,wp.si/1e3,'intv',1);
        tmpCutout=tsl2exc(d(:,detEvCol),'idx',wp.evtTsl,'win',tmpwinEvtCutout);
        % call detPSCAmp - evtCutout must have been produced from partly
        % preconditioned data (see error checking above)
        [evt.amp,evt.tRise,tmpBase]=detPSCAmp(evtCutout{1},tmpCutout,1-tmpwinEvtCutout(1),...
          wp.si,wp.evtTsl{1},'d',d(:,3),'thresh',ap.thresh,'fh',tmpFh,...
          'nPlotEv',min(1000,numel(wp.evtTsl{1})),'plotOverview',true);
        % convert to cell 
        evt.amp={evt.amp};
        evt.tRise={evt.tRise};
        % evt.base={tmpBase};
      end
      job(1)=[];

    case 'experimental'
      
      warndlg('experimental procedures inactive');
      return

      tmpFh=findobj('Tag','ExpAnPlot','type','figure');
      if isempty(tmpFh)
        tmpFh=mkfig([],'b');
        set(tmpFh,'Tag','ExpAnPlot','name','Experimental Analysis Plot');
      else
        figure(tmpFh);
        clf
      end
      
      anCol=2;
      % range of thresholds to test
      p=prctile(d(:,anCol),[70 99.9]);
      nThresh=20;
      threshArr=flipud(linspace(p(1),p(2),nThresh)');

      resArr=nan(nThresh,2);
      waitBarH = waitbar(0,'Computing...');
      disp('computing...')
      for k=1:nThresh
          waitbar(k/nThresh,waitBarH);
          tmpEtsl=tbt(d(:,anCol),...
            'idx',threshArr(k),...
            'minActive',0,'minInactive',100,...
            'elimOrder','inactive');
          activeIx=etsl2logical(wp.rawNRow,tmpEtsl);
%         activeIx=d(:,anCol)>=threshArr(k);
%         % post-ts 
%         for ii=1:50
%           activeIx=activeIx | circshift(activeIx,1);
%         end
%         % pre-ts 
%         for ii=1:50
%           activeIx=activeIx | circshift(activeIx,-1);
%         end
      end
      delete(waitBarH);
      if diff(nanmean(resArr))>3
        yyaxis left
        plot(threshArr,resArr(:,1),'o-');
        ylabel('all')
        yyaxis right
        plot(threshArr,resArr(:,2),'o-');
        ylabel('inactive')
      else
        plot(threshArr,resArr,'o-');
      end
      
      job(1)=[];
      
    case 'plotEvt'
      figure(findobj('Tag','threshdetgui','type','figure'));
        if any(ishandle(sp.rawOv.evMarkH))
          delete(sp.rawOv.evMarkH);
          sp.rawOv.evMarkH=nan;
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
      
      if ~iscell(evt.tsl) || isempty(evt.tsl{1})
        subplot(sp.cutout.axH); hold on
        axis([0 2 0 2]);
        th=text(1,1,':-(');
        set(th,'fontsize',60,'HorizontalAlignment','center','Rotation',-90);
        urtext('0 events');
      else
        % helper vars for plotting markers of spx
        nTs=0;
        sp.cutout.axH.XTickMode='auto';
        sp.cutout.axH.YTickMode='auto';
        % ** watch out: for multi-sweep data we cannot use
        % ds.fileInfo.recTime to know the x extent of the overview plot
        % window because the sweeps were concatenated without gaps even if
        % in real time there were gaps **
        tmpFac=diff(get(sp.rawOv.axH,'xlim'))/(wp.rawNRow/1e6*wp.si);
        tmpyl=get(sp.rawOv.axH,'ylim');
        for ui=1:length(evt.tsl)
          nTs(ui)=length(evt.tsl{ui});
          pCol=wp.tsMarkerCol{mod(ui-1,length(wp.tsMarkerCol))+1};
          % plot markers for detected events in overview plot
          tmpy=repmat(tmpyl(1)+diff(tmpyl)*sp.rawOv.EvMarkRelHeight, [nTs(ui) 1]);
          subplot(sp.rawOv.axH);
          sp.rawOv.evMarkH(ui)=plot(evt.tsl{ui}/1000*tmpFac,tmpy,[pCol '.']);
          % plot specified number of cutouts
          tmpv1=min(wp.maxNPlotCutout,nTs(ui));
          tmpix=unique(ceil((1:tmpv1)/tmpv1*nTs(ui)));
          time=discrete2cont(1:size(evtCutout{ui},1),wp.si/1000,'intv',0)+ap.winEvtCutout(1);
          subplot(sp.cutout.axH); hold on
          plot(time,evtCutout{ui}(:,tmpix),pCol);
        end
        axis tight
        set(sp.cutout.axH,'color',wp.stdAxCol);
        xlabel('time (ms)');
        % place information on firing rates in legend
        tmps=cat(2,num2str(nTs'/diff(ds.fileInfo.recTime),'%3.3f'),...
          repmat(' Hz (', length(evt.tsl),1),...
          int2str(nTs'),...
          repmat(' events)', length(evt.tsl),1)...
          );
        wp.legH=legend(sp.cutout.axH,sp.rawOv.evMarkH,tmps);
        set(wp.legH,'fontsize',11);
        % replot raw excerpt (in order to display markers)
        callbackfn_rawovplot(sp.rawExc.axH,[],d,'x',wp.curPtX);
        % compute ieih
        h=[];
        for ui=1:length(evt.tsl)
          h(1:numel(wp.ieiBin),ui)=tslisi(evt.tsl{ui},'bins',wp.ieiBin);
        end
        % plot iei
        subplot(sp.iei.axH),
        bh=bar(wp.ieiBin,h,1.0,'stacked');
        for ui=1:length(evt.tsl)
          pCol=wp.tsMarkerCol{mod(ui-1,length(wp.tsMarkerCol))+1};
          set(bh(ui),'facecolor',pCol);
        end
        niceyuax;
        set(sp.iei.axH,'color',wp.stdAxCol);
        xlabel('iei (ms)');
        ylabel('N');
        % finally, if in batch mode, collect event number and rate (** note
        % that we're collecting the data in reverse order here, so they
        % have to be inverted before display. Note also that in batch mode
        % we cannot have more than one unit per channel)
        if wp.isBatchMode
          wp.batchEvNumber(numel(wp.batchFileList)+1)=nTs;
          wp.batchEvRate(numel(wp.batchFileList)+1)=nTs/diff(ds.fileInfo.recTime);
        end
      end
      job(1)=[];

    case 'detBurst'
      doBuDet=true;
      if ~isnan(wp.rawNEpisode) && wp.rawNEpisode>1 
        % episodic data - allow burst detection only if there are no gaps
        % between episodes (but let's not worry about gaps of less than a
        % sample point as such values most likely reflect rounding errors
        % anyways)
        tmp=unique(diff(ds.fileInfo.sweepStartInPts));
        if numel(tmp)==1 && isalmost(tmp/ap.sampFac,wp.rawNEpisodePt,.1)
          doBuDet=true;
        else
          doBuDet=false;
        end
      end
      if doBuDet
        % ******************************************************************
        % etsl must be emptied because it may exist from former session
        bu.etsl=zeros(0,etslc.nCol);
        bu.silentEtsl=zeros(0,etslc.nCol);
        wp.buEtsl=zeros(0,etslc.nCol);
        wp.buSilentEtsl=zeros(0,etslc.nCol);
        buCutout=[];
        wp.manualDeletedSilentPerIx=[];
        wp.transIA=zeros(0,1);
        wp.transAI=zeros(0,1);
        if any(ishandle(sp.rawOv.burstLh))
          delete(sp.rawOv.burstLh);
        end
        if any(ishandle(sp.rawExc.burstLh))
          delete(sp.rawExc.burstLh);
        end
        if any(ishandle(sp.rawOv.silentPerLh))
          delete(sp.rawOv.silentPerLh);
        end
        if any(ishandle(sp.rawExc.silentPerLh))
          delete(sp.rawExc.silentPerLh);
        end
        if any(ishandle(sp.rawExc.burstStartMarkH))
          delete(sp.rawExc.burstStartMarkH);
          sp.rawExc.burstStartMarkH=nan;
        end
        % ******************************************************************
        if ~isnan(ap.thresh) && ~isempty(d)
          % wipe iei & cutouts
          subplot(sp.cutout.axH), cla
          th=ultext('se patiencier, svp...',0.02,'fontsize',12,'color','r');
          subplot(sp.iei.axH), cla
          % drawnow;
          % column of d in which to detect bursts
          if wp.precondFlag
            tmpCol=2;
          else
            tmpCol=1;
          end
          % extended time stamp lists in ms (the first etsl describes the
          % bursts, the second the silent periods in between)
          % ** note: the order in which active and inactive periods of less
          % than minimal length are purged matters! Here, we first
          % eliminate gaps between active periods and THEN get rid of too
          % short active events, because otherwise detecting bursts of
          % high-frequency events would be difficult
          [bu.etsl,bu.silentEtsl,~,wp.transIA,wp.transAI]=tbt(d(:,tmpCol),...
            wp.si/1000,ap.thresh*[1; ap.threshFract],...
            'minActive',ap.minEventWidth,'minInactive',ap.maxBurstGapWidth,...
            'elimOrder','inactive','threshLag',ap.threshLag);
          delete(th);
          job(1)={'afterBurst'};
        else
          warndlg('no data loaded or threshold not set');
          job(1)=[];
        end
      else
        errordlg('burst detection in current episodic data not possible because there are gaps or episodes of differing lengths');
        job(1)=[];
      end
      
    case 'afterBurst'
      % think of 'afterburner' - this job computes burst statistics, plots
      % the bursts, displays statistics, etc.. It is called either
      % immediately following 'detBurst' or after bursts were manipulated
      % manually. It is assumed that plausibility checks (attempt to delete
      % nonexisting silent period) have been accomplished before
      % *** statistics
      bu.stats=etslstats(bu.etsl,bu.silentEtsl);
      nTs=size(bu.etsl,1);
      nSp=size(bu.silentEtsl,1);
      % convert event list of detected bursts to ticks for markers on
      % excerpt plot (time is ms) and for cutting out excerpts
      wp.buEtsl=zeros(size(bu.etsl));
      wp.buEtsl(:,[etslc.tsCol etslc.durCol])=cont2discrete(bu.etsl(:,[etslc.tsCol etslc.durCol]),wp.si/1e3,'intv',0);
      % same for silent periods
      wp.buSilentEtsl=zeros(size(bu.silentEtsl));
      wp.buSilentEtsl(:,[etslc.tsCol etslc.durCol])=cont2discrete(bu.silentEtsl(:,[etslc.tsCol etslc.durCol]),wp.si/1e3,'intv',0);
      if nTs>0
        % column of d from which to extract cutouts
        if wp.precondFlag
          tmpCol=3-double(~ap.doUncondCutout);
        else
          tmpCol=1;
        end
        % ** place variably-sized cutouts in cell array
        tmpwinBuCutout=cont2discrete(ap.winBuCutout,wp.si/1e3,'intv',1);
        [buCutout,isCutout]=etsl2exc(d(:,tmpCol),'idx',wp.buEtsl,'win',tmpwinBuCutout);
        % if no cutout survives etsl2exc will issue a warning, but we
        % need a warndlg nonetheless
        if isempty(isCutout)
          warndlg('with the current cutout window extension no single event survived (too close to either border)')
        elseif any(~isCutout)
          bu.etsl=bu.etsl(isCutout,:);
          wp.buEtsl=wp.buEtsl(isCutout,:);
          % §§§ what about adjacent gaps? they should be kicked out, too
          warning('given current interval for variable-length cutouts at least one event could not be retained because it is too close to either border');
          % §§§ wp.buSilentEtsl should be changed here, too
        end
      else
        if nSp
          warndlg('no event detected/survived - ''minimal event width'' is too short or current threshold is beyond the range of data amplitudes');
        else
          warndlg('no event and no silent period detected/survived - ''minimal event width'' is too short or current threshold is beyond the range of data amplitudes');          
        end
      end
      nTs=size(bu.etsl,1);
      if nTs>0
        % ----- compute & display a few basic parameters:
        % - burst length (median & quartiles)
        bLenQ=prctile(bu.etsl(:,etslc.durCol),[25 50 75]);
        % - length of silent states (median & quartiles)
        spLenQ=prctile(bu.silentEtsl(:,etslc.durCol),[25 50 75]);
        % mis-use cutouts subplot for display of results..
        subplot(sp.iei.axH), cla
        set(gca,'xlim',[0 1],'ylim',[0 1]);
        txt={...
          ['burst rate: ' num2str(bu.stats.burstRate,'%1.4f') ' Hz'],...
          ['burst length: ' num2str(bLenQ(2),'%5.0f') ' [' num2str(bLenQ([1 3]),'%7.0f') '] ms'],...
          ['rel. active time: ' num2str(bu.stats.relTimeInBurst,'%1.4f')],...
          ['silent period length: ' num2str(spLenQ(2),'%5.0f') ' [' num2str(spLenQ([1 3]),'%7.0f') '] ms'],...
          };
        th=text(0.05,0.5,txt,'fontsize',12,'color',[.1 .6 .1]);
        % also, dump numbers on screen, separated by tabs, so they can be
        % copied and pasted to e.g. spreadsheets
        disp(num2str([bu.stats.burstRate bLenQ(2) bLenQ([1 3]) bu.stats.relTimeInBurst spLenQ(2) spLenQ([1 3])],'%1.4f\t %5.2f\t %5.2f\t %1.4f\t %5.2f\t %5.2f'));
        % draw lines marking detected bursts and silent periods in
        % overview plot:
        tmpFac=diff(get(sp.rawOv.axH,'xlim'))/diff(ds.fileInfo.recTime);
        tmpy=get(sp.rawOv.axH,'ylim');
        subplot(sp.rawOv.axH);
        % - bursts
        tmpyExt=repmat(tmpy(1)+diff(tmpy)*sp.rawOv.BurstLineRelHeight, [2 size(bu.etsl,1)]);
        sp.rawOv.burstLh=line(cumsum((bu.etsl(:,[etslc.tsCol etslc.durCol]))'/1000,1)*tmpFac,tmpyExt,'color','r','linewidth',2);
        % silent periods:
        tmpyExt=repmat(tmpy(1)+diff(tmpy)*sp.rawOv.BurstLineRelHeight, [2 size(bu.silentEtsl,1)]);
        sp.rawOv.silentPerLh=line(cumsum((bu.silentEtsl(:,[etslc.tsCol etslc.durCol]))'/1000,1)*tmpFac,tmpyExt,'color','b','linewidth',2);
        % replot excerpt (in order to display markers)
        callbackfn_rawovplot(sp.rawExc.axH,[],d,'x',wp.curPtX);
        % finally, plot specified number of cutouts (first 500 ms including pretrigger win)
        tmpv1=min(wp.maxNPlotCutout,nTs);
        tmpBuIx=unique(ceil((1:tmpv1)/tmpv1*nTs));
        time=discrete2cont(1:500/(wp.si/1000),wp.si/1000,'intv',0)+ap.winBuCutout(1);
        maxNPt=numel(time);
        tmpBu=nan(maxNPt,numel(tmpBuIx));
        for ix=1:numel(tmpBuIx)
          nPt=numel(buCutout{tmpBuIx(ix)});
          nPt=min(maxNPt,nPt);
          tmpBu(1:nPt,ix)=buCutout{tmpBuIx(ix)}(1:nPt);
        end
        subplot(sp.cutout.axH);
        cla;
        plot(time,tmpBu,'k');
        axis tight
        set(sp.cutout.axH,'color',wp.stdAxCol);
        xlabel('time (ms)');
      end
      job(1)=[];

    case 'saveResults'
      % first thing to do: build results file names
      tmpix=strfind(ds.dataFn,'.');
      % names of .mat files
      ap.resFn=[ds.dataPath ds.dataFn(1:tmpix(end)-1) '_'  wp.dataChanName{1}(~isspace(wp.dataChanName{1})) '_' wp.resFnString '_res'];
      ap.rawFn=[ds.dataPath ds.dataFn(1:tmpix(end)-1) '_'  wp.dataChanName{1}(~isspace(wp.dataChanName{1})) '_' wp.resFnString '_raw'];
      ap.evtCutoutFn=[ds.dataPath ds.dataFn(1:tmpix(end)-1) '_'  wp.dataChanName{1}(~isspace(wp.dataChanName{1})) '_' wp.resFnString '_evtCutout'];
      ap.buCutoutFn=[ds.dataPath ds.dataFn(1:tmpix(end)-1) '_'  wp.dataChanName{1}(~isspace(wp.dataChanName{1})) '_' wp.resFnString '_buCutout'];
      % name of ascii files
      ap.asciiEvtFn=[ds.dataPath ds.dataFn(1:tmpix(end)-1) '_'  wp.dataChanName{1}(~isspace(wp.dataChanName{1})) '_' wp.resFnString];
      % dump ap, ds, and wp into results file so parameters can be
      % retrieved post-analysis
      head.ds=ds;
      head.ap=ap;
      head.wp=wp;
      % ** important: remove all graphics handles from variables (wp)
      % because saving graphics handles in Matlab versions 2014b and later
      % implies writing the whole graphics object to file!
      % - first get rid of the major culprit, a struct of handles
      head.wp=rmfield(head.wp,'mainGuiHandles');
      % - eliminate others
      fiNa=fieldnames(head.wp);
      fiKiIx=false(size(fiNa));
      for fiInd=1:numel(fiNa)
        if any(isgraphics(head.wp.(fiNa{fiInd})(:) )) || any(ishandle(head.wp.(fiNa{fiInd})(:)))
          fiKiIx(fiInd)=true;
        end
      end
      head.wp=rmfield(head.wp,fiNa(fiKiIx));
      % ------------------------------------------------------------------
      % ** if more than one unit was detected and separated write both
      % time stamps and cutouts in separate files. The file names will be
      % modified accordingly. This may appear cumbersome, but gives the
      % user better control over post-event detection analyses
      if wp.saveEvt
        % *** first thing: copy evt into temporary var (because evt may be
        % 'formatted' below, and if then the data are accidentally saved
        % again an error may occur or the reformatting may mess evt up)
        tmpEvt0=evt;
        % a case in point: take care of episodic data by reformatting evt
        if ~isnan(wp.rawNEpisode) && wp.rawNEpisode>1
          % if data is episodic append field 'sweepIx' to evt which
          % contains the index to the sweep to which the time stamp belongs
          % and transform the time stamps back to their local time frame
          for ui=1:length(evt.tsl)
            evt.tsl{ui}=cont2discrete(evt.tsl{ui},wp.si/1e3,'intv',0);
            evt.sweepIx{ui}=floor((evt.tsl{ui}-1)/wp.rawNEpisodePt)+1;
            evt.tsl{ui}=discrete2cont(mod(evt.tsl{ui}-1,wp.rawNEpisodePt)+1,wp.si/1e3,'intv',0);
          end
          % §§ the same would not make an awful lot of sense for burst
          % etsls because burst have nonzero length and could thus straddle
          % sweeps
        end
        % now copy the present version of evt into another temporary var
        tmpEvt=evt;
        % if field .amp is empty, get rid of it
        if isempty(evt.amp)
          evt=rmfield(evt,'amp');
          evt=rmfield(evt,'tRise');
        end
        tmpResFn=[];
        if length(evt.tsl)>1
          % for each unit intercalate a '.1', '.2' etc. and place result in
          % cell array
          sstring=wp.dataChanName{1}(~isspace(wp.dataChanName{1}));
          for ui=1:length(evt.tsl)
            tmpResFn{ui}=strrep(ap.resFn,sstring,[sstring '.' int2str(ui)]);
          end
        else
          % use results file name as is, but wrap up nicely in a cell 
          tmpResFn={ap.resFn};
        end
        evtFieldName=fieldnames(evt);
        for ui=1:length(evt.tsl)
          for h=1:numel(evtFieldName)
            % separate elements of fields of evt
            evt.(evtFieldName{h})=tmpEvt.(evtFieldName{h})(ui);
          end
          if exist([tmpResFn{ui} '.mat'],'file')
            % overwrite everything except bu
            save([tmpResFn{ui} '.mat'],'evt','head','-mat','-append');
          else
            % save everything
            save([tmpResFn{ui} '.mat'],'evt','bu','head','-mat');
          end
        end
        % revert evt to original post-detection format
        evt=tmpEvt0;
        clear tmp*
      end
      % ------------------------------------------------------------------      
      if wp.saveBu
        % burst data as generated here are by definition from one 'unit',
        % so there is no need to write results in different files
        if exist([ap.resFn '.mat'],'file')
          % overwrite everything except evt
          save([ap.resFn '.mat'],'bu','head','-mat','-append');
        else
          % save everything 
          save([ap.resFn '.mat'],'evt','bu','head','-mat');
        end
      end
      if wp.saveEvtCutout
        % ** saving event cutouts: same logic as with tsl applies
        tmpCutoutFn=[];
        if length(evtCutout)>1
          % for each unit intercalate a '.1', '.2' etc. and place result in
          % cell array
          sstring=wp.dataChanName{1}(~isspace(wp.dataChanName{1}));
          for ui=1:length(evtCutout)
            tmpCutoutFn{ui}=strrep(ap.evtCutoutFn,sstring,[sstring '.' int2str(ui)]);
          end
        else
          % use cutouts file name as is, but wrap up nicely in a cell
          tmpCutoutFn={ap.evtCutoutFn};
        end
        % now copy cutouts into temporary var
        tmpEvtCutout=evtCutout;
        for ui=1:length(evtCutout)
          evtCutout=tmpEvtCutout(ui);
          save([tmpCutoutFn{ui} '.mat'],'evtCutout','head','-mat');
        end
        % reverse
        evtCutout=tmpEvtCutout;
        clear tmp*
      end
      % ------------------------------------------------------------------
      if wp.saveBuCutout
        save([ap.buCutoutFn '.mat'],'buCutout','head','-mat');
      end
      % ------------------------------------------------------------------
      if wp.savePrecondData
        save([ap.rawFn '.mat'],'d','head','-mat');
      end
      disp('**** data saved');
      % finally, if we're in batch mode and if wp.batchFileList is empty
      % set kill batch job
      if wp.isBatchMode && isempty(wp.batchFileList)
        job{1}='killBatchJob';
      else
        job(1)=[];
      end

    case 'done'
      disp('bye...');
      clear global
      job(1)=[];

    otherwise
      error(['internal:illegal job:' partJob]);
      
  end
  drawnow
  jobsToDo= ~isempty(job);
end


