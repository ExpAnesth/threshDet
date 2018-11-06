# threshDet

is a graphical user interface for the detection of neuronal signals in time series. Signals to be detected can be either events of fixed duration (e.g. action potentials) or bursts (of variable duration, e.g. Up States). The code makes use of in-house formats for time information, namely 'time stamp lists' (tsl) and 'extended time stamp lists' (etsl). 

![screenshot](/doc/screenshot_threshDet.png)

### Features: 
* reading of files in diverse formats (presently, Axon Binary Format (abf), HDF5 as resulting from conversion of MultiChannelSystems (*.mcd) files, custom Matlab format)
* multiple ways of pre-processing the time series (removal of artifacts; highpass, lowpass, notch and differentiator filters; up- and downsampling; interface for custom commands)
* threshold-based algorithms of event and burst detection
* optional automatic adjustment of threshold level to baseline noise
* inspection and removal of individual events
* sorting of events via PCA
* batch mode of event detection
* saving to disk of time stamps, excerpts of detected events/bursts, and pre-processed time series

See **manual_threshdetgui.pdf** for an intro with illustrations.

Please note that the code in this repository is not self-sufficient, you'll additionally need the following repositories:
* fileIO
* etslfunc
* graphics
* sampledSeries
* utilities

Matlab toolboxes required:
* Signal Processing
* Statistics and Machine Learning


## General note on repositories in the ExpAnesth organization
The code in these repositories provides basic tools for the analysis of electrophysiological time series to members of the Section of Experimental Anesthesiology, Department of Anesthesiology, University Hospital of Tuebingen. Except where noted, code was written by Harald Hentschke. It has been designed primarily for in-house use by individuals who were instructed on its scope and limitations. Also, a substantial proportion of the code has been developed and extended over a time span of >10 years. In detail,

* the implementation of algorithms reflects the evolution of Matlab itself, that is, code that had been initially developed on older versions of Matlab does not necessarily feature newer techniques such as the new automatic array expansion as introduced in Matlab Release 2016b
* nonetheless, all code has been tested to run on Matlab R2018b
* while most m-files contain ample comments, documentation exists only for a few repositories
* checks of user input are implemented to varying degrees
* the code will be improved, updated and documented when and where the need arises