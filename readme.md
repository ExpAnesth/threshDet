# threshDet

is a graphical user interface for the detection of neuronal signals in time series. Signals to be detected can be either events of fixed duration (e.g. action potentials) or bursts (of variable duration, e.g. Up States). The code makes use of in-house formats for time information, namely 'time stamp lists' (tsl) and 'extended time stamp lists' (etsl). 

![screenshot](/doc/screenshot.png)

### Features: 
* reading of files in diverse formats (presently, Axon Binary Format (abf), HDF5 as resulting from conversion of MultiChannelSystems (*.mcd) files, custom Matlab format)
* multiple ways of pre-processing the time series (removal of artifacts; highpass, lowpass, notch and differentiator filters; up- and downsampling; interface for custom commands)
* threshold-based algorithms of event and burst detection
* inspection and removal of individual events
* sorting of events via PCA
* saving to disk of time stamps, excerpts of detected events/bursts, and pre-processed time series

See **manual_threshdetgui.pdf** for an intro with illustrations.

Please note that the code in this repository is not self-sufficient, you'll additionally need the following repositories:
* fileIO
* etslfunc
* graphics
* sampledSeries
* utilities


## General note on repositories in the ExpAnesth organization
Except where noted, code was written by Harald Hentschke, Section of Experimental Anesthesiology, Department of Anesthesiology, University Hospital of Tuebingen. It has been designed primarily for in-house use by individuals who were instructed on its purpose and limitations. Also, a substantial proportion of the code has been developed over a time span of >10 years. Therefore,

* documentation is anywhere from non-existent to incomplete
* due to its gradual development over years, design and implementation do not necessarily reflect best programming practice and techniques. For example, in terms of design, many code files are quite large and had better be broken down into smaller units. In terms of techniques, e.g. the new automatic array expansion as introduced in Matlab Release 2016b features only rarely. The code will be improved and updated when and where the need arises.