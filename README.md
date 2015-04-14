# UCSDcellTrack
This is a MATLAB toolbox designed for automated cell tracking, analysis, and visualization from sequences live-cell microscopy images.
Author: Brooks Taylor, b4taylor@ucla.edu

## Download & Installation
From command line:

```
$ git clone https://github.com/brookstaylorjr/UCSDcellTrack.git
...
$ cd UCSDcellTrack
$ sh install.sh
...
```

## Basic usage using the GUI
1. To launch, type `UCSDcellTrack` at the MATLAB command prompt.
2. A window will open, allowing you to personalize locations for your computer: in order to allow tracking on a computational server, then later analysis of the results on another machine, UCSDcellTrack creates composite paths for both the input microscopy images and the output data: you choose the "mount" location for (e.g.) a NAS containing images or data, and then specify a subdirectory under that location. If images or data are stored locally, just specify your file root (`/` or `C:\`).
3. Set parameters, then test and initiate tracking/measurement from the UCSDcellTrack GUI. Measured trajectories will be saved to `AllMeasurements.mat`
4. In order to explore results, a general-purpose tool is included. Launch this by typing `UCSDcellQuery` at the MATLAB command prompt, then navigate to, and load, your newly-saved `AllMeasurements.mat` file.

## Analyzing experiments from a spreadsheet
1. Establish a set of parameters that works for the set(s) to be tracked using the GUI, and save in the `Parameters` folder.
2. Publish a Google Spreadsheet like the one at `https://docs.google.com/spreadsheets/d/10o_d9HN8dhw8bX4tbGxFBJ63ju7tODVImZWNrnewmwY/edit#gid=0` and put the link in the `locations.mat` file. The essential column headings are in blue (the spreadsheet is parsed by these headings, so yours need to match mine exactly).
3. Run the command (e.g.) `runID(1:4)` to track and measure the first 4 entries of your spreadsheet.

## Included scripts:

## Algorithms


## TODO:
* Make it a matlab toolbox, ref: potterswheel or https://github.com/nwh/matlab
