# MACKtrack
This is a MATLAB toolbox designed for automated cell tracking, analysis, and visualization from sequences of live-cell microscopy images.
Author: Brooks Taylor, b4taylor@ucla.edu

## Download & Installation
From command line:

```
$ git clone https://github.com/brookstaylorjr/MACKtrack.git
...
$ cd MACKtrack
$ sh install.sh
...
```

## Basic usage using the GUI
1. Ensure that MACKtrack (and all subfolders) have been added to your MATLAB path.
2. To launch, type `MACKtrack` at the MATLAB command prompt.
3. A window will open, allowing you to personalize locations for your computer: in order to allow tracking on a computational server, then later analysis of the results on your personal machine, MACKtrack creates composite paths for both the input microscopy images and the output data: you choose the "mount" location for (e.g.) a NAS containing images or data, and then specify a subdirectory under that location. If images or data are stored locally, just specify your system's file root (`/` or `C:\`).
4. Set parameters, then test and initiate tracking/measurement from the MACKtrack GUI. Measured trajectories will be saved to `AllMeasurements.mat`
5. In order to explore results, a general-purpose tool is included. Launch this by typing `MACKquery` at the MATLAB command prompt, then navigate to, and load, your newly-saved `AllMeasurements.mat` file.

## Analyzing experiments from a spreadsheet
1. Establish a set of parameters that works for the set(s) to be tracked using the GUI, and save in the `Parameters` folder.
2. Publish a Google Spreadsheet like the one [here](https://docs.google.com/spreadsheets/d/10o_d9HN8dhw8bX4tbGxFBJ63ju7tODVImZWNrnewmwY/edit#gid=0) and put the link in the `locations.mat` file. The essential column headings are in blue (the spreadsheet is parsed by these headings, so yours need to match mine exactly).
3. Run the command (e.g.) `runID(1:4)` to track and measure the first 4 entries of your spreadsheet.


## TODO:
* [ ] Make package into a matlab toolbox (e.g.: potterswheel or https://github.com/nwh/matlab)
  + During install, move all the files into /path-to-matlab/toolbox/MACKtrack; make that discoverable during startup 
  + It then can be used by any users in the same computer or cluster 
* [ ] Make a wiki or MACKtrack.github.io pages for better documentation. 
  + eg: https://select2.github.io/
