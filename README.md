# UCSDcellTrack
This is the UCSD cell tracking and quantification software. 

## Download & Installation
From command line:

```
$ git clone https://github.com/brookstaylorjr/UCSDcellTrack.git
...
$ cd UCSDcellTrack
$ sh install.sh
...
```

## Usage
1. Establish a set of parameters that works for the set(s) to be tracked using the GUI, and save in the `Parameters` folder.
2. Publish a Google Spreadsheet like mine and put the link in the `locations.mat` file. The essential column headings are in blue (the spreadsheet is parsed by these headings, so yours need to match mine exactly).
3. Run the command (e.g.) `runID(1:4)` to track and measure the first 4 entries of your spreadsheet.

## Included scripts:

## Algorithms


## TODO:
* Make it a matlab toolbox, ref: potterswheel or https://github.com/nwh/matlab