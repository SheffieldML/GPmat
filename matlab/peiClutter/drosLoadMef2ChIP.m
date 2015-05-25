% DROSLOADMEF2CHIP Load Mef2 Chip Data.


DATAPATH = sprintf('./data/eileen_data/');

chipdata = importdata([DATAPATH, 'chip2/mef2_chipdata_fbgn.txt']);
drosmef2chip.labels = {'mef2_2-4', 'mef2_4-6', 'mef2_6-8', 'mef2_8-10', 'mef2_10-12'};
drosmef2chip.data = chipdata.data;
drosmef2chip.genes = chipdata.rowheaders;

clear chipdata
