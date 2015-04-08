importTool('optimi')
importTool('kern')
importTool('noise')
importTool('ndlutil')
cd ~/mlprojects/ivm/matlab 
vivmRunDataSet('usps0', 'translate', 1, 'rentropy', 1000, 1e5);
exit
