importTool('optimi')
importTool('kern')
importTool('noise')
importTool('ndlutil')
cd ~/mlprojects/ivm/matlab 
ivmRunDataSetRegression('pumadynSeeger', 1, {'rbfard', 'white', 'bias'}, 'gaussian', 'rentropy', 50, 1e5);
exit
