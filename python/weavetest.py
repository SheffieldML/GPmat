import numpy as np
from scipy.weave import inline
	
def testfn(x1,x2):
	code = """
	double *p1, *p2, *p3, *p4;
	PyObject *itr;
	itr = PyArray_MultiIterNew(4, x1_array, x2_array,x3_array,x4_array);
	while(PyArray_MultiIter_NOTDONE(itr)) {
		p1 = (double *) PyArray_MultiIter_DATA(itr, 0);
		p2 = (double *) PyArray_MultiIter_DATA(itr, 1);
		p3 = (double *) PyArray_MultiIter_DATA(itr, 2);
		p4 = (double *) PyArray_MultiIter_DATA(itr, 3);

		*p4 = 1;
		*p3 = 1;

		PyArray_MultiIter_NEXT(itr);
	}
	Py_DECREF(itr);
	"""
	x1 = np.atleast_1d(x1)
	x2 = np.atleast_1d(x2)
	bcs = np.broadcast(x1,x2).shape
	x3 = np.empty(bcs,dtype=np.float64)
	x4 = np.empty(bcs,dtype=np.float64)
	inline(code,['x1','x2','x3','x4'])
	return (x3,x4)

if __name__=="__main__":
	for i in range(100000):
		a = testfn(np.random.randn(1),np.random.randn(1))
