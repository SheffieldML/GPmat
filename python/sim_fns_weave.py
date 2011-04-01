import numpy as np
import pylab as pb
import pdb
from scipy.special import erf, erfc
from erfcx import erfcx
from scipy.weave import inline

c_erfcx = """
double erfcx(double x){
	double t, u, y;
	t = 3.97886080735226 / ( fabs(x)+ 3.97886080735226);
	u = t-0.5;
	y = (((((((((u * 0.00127109764952614092 + 1.19314022838340944e-4) * u 
	    - 0.003963850973605135)   * u - 8.70779635317295828e-4) * u +     
	      0.00773672528313526668) * u + 0.00383335126264887303) * u -     
	      0.0127223813782122755)  * u - 0.0133823644533460069)  * u +     
	      0.0161315329733252248)  * u + 0.0390976845588484035)  * u +     
	      0.00249367200053503304;
	y = ((((((((((((y * u - 0.0838864557023001992) * u -           
	      0.119463959964325415) * u + 0.0166207924969367356) * u + 
	      0.357524274449531043) * u + 0.805276408752910567)  * u + 
	      1.18902982909273333)  * u + 1.37040217682338167)   * u + 
	      1.31314653831023098)  * u + 1.07925515155856677)   * u + 
	      0.774368199119538609) * u + 0.490165080585318424)  * u + 
	      0.275374741597376782) * t;
	
	if (x<0){
		y = 2.0*exp(x*x) - y;
		}
	return y;
	}
"""
c_lndifferfs = """
double lndifferfs(double *x1, double *x2){
	if (((*x1)*(*x2))<0.0){
	    return log(erf(*x1)-erf(*x2));
	    }
	else if (((*x1)>0.0) & ((*x2)>0.0)){
	    return log(erfcx(*x2) - erfcx(*x1) * exp((*x2)*(*x2) - (*x1)*(*x1))) - (*x2)*(*x2);
	    }
	else if ((*x1)==(*x2)){
	    return -1e10000;
	    }
	else{
	    return log(erfcx(-*x1) - erfcx(-*x2) * exp((*x1)*(*x1) - (*x2)*(*x2))) - (*x1)*(*x1);
	    }
	}
"""

def erfcx(x):
	x = float(x)
	if x > 2.53e307:
		return np.inf
	elif x < -26.228:
		return 0.

	support_code = """
	#include <math.h>;
	double erfcx(double x){
	double t, u, y;
	t = 3.97886080735226 / ( fabs(x)+ 3.97886080735226);
	u = t-0.5;
	y = (((((((((u * 0.00127109764952614092 + 1.19314022838340944e-4) * u 
	    - 0.003963850973605135)   * u - 8.70779635317295828e-4) * u +     
	      0.00773672528313526668) * u + 0.00383335126264887303) * u -     
	      0.0127223813782122755)  * u - 0.0133823644533460069)  * u +     
	      0.0161315329733252248)  * u + 0.0390976845588484035)  * u +     
	      0.00249367200053503304;
	y = ((((((((((((y * u - 0.0838864557023001992) * u -           
	      0.119463959964325415) * u + 0.0166207924969367356) * u + 
	      0.357524274449531043) * u + 0.805276408752910567)  * u + 
	      1.18902982909273333)  * u + 1.37040217682338167)   * u + 
	      1.31314653831023098)  * u + 1.07925515155856677)   * u + 
	      0.774368199119538609) * u + 0.490165080585318424)  * u + 
	      0.275374741597376782) * t;
	
	if (x<0){
		y = 2.0*exp(x*x) - y;
	}

	return y;
	}
	"""
	code = """
	return_val = erfcx(x);
	"""
	return inline(code,['x'],support_code=support_code)

	
def lnDiffErfs(x1,x2):
	"""this only works for x1.shape==x2.shape for the moment"""
	support_code = """
	#include <math.h>;
	"""+c_erfcx+c_lndifferfs
	code = """
	double *p1, *p2, *p3, *p4;
	PyObject *itr;
	itr = PyArray_MultiIterNew(4, x1_array, x2_array,x3_array,x4_array);
	while(PyArray_MultiIter_NOTDONE(itr)) {
		p1 = (double *) PyArray_MultiIter_DATA(itr, 0);
		p2 = (double *) PyArray_MultiIter_DATA(itr, 1);
		p3 = (double *) PyArray_MultiIter_DATA(itr, 2);
		p4 = (double *) PyArray_MultiIter_DATA(itr, 3);

		if((*p1-*p2)<0){
		*p4 =-1;
		*p3 = lndifferfs(p2, p1);
		}
		else{
		*p4 = 1;
		*p3 = lndifferfs(p1, p2);
		}

		PyArray_MultiIter_NEXT(itr);

	}
	Py_DECREF(itr);
	"""
	x1 = np.atleast_1d(x1)
	x2 = np.atleast_1d(x2)
	bcs = np.broadcast(x1,x2).shape
	x3 = np.empty(bcs,dtype=np.float64)
	x4 = np.empty(bcs,dtype=np.float64)
	inline(code,['x1','x2','x3','x4'],support_code=support_code)
	return (x3,x4)

def simComputeH(t1, t2, D_i, D_j, delta_i, delta_j, sigma, compute_gradDdecay1=False, compute_gradDdecay2=False, compute_gradsigma=False):
	"""
	Creates structures ant then uses the in-place version.

	See also docsting for non weave version
	"""

	assert(((t1.shape[1]==1) and (t2.shape[1]==1)),'Input can only have one column')

	dim1 = t1.shape[0];
	dim2 = t2.shape[0];
	t1 = t1 - delta_i;
	t2 = t2 - delta_j;

	t1Mat = t1*np.ones((1,dim2))
	t2Mat = t2*np.ones((1,dim1))
	t2Mat = t2Mat.T
	h = np.empty(t1Mat.shape)
	simComputeH_inplace(t1Mat,t2Mat,D_i,D_j,sigma,h)
	return h



def simComputeH_inplace(t1mat,t2mat,Di,Dj,sigma,H):
	support_code = "#include <math.h>\n"+c_erfcx+c_lndifferfs
	code= """
	double *pt1, *pt2, *pH;
	double halfsigmaDi = 0.5*sigma*Di;
	double halfsigmaDi2 = halfsigmaDi*halfsigmaDi;
	double logDiDj = log(Di+Dj);
	double arg1, arg2;
	PyObject *itr;
	itr = PyArray_MultiIterNew(3, t1mat_array, t2mat_array,H_array);
	while(PyArray_MultiIter_NOTDONE(itr)) {
		pt1 = (double *) PyArray_MultiIter_DATA(itr, 0);
		pt2 = (double *) PyArray_MultiIter_DATA(itr, 1);
		pH = (double *) PyArray_MultiIter_DATA(itr, 2);

		//do the first part
		arg1 = halfsigmaDi+*pt2/sigma;
		arg2 = halfsigmaDi-(1.0/sigma)*(*pt1-*pt2);
		if((arg1-arg2)<0){
			*pH = -1*exp(halfsigmaDi2 - Di*(*pt1-*pt2) +lndifferfs(&arg2,&arg1) - logDiDj);
		}
		else{
			*pH = exp(halfsigmaDi2 - Di*(*pt1-*pt2) +lndifferfs(&arg1,&arg2) - logDiDj);
		}
		
		//do the second part
		arg1 = halfsigmaDi;
		arg2 = halfsigmaDi - *pt1/sigma;
		if((arg1-arg2)<0){
			*pH += exp(halfsigmaDi2 - Di*(*pt1) -Dj*(*pt2) +lndifferfs(&arg2,&arg1) -logDiDj);
		}
		else{
			*pH += -exp(halfsigmaDi2 - Di*(*pt1) -Dj*(*pt2) +lndifferfs(&arg1,&arg2) -logDiDj);
		}
		PyArray_MultiIter_NEXT(itr);

	}
	Py_DECREF(itr);
	"""
	Di = float(Di)
	Dj = float(Dj)
	sigma = float(sigma)
	inline(code,['t1mat','t2mat','Di','Dj','sigma','H'],support_code = support_code)
	

