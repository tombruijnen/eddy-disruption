
#include "mex.h" 
#include <stdio.h>
#include <math.h>

#define GAMMA  26753.0
#define TWOPI	6.283185

/*
#define DEBUG
*/

void multmatvec(double *mat, double *vec, double *matvec)

	/* Multiply 3x3 matrix by 3x1 vector. */

{
*matvec++ = mat[0]*vec[0] + mat[3]*vec[1] + mat[6]*vec[2];
*matvec++ = mat[1]*vec[0] + mat[4]*vec[1] + mat[7]*vec[2];
*matvec++ = mat[2]*vec[0] + mat[5]*vec[1] + mat[8]*vec[2];
}



void addvecs(double *vec1, double *vec2, double *vecsum)

	/* Add two 3x1 Vectors */

{
*vecsum++ = *vec1++ + *vec2++;
*vecsum++ = *vec1++ + *vec2++;
*vecsum++ = *vec1++ + *vec2++;
}




void adjmat(double *mat, double *adj)

/* ======== Adjoint of a 3x3 matrix ========= */

{
*adj++ = (mat[4]*mat[8]-mat[7]*mat[5]);	
*adj++ =-(mat[1]*mat[8]-mat[7]*mat[2]);
*adj++ = (mat[1]*mat[5]-mat[4]*mat[2]);
*adj++ =-(mat[3]*mat[8]-mat[6]*mat[5]);
*adj++ = (mat[0]*mat[8]-mat[6]*mat[2]);
*adj++ =-(mat[0]*mat[5]-mat[3]*mat[2]);
*adj++ = (mat[3]*mat[7]-mat[6]*mat[4]);
*adj++ =-(mat[0]*mat[7]-mat[6]*mat[1]);
*adj++ = (mat[0]*mat[4]-mat[3]*mat[1]);
}


void zeromat(double *mat)

/* ====== Set a 3x3 matrix to all zeros	======= */

{
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
*mat++=0;
}


void eyemat(double *mat)

/* ======== Return 3x3 Identity Matrix  ========= */

{
zeromat(mat);
mat[0]=1;
mat[4]=1;
mat[8]=1;

}

double detmat(double *mat)

/* ======== Determinant of a 3x3 matrix ======== */

{
double det;

det = mat[0]*mat[4]*mat[8];
det+= mat[3]*mat[7]*mat[2];
det+= mat[6]*mat[1]*mat[5];
det-= mat[0]*mat[7]*mat[5];
det-= mat[3]*mat[1]*mat[8];
det-= mat[6]*mat[4]*mat[2];

return det;
}


void scalemat(double *mat, double scalar)

/* ======== multiply a matrix by a scalar ========= */

{
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
*mat++ *= scalar;
}


void invmat(double *mat, double *imat)

/* ======== Inverse of a 3x3 matrix ========= */
/*	DO NOT MAKE THE OUTPUT THE SAME AS ONE OF THE INPUTS!! */

{
double det;
int count;

det = detmat(mat);	/* Determinant */
adjmat(mat, imat);	/* Adjoint */

for (count=0; count<9; count++)
	*imat = *imat++ / det;		
}


void addmats(double *mat1, double *mat2, double *matsum)

/* ====== Add two 3x3 matrices.	====== */

{
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
*matsum++ = *mat1++ + *mat2++;
}


double multmats(double *mat1, double *mat2, double *matproduct)

/* ======= Multiply two 3x3 matrices. ====== */
/*	DO NOT MAKE THE OUTPUT THE SAME AS ONE OF THE INPUTS!! */

{
*matproduct++ = mat1[0]*mat2[0] + mat1[3]*mat2[1] + mat1[6]*mat2[2];
*matproduct++ = mat1[1]*mat2[0] + mat1[4]*mat2[1] + mat1[7]*mat2[2];
*matproduct++ = mat1[2]*mat2[0] + mat1[5]*mat2[1] + mat1[8]*mat2[2];
*matproduct++ = mat1[0]*mat2[3] + mat1[3]*mat2[4] + mat1[6]*mat2[5];
*matproduct++ = mat1[1]*mat2[3] + mat1[4]*mat2[4] + mat1[7]*mat2[5];
*matproduct++ = mat1[2]*mat2[3] + mat1[5]*mat2[4] + mat1[8]*mat2[5];
*matproduct++ = mat1[0]*mat2[6] + mat1[3]*mat2[7] + mat1[6]*mat2[8];
*matproduct++ = mat1[1]*mat2[6] + mat1[4]*mat2[7] + mat1[7]*mat2[8];
*matproduct++ = mat1[2]*mat2[6] + mat1[5]*mat2[7] + mat1[8]*mat2[8];
}


double calcrotmat(double nx, double ny, double nz, double *rmat)

	/* Find the rotation matrix that rotates |n| radians about
		the vector given by nx,ny,nz				*/
{
double ar, ai, br, bi, hp, cp, sp;
double arar, aiai, arai2, brbr, bibi, brbi2, arbi2, aibr2, arbr2, aibi2;
double phi;

phi = sqrt(nx*nx+ny*ny+nz*nz);

if (phi == 0.0)
	{
	*rmat++ = 1;
	*rmat++	= 0;
	*rmat++ = 0;
	*rmat++ = 0;
	*rmat++ = 1;
	*rmat++	= 0;
	*rmat++ = 0;
	*rmat++ = 0;
	*rmat++ = 1;
	}

/*printf("calcrotmat(%6.3f,%6.3f,%6.3f) -> phi = %6.3f\n",nx,ny,nz,phi);*/

else
	{
	/* First define Cayley-Klein parameters 	*/
	hp = phi/2;		
	cp = cos(hp);
	sp = sin(hp)/phi;	/* /phi because n is unit length in defs. */
	ar = cp;
	ai = -nz*sp;
	br = ny*sp;
	bi = -nx*sp;

 	/* Make auxiliary variables to speed this up	*/

	arar = ar*ar;
	aiai = ai*ai;
	arai2 = 2*ar*ai;
	brbr = br*br;
	bibi = bi*bi;
	brbi2 = 2*br*bi;
	arbi2 = 2*ar*bi;
	aibr2 = 2*ai*br;
	arbr2 = 2*ar*br;
	aibi2 = 2*ai*bi;


	/* Make rotation matrix.	*/

	*rmat++ = arar-aiai-brbr+bibi;
	*rmat++ = -arai2-brbi2;
	*rmat++ = -arbr2+aibi2;
	*rmat++ =  arai2-brbi2; 
	*rmat++ = arar-aiai+brbr-bibi;
	*rmat++ = -aibr2-arbi2;
	*rmat++ =  arbr2+aibi2;
	*rmat++ =  arbi2-aibr2;
	*rmat++ = arar+aiai-brbr-bibi;
	}
}



void zerovec(double *vec)

/*	Set a 3x1 vector to all zeros	*/

{
*vec++=0;
*vec++=0;
*vec++=0;
}



int blochsim(double *b1real, double *b1imag, 
		double *xgrad, double *ygrad, double *zgrad, double *tsteps, 
		int ntime, double *e1, double *e2, double df, 
		double dx, double dy, double dz, 
		double *mx, double *my, double *mz, int mode)

	/* Go through time for one df and one dx,dy,dz.		*/

{
int count;
int tcount;
double gammadx;
double gammady;
double gammadz;
double rotmat[9];
double amat[9], bvec[3];	/* A and B propagation matrix and vector */
double arot[9], brot[3];	/* A and B after rotation step. */
double decmat[9];		/* Decay matrix for each time step. */
double decvec[3];		/* Recovery vector for each time step. */
double rotx,roty,rotz;		/* Rotation axis coordinates. */
double mstart[3];
double mfinish[3];
double imat[9], mvec[3];

eyemat(amat); 		/* A is the identity matrix.	*/
eyemat(imat); 		/* I is the identity matrix.	*/

zerovec(bvec);
zerovec(decvec);
zeromat(decmat);

gammadx = dx*GAMMA;	/* Convert to Hz/cm */
gammady = dy*GAMMA;	/* Convert to Hz/cm */
gammadz = dz*GAMMA;	/* Convert to Hz/cm */

for (tcount = 0; tcount < ntime; tcount++)
	{
		/*	Rotation 	*/

	rotz = (*xgrad++ * gammadx + *ygrad++ * gammady + *zgrad++ * gammadz +
								df*TWOPI ) * *tsteps;
	rotx = (- *b1imag++ * GAMMA * *tsteps);
	roty = (*b1real++ * GAMMA * *tsteps++);
	calcrotmat(rotx, roty, rotz, rotmat);


	multmats(rotmat,amat,arot);
	multmatvec(rotmat,bvec,brot);


		/* 	Decay	*/

	decvec[2]= 1- *e1;
	decmat[0]= *e2;
	decmat[4]= *e2++;
	decmat[8]= *e1++;
	
	multmats(decmat,arot,amat);
	multmatvec(decmat,brot,bvec);
	addvecs(bvec,decvec,bvec);

	/*
	printf("rotmat = [%6.3f  %6.3f  %6.3f ] \n",rotmat[0],rotmat[3],
	  			rotmat[6]);
	printf("         [%6.3f  %6.3f  %6.3f ] \n",rotmat[1],rotmat[4],
				rotmat[7]);
	printf("         [%6.3f  %6.3f  %6.3f ] \n",rotmat[2],rotmat[5],
				rotmat[8]);
	printf("A = [%6.3f  %6.3f  %6.3f ] \n",amat[0],amat[3],amat[6]);
	printf("    [%6.3f  %6.3f  %6.3f ] \n",amat[1],amat[4],amat[7]);
	printf("    [%6.3f  %6.3f  %6.3f ] \n",amat[2],amat[5],amat[8]);
	printf(" B = <%6.3f,%6.3f,%6.3f> \n",bvec[0],bvec[1],bvec[2]);
	printf("<mx,my,mz> = <%6.3f,%6.3f,%6.3f> \n",
		amat[6] + bvec[0], amat[7] + bvec[1], amat[8] + bvec[2]);

	printf("\n");
	*/
	}

	/* Assume we started at zero magnetization. */

if (mode==0)		/* Indicates start at given m, or m0. */
	{
	mstart[0] = *mx;
	mstart[1] = *my;
	mstart[2] = *mz;

	multmatvec(amat,mstart,mfinish);
	addvecs(mfinish,bvec,mfinish);

	*mx = mfinish[0];
	*my = mfinish[1];
	*mz = mfinish[2];
	}
else if (mode==1)	/* Indicates to find steady-state magnetization */
	{
	scalemat(amat,-1.0);		/* Negate A matrix 	*/
	addmats(amat,imat,amat);	/* Now amat = (I-A)		*/
	invmat(amat,imat);		/* Reuse imat as inv(I-A) 	*/
	multmatvec(imat,bvec,mvec);	/* Now M = inv(I-A)*B		*/
	*mx = mvec[0];
	*my = mvec[1];
	*mz = mvec[2];
	}


}



int blochsimfz(double *b1real, double *b1imag, double *xgrad, double *ygrad, double *zgrad, 
		double *tsteps, 
		int ntime, double t1, double t2, double *dfreq, int nfreq,
		double *dxpos, double *dypos, double *dzpos, int npos, 
		double *mx, double *my, double *mz, int mode)


{
int count;
int poscount;
int fcount;
int totpoints;
int totcount = 0;

double *e1;
double *e2;
double *e1ptr;
double *e2ptr;
double *tstepsptr;
double *dxptr, *dyptr, *dzptr;

	/* First calculate the E1 and E2 values at each time step. */

e1 = (double *) malloc(ntime * sizeof(double));
e2 = (double *) malloc(ntime * sizeof(double));

e1ptr = e1;
e2ptr = e2;
tstepsptr = tsteps;

for (count=0; count < ntime; count++)
	{
	*e1ptr++ = exp(- *tstepsptr / t1);
	*e2ptr++ = exp(- *tstepsptr++ / t2);
	}

totpoints = npos*nfreq;

for (fcount=0; fcount < nfreq; fcount++)
    {
    dxptr = dxpos;
    dyptr = dypos;
    dzptr = dzpos;
    for (poscount=0; poscount < npos; poscount++)

	{
	blochsim(b1real, b1imag, xgrad, ygrad, zgrad, tsteps, ntime, 
			e1, e2, *dfreq, *dxptr++, *dyptr++, 
			*dzptr++, mx++, my++, mz++, mode);
		
	totcount++;
	if ((totpoints > 40000) && ( ((10*totcount)/totpoints)> (10*(totcount-1)/totpoints) ))
		printf("%d%% Complete.\n",(100*totcount/totpoints));
	}
    dfreq++;
    }

free(e1);
free(e2);

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

	/* bloch(b1,gradxyz,dt,t1,t2,df,dx,dy,dz,mode,mx,my,mz) */
{
double *b1r;	/* Real-part of B1 field.	*/
double *b1i;	/* Imag-part of B1 field.	*/
double *gx;	/* X-axis gradient. 		*/
double *gy;	/* Y-axis gradient. 		*/
double *gz;	/* Z-axis gradient. 		*/
double *tp;	/* Time steps (s)		*/
double t1;	/* T1 time constant (s)	 	*/
double t2;	/* T2 time constant (s)		*/
double *df;	/* Off-resonance Frequencies (Hz)	*/
double *dx;	/* X Positions (cm)			*/
double *dy;	/* Y Positions (cm)			*/
double *dz;	/* Z Positions (cm)			*/
double md;	/* Mode - 0=from M0, 1=steady-state	*/
double tstep;	/* Time step, if single parameter */
double *mxin;	/* Input points */
double *myin;
double *mzin;

double *mxout;	/* Input points */
double *myout;
double *mzout;

double *mx;	/* Output Arrays */
double *my;
double *mz;

int gyaflag=0;	/* 1 if gy was allocated. */ 
int gzaflag=0;	/* 1 if gy was allocated. */ 
int dyaflag=0;	/* 1 if dy was allocated. */ 
int dzaflag=0;	/* 1 if dy was allocated. */ 

int ntime;	/* Number of time points. 	 */
int ngrad;	/* Number of gradient dimensions */
int nf;
int npos;	/* Number of positions.  Calculated from nposN and nposM, depends on them. */
int nposM;	/* Height of passed position matrix. */
int nposN;	/* Width of passed position matrix. */
int nfnpos;	/* Number of frequencies * number of positions. */
int count;



#ifdef DEBUG
  printf("---------------------------------------\n");
  printf("3D-position + frequency Bloch Simulator\n");
  printf("---------------------------------------\n\n");
#endif

ntime = mxGetM(prhs[0]) * mxGetN(prhs[0]);	/* Number of Time, RF, and Grad points */


/* ====================== RF (B1) =========================
 * :  If complex, split up.  If real, allocate an imaginary part. ==== */
if (mxIsComplex(prhs[0]))
	{
	b1r = mxGetPr(prhs[0]);
	b1i = mxGetPi(prhs[0]);
	}
else
	{
	b1r = mxGetPr(prhs[0]);
	b1i = (double *)malloc(ntime * sizeof(double));
	for (count=0; count < ntime; count++)
		b1i[count]=0.0;
	}
#ifdef DEBUG
  printf("%d B1 points.\n",ntime);
#endif


/* ======================= Gradients ========================= */

ngrad = mxGetM(prhs[1]) * mxGetN(prhs[1]);	/* Number of Time, RF, and Grad points */
gx = mxGetPr(prhs[1]);				/* X-gradient is first N points. */

if (ngrad < 2*ntime)		/* Need to allocate Y-gradient. */
	{
	#ifdef DEBUG
	  printf("Assuming 1-Dimensional Gradient\n");
	#endif
	gy = (double *)malloc(ntime * sizeof(double));
	gyaflag=1;
	for (count=0; count < ntime; count++)
		gy[count]=0.0;
	}
else
	{
	#ifdef DEBUG
	  printf("Assuming (at least) 2-Dimensional Gradient\n");
	#endif
	gy = gx + ntime;	/* Assign from Nx3 input array. */
	}

if (ngrad < 3*ntime)		/* Need to allocate Z-gradient. */
	{
	gz = (double *)malloc(ntime * sizeof(double));
	gzaflag=1;
	for (count=0; count < ntime; count++)
		gz[count]=0.0;
	}
else
	{
	#ifdef DEBUG
	  printf("Assuming 3-Dimensional Gradient\n");
	#endif
	gz = gx + 2*ntime; 	/* Assign from Nx3 input array. */
	}

	/* Warning if Gradient length is not 1x, 2x, or 3x RF length. */

	
#ifdef DEBUG
  printf("%d Gradient Points (total) \n");
#endif
if ( (ngrad != ntime) && (ngrad != 2*ntime) && (ngrad != 3*ntime) )
		printf("Gradient length differs from B1 length\n");


if (gx == NULL) 
	printf("ERROR:  gx is not allocated. \n");
if (gy == NULL) 
	printf("ERROR:  gy is not allocated. \n");
if (gz == NULL) 
	printf("ERROR:  gz is not allocated. \n");



/* === Time points ===== */

tp = mxGetPr(prhs[2]);
if (mxGetM(prhs[2]) * mxGetN(prhs[2]) == 1)
	{
	tp = (double *)malloc(ntime * sizeof(double));
	tstep = *(mxGetPr(prhs[2]));
	for (count =0; count < ntime; count++)
		tp[count]=tstep;
	}
else if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != ntime)
	printf("Time-point length differs from B1 length\n");



/* === Relaxation Times ===== */

t1 = *mxGetPr(prhs[3]);
t2 = *mxGetPr(prhs[4]);

/* === Frequency Points ===== */

df = mxGetPr(prhs[5]);
nf = mxGetM(prhs[5]) * mxGetN(prhs[5]);
	
#ifdef DEBUG
  printf("%d Frequency points.\n",nf);
#endif


/* === Position Points ===== */

nposM = mxGetM(prhs[6]);
nposN = mxGetN(prhs[6]);

#ifdef DEBUG
  printf("Position vector is %d x %d. \n",nposM,nposN);
#endif

if (nposN==3)			/* Assume 3 position dimensions given */
	{
	npos = nposM;
	#ifdef DEBUG
	  printf("Assuming %d 3-Dimensional Positions\n",npos);
	#endif
	dx = mxGetPr(prhs[6]);
	dy = dx + npos;
	dz = dy + npos;
	}

else if (nposN==2)		/* Assume only 2 position dimensions given */
	{
	npos = nposM;
	#ifdef DEBUG
	  printf("Assuming %d 2-Dimensional Positions\n",npos);
	#endif
	dx = mxGetPr(prhs[6]);
	dy = dx + npos;
	dz = (double *)malloc(npos * sizeof(double));
	dzaflag=1;
	for (count=0; count < npos; count++)
		dz[count]=0.0;
	}

else				/* Either 1xN, Nx1 or something random.  In all these
				   cases we assume that 1 position is given, because it
				   is too much work to try to figure out anything else! */
	{
	npos = nposM * nposN;
	#ifdef DEBUG
	  printf("Assuming %d 1-Dimensional Positions\n",npos);
	#endif
	dx = mxGetPr(prhs[6]);
	dy = (double *)malloc(npos * sizeof(double));
	dz = (double *)malloc(npos * sizeof(double));
	dyaflag=1;
	dzaflag=1;
	for (count=0; count < npos; count++)
		{
		dy[count]=0.0;
		dz[count]=0.0;
		}
	#ifdef DEBUG
	  if ((nposM !=1) && (nposN!=1))		
		{
		printf("Position vector should be 1xN, Nx1, Nx2 or Nx3. \n");
		printf(" -> Assuming 1 position dimension is given. \n");
		}	
	#endif
	}

if (dx == NULL) 
	printf("ERROR:  dx is not allocated. \n");
if (dy == NULL) 
	printf("ERROR:  dy is not allocated. \n");
if (dz == NULL) 
	printf("ERROR:  dz is not allocated. \n");

nfnpos = nf*npos;	/* Just used to speed things up below. 	*/ 


/* ===== Mode, defaults to 0 (simulate from starting postion). ====== */

if (nrhs > 7)
	md = *mxGetPr(prhs[7]);		/* (Mode of 1 simulates steady state. */
else
	md = 0;

#ifdef DEBUG
if (md==0)
	printf("Simulation from Initial Condition.\n");
else
	printf("Simulation of Steady-State.\n");
#endif


/* ===== Allocate Output Magnetization vectors arrays.	*/

plhs[0] = mxCreateDoubleMatrix(npos,nf,mxREAL);		/* Mx, output. */
plhs[1] = mxCreateDoubleMatrix(npos,nf,mxREAL);		/* My, output. */
plhs[2] = mxCreateDoubleMatrix(npos,nf,mxREAL);		/* Mz, output. */

mx = mxGetPr(plhs[0]);
my = mxGetPr(plhs[1]);
mz = mxGetPr(plhs[2]);

mxout = mx;
myout = my;
mzout = mz;

/* ===== If Initial Magnetization is given... */

if ( (nrhs > 10) &&	 
	(mxGetM(prhs[8]) * mxGetN(prhs[8]) == nfnpos) &&
        (mxGetM(prhs[9]) * mxGetN(prhs[9]) == nfnpos) &&
        (mxGetM(prhs[10]) * mxGetN(prhs[10]) == nfnpos)  )

		/* Set output magnetization to that passed. */

		
	{
	#ifdef DEBUG
  	  printf("Using Speicified Initial Magnetization.\n");
	#endif

	mxin = mxGetPr(prhs[8]);
	myin = mxGetPr(prhs[9]);
	mzin = mxGetPr(prhs[10]);
	for (count =0; count < nfnpos; count++)
		{
		*mxout++ = *mxin++;
		*myout++ = *myin++;
		*mzout++ = *mzin++;
		}
	}
else 
	{
	#ifdef DEBUG
	if (nrhs > 10) 	 /* Magnetization given, but wrong size! */
		{
		printf("Initial magnetization passed, but not Npositions x Nfreq. \n");
		}
	  printf(" --> Using [0; 0; 1] for initial magnetization. \n");
	#endif
	for (count =0; count < nfnpos; count++)
		{
		*mxout++ = 0;	/* Set magnetization to Equilibrium */
		*myout++ = 0;
		*mzout++ = 1;
		}
	}	


/* ======= Do The Simulation! ====== */

#ifdef DEBUG
  printf("Calling blochsimfz() function in Mex file.\n");
#endif

blochsimfz(b1r,b1i,gx,gy,gz,tp,ntime,t1,t2,df,nf,dx,dy,dz,npos,mx,my,mz,(int)(md));



/* ====== Free up allocated memory, if necessary. ===== */

if (!mxIsComplex(prhs[0]))
	free(b1i);	/* We had to allocate this before. */

if (mxGetM(prhs[2]) * mxGetN(prhs[2]) == 1)
	free(tp);	/* We had to allocate this. */

if (dyaflag==1)
	free(dy);
if (dzaflag==1)
	free(dz);
if (gyaflag==1)
	free(gy);
if (gzaflag==1)
	free(gz);

}





