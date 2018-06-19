#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#define _LINUX
//#include "libvec.h"
#include "giv.h"
#define ETA DBL_EPSILON
#define THETA DBL_MAX

#define MAX(x, y) ( ( (x) >= (y) )? (x) : (y) )
#define SQRE(x) (x) * (x)
#define SIGN(x, y) ( ((y) >= 0) ? (x) : -(x) )
#define MIN(x, y) ( ((x) <= (y)) ? (x) : (y) )
#define MAXITER 30
#define ONE 1

int idif;
double fact;


int number_below(double cut_off, double *root, int n)
				/* counts eigenvalues below cut_off. */
				/* introduced on August 12, 1993 */
{
  while(n && root[n] > cut_off)
    n--;
  return n;
}

int rtcomp(double **proot1, double **proot2)
     /**********************************************************************
      *
      *	function rtcomp
      *		compares two elements
      *
      *	input:	proot1	pointer to pointer to element 1
      *		proot2	pointer to pointer to element 2
      *
      *	returns: -1, 0, or 1 if element 1 is less than, equal to, or
      *		 larger than element 2, respectively.
      *************************************************************************/
     
{
  if(**proot1 < **proot2)
    return -1;
  else if(**proot1 > **proot2)
    return 1;
  else
    return 0;
}

void vecopy(double dest[], double sourc[], int n)
     /********************************************************************
      *
      *	function vecopy
      *		copies a linear vector into another
      *
      *	input:	sourc	source vector
      *		n	dimension of vectors
      *	
      *	output:	dest	destination vector
      ***********************************************************************/
{
  register int i;
  
  for(i = 1; i <= n; i++)
    dest[i] = sourc[i];
}

#define GSIZE(x)  ( (unsigned)((x) * sizeof(double)) )
#define ISIZE(x)  ( (unsigned)((x) * sizeof(int)) )

int allocg(double **a, int n, int dim)
{
  while(n >= 0){
    if( (a[n] = (double*)calloc(dim, sizeof(double)) ) == NULL){
      printf("\nInsufficient memory to allocm\n");
      return -1;
    }
    --n;
  }
  return 0;
}

void freeg(double **a, int n)
{	
  while(n >= 0){
    free(a[n]);
    n--;
  }
}

int alloci(int **a, int n, int dim)
{	
  register int i;
  
  while(n >= 0){
    if( (a[n] = (int*)malloc(ISIZE(dim)) ) == NULL){
      printf("\nInsufficient memory for allocm\n");
      return -1;
    }
    for(i = 0; i < dim; i++)
      a[n][i] = 0;
    --n;
  }
  return 0;
}

void freei(int **a, int n)
{	
  while(n >= 0)
    free(a[n--]);
}

#define quit {printf("\nInsufficient memory (sort)\n"); return -1;}
#ifndef __CDECL
#define __CDECL 
#endif
void sort(int n, int nroot, double *root, double *vect[])
     /*******************************************************************
      *
      *	function sort
      *		sorts eigenvalues into increasing order.
      *		orders corresponding eigenvectors accordingly.
      *
      *	input: 	n 	dimension of eigenvectors
      *		nroot 	number of eigenvalues to be sorted
      *		root	vector containing eigenvalues
      *		vect	matrix (nroot X n) containing eigenvectors
      *
      *	output:	root	ordered eigenvalues
      *		vect	ordered eigenvectors
      *
      *	calls:	qsort
      *		malloc
      *		free
      *		printf
      *
      *		rtcomp
      *
      *		written by	L.N.Oliveira
      ***********************************************************************/
{
  register int i;
  
				/* following loop copies eigenvalues */
				/* into matrix vect	*/
  for(i = 1; i <= nroot; i++)
    *vect[i] = root[i];
  
				/* next sort lines of vect in increasing */
				/* order of their first elements */
				/* (the eigenvalues)	*/
  qsort((void*) (vect + 1), (unsigned)nroot, sizeof(double *),
			 (int (__CDECL *) (const void*, const void*)) rtcomp);
   
				/* copy sorted eigenvalues back into root */
  for(i = 1; i <= nroot; i++)
    root[i] = *vect[i];
  
  
  return;
}

int norm_a(int n, int nsize, int nroot, double *anormp,
		 double *tracep, double *a, double *root, double *vect[])
{
  register int i, j;
  int nm1, nm2, np1, jump = 1;
  double floatn, factor = 0.0;
  double scale;
  
  floatn = (double)n;
  nm1=n-1;
  nm2=n-2;
  np1=n+1;
  
  for(j=2; j<=np1; j++){
    *tracep += a[jump];
    jump += j;
  }
  jump=1;
  *tracep /= floatn;
  for(j = 2; j <= np1; j++){
    a[jump] -= *tracep;
    jump += j;
  }
  
  for(i = 1; i <= nsize; i++){
    scale = fabs(a[i]);
    factor = MAX(factor, scale);
  }
  
  if( factor == 0. ){
    for(i=1; i <= abs(nroot); i++){
      if(nroot >= 0){
	for(j = 1; j <= n; j++)
	  vect[i][j] = 0.0;
	vect[i][i] = 1.0;
      }
      root[i] = 0;
    }
    return -1;
  }
  
  
  scale = 1./factor;
  for(i=1; i<=nsize; i++)
    *anormp += SQRE(a[i]) * SQRE(scale);
  *anormp += *anormp;
  jump = 1;
  for(j = 2; j <= np1; j++ ){
    *anormp -= SQRE(a[jump]) * SQRE(scale);
    jump += j;
  }
  *anormp= factor * sqrt(*anormp);
  scale = 1./ *anormp;
  for( i=1; i<=nsize; i++ )
    a[i] *= scale;
  return 0;
}

void rstroot(double anorm, int n, double trace, double *root)
				// restores roots by scaling them by
				// anorm and shifting them by trace
{	
  do{
    root[n] *= anorm;
    root[n--] += trace;
  }
  while(n);
}

double dsum(double *b, double *a, int ip1, int limit)
{
  register int ii = ip1 -1;
  double dsumc=0.0;
  
  while(ii < limit){
    dsumc += b[ii] * *a;
    a += ++ii;
  }
  return dsumc;
}

void rotate(double *v[], double c, double s, int njx, int jtop)
{
  double temp, *pv, *pv1, *plimit;
  
  pv = v[0] + 1;
  pv1 = v[1] + 1;
  plimit = pv + jtop;
  
  while(pv <= plimit){
    temp = *pv;
    *(pv++) = temp * c + *pv1 * s;
    *pv1 = *pv1 * c - temp * s;
    pv1++;			/* pv1 must be incremented after rotation, */
				/* for either side of previous assignment  */
				/* may be evaluated first; bug caught on */
				/* 2 April '93 */
  }
}

double dot(double *a, double *b)
{
  double dotc;
#ifdef _INCLUDE_HPUX_SOURCE
  int itop = idif + 1;
  dotc = vec_$ddot(a, b, &itop);
#else
  double *itop;
  
  itop = a + idif ;
  dotc = 0.0;
  while( a <= itop)
    dotc += *(a++) * *(b++);
#endif
  return dotc;
}

void vecsum(double *a, double *b)
{
#ifdef _INCLUDE_HPUX_SOURCE
  int itop = idif + 1;

  vec_$dmult_add(a, b, &itop, &fact, a);
#else
  double *itop;
  itop = a + idif;
  while( a <= itop)
    *(a++) += *(b++) * fact;
#endif
}

double inivec(int i, int n, double *b[], double *vect[],
		 double delta, double delbig, double temp)
				// returns temp, used in bkiter
				// bug caught 13 Jan 94
{	
  register int k, l;
  double ret = temp;
  double elim1;
  
  for(l = n; l >= 1; l--){
    elim1=vect[i][l];
    if( l!=n ){
      elim1-=vect[i][l+1]*b[3][l];
      if( l!=n-1 )
	elim1-=vect[i][l+2]*b[4][l];
    }
    
    if( fabs(elim1) > delbig ){
      for(k=1; k<=n; k++)
	vect[i][k]/=delbig;
      l++;
      continue;			// do l loop over
    }
    
    if( fabs(ret = b[2][l]) < delta )
      ret = delta;
    vect[i][l] = elim1 / ret;
  }
  return ret;
}

void bkiter(int n, int nroot, double *a, double *b[], double *root,
	    double *vect[], double delta, double delbig,
	    double del1, double theta1)
{	
  register int i, j;
  int nm1 = n - 1;
  int jump;
  double aroot, elim1, elim2;
  double temp = 0.0;

  for(i = 1; i<=nroot; i++){
    for(j = 1; j<=n; j++)
      vect[i][j] = 1.0;
  }
  
  for(i = 1; i<=nroot; i++){
    aroot = root[i];
    elim1 = a[1]-aroot;
    elim2 = b[1][1];
    jump = 1;
    for(j = 1; j<=nm1; j++){
      jump += j+1;
      if( fabs(elim1) > fabs(b[1][j]) ){
	b[2][j] = elim1;
	b[3][j] = elim2;
	b[4][j] = 0.0;
	temp = b[1][j]/elim1;
	elim1 = a[jump]-aroot-temp*elim2;
	elim2 = b[1][j+1];
      }
      
      else{
	b[2][j] = b[1][j];
	b[3][j] = a[jump] - aroot;
	b[4][j] = b[1][j+1];
	temp = 1.0;
	if(fabs(b[1][j]) > theta1)
	  temp = elim1 / b[1][j];
	elim1 = elim2 - temp * b[3][j];
	elim2 = -temp * b[1][j+1];
      }
      
      b[0][j]=temp;
    }
    
    b[2][n] = elim1;
    b[3][n] = b[4][n] = b[4][nm1] = 0.0;
    
    temp = inivec(i, n, b, vect, delta, delbig, temp);
    
    elim1 = vect[i][1];
    for(j = 1; j<=nm1; j++){
      if(b[2][j] != b[1][j]){
	vect[i][j] = elim1;
	elim1 *= -b[0][j];
	elim1 += vect[i][j+1];
      }
      else{
	vect[i][j] = vect[i][j+1];
	elim1 -= vect[i][j+1]*temp;
      }
    }
    vect[i][n] = elim1;
    
    inivec(i, n, b, vect, delta, delbig, temp);
    
    elim1 = 0.0;
    for(j = 1; j<=n; j++)
      elim1 = MAX(fabs(vect[i][j]),elim1);
    temp = 0.0;
    for(j = 1; j<=n; j++){
      elim2 = vect[i][j]/elim1;
      temp+=SQRE(elim2);
    }
    
    temp = 1./sqrt(temp)/elim1;
    for( j = 1; j <= n; j++){
      vect[i][j] *= temp;
      if(fabs(vect[i][j]) < del1)
	vect[i][j] = 0.0;
    }
  }
}

int qrit(int *nx, int *ni, int *k, double *a, double *b, double *sh,
		double tol, int m, int ntop)
{
  int i;
  double at, cx, cxs, disc;
  double g, ps, rs, st, sx, u;
  double b_i;
  
  
  (*k)++;
  if( *k >= m ){
    printf("NO CONVERGENCE OF QR ALGORITHM IN %d ITERATIONS", *k);
    return -1;
  }

				/*
    Calculate partial shift st and add to total shift *sh
    st is the smaller eigenvalue of the 2 X 2 matrix whose
    diagonal elements are a[*nx] and a[*nx -1]
			       */
  at= a[*nx]+a[(*nx) - 1];
  st=at/2.0;
  disc=SQRE(at)-4.*( a[*nx]*a[(*nx) - 1] - b[(*nx) - 1] );
  if( disc>0. )
    st -= ((st >= 0) ? sqrt(disc) : -sqrt(disc))/2.0;
  *sh += st;
  for(i=1; i<=*nx; i++)
    a[i]-=st;
				/*
    Now proceed to QR transform
    Algorithm described in Wilkinson (Algebraic Eigenvalue Problem)
    p. 567
			        */
  
  u = 0.0;
  cxs = 1.0;
  cx = 1.0;
  sx = 0.0;
  g = 0.0;
  for(i = *ni; i <= *nx; i++){
    g = a[i] * cx - g * sx;
    if(cx <= tol)
      ps = b[i - 1] * cxs;	/* never reached when i= *ni */
    else if(fabs(g) > tol)
      ps = SQRE(g) / cx;
    else{
      ps = 0.0;
      if(i == *nx){		/* diversion to avoid */
	b[i - 1] = 0.0;		/* division by rs = 0 */
	sx = 0.0;
	cxs = cx;
	cx = 1.0;
	u = 0.0;
	a[i] = g;
	continue;
      }
    }
    
    b_i = (i == *nx) ? 0.0 : b[i];
    rs = ps + b_i;
    if(i > *ni)
      b[i - 1] = sx * rs;
    sx = b_i / rs;
    cxs = cx;
    cx = ps / rs;
    u = (i == *nx) ? 0.0 : sx * (g + a[i+1]); // conditional introduced
					      // on 13 Jan 94
					      // to avoid reference to
					      // a[*nx + 1], which is
					      // undefined
    a[i] = g + u;
  }
  
  return 0;
}

int swtch1(int *nx, int *k, double *a, double *sh)
{
  
  a[(*nx)--]+=*sh;
  *k = 1;
  return 0;
}

int swtch2(int *nx, int *k, double *a, double *b, double *sh,
	   double small)
{
  double
    al,
    am,
    ams,
    an,
    cx,
    sx,
    sam,
    ta,
    tb,
    tc;
  
  al=b[(*nx) - 1];
  am=(a[(*nx) - 1]-a[*nx])/2.0;
  if(fabs(al) > small * fabs(am)){
    ams=SQRE(am);
    sam= (am>=0)? 1. : -1;
    an=sqrt( al+SQRE(am) );
    cx=( an+fabs(am) )/2./an;
    sx=b[(*nx) - 1]/(4.*SQRE(an)*cx);
    ta=a[(*nx) - 1];
    tb=a[*nx];
    tc=b[(*nx) - 1];
  
    a[(*nx) - 1]=ta*cx+tb*sx+tc*sam/an+*sh;
    a[*nx]=ta*sx+tb*cx-tc*sam/an+*sh;
    b[(*nx) - 1]=4.*ams*cx*sx-fabs(am)*tc/an+tc*SQRE(cx-sx);
  }
  else{				/* following detour avoids division by zero */
				/* when b[(*nx) -1] is very small */
				/* introduced on 4 June '93 */
    a[(*nx) - 1] += *sh;
    a[*nx] += *sh;
    b[(*nx) -1] = 0.0;
  }
  
  *k=1;
  *nx-=2;
  
  return 0;
}

int swtchg(int *nx, int *ni, int *k, double *a, double *b,
				double *sh, double tol, int m, int ntop, int *it)
{
  do{				/* this loop repeated until |b[*it]|<=tol */
				/* or error condition occurs in function qrit*/
    for(*it = (*nx) - 1; *it > *ni; --(*it)){
      if( fabs(b[*it]) <= tol )
	return 0;
    }
  }
  while( qrit(nx, ni, k, a, b, sh, tol, m, ntop)+ONE );
				/* here if error condition in qrit */
  return -1;
}

int evqr(double *a, double *b, int n, int m, double small)
{
  int nx, ni, k, ntop = 0, it = 0;
  int size, itloop;
  double sh;
  nx=n;
  ni=1;
  sh=0.0;
  k=0;
  for(itloop=0; ; ){
    if(!itloop){
      switch(nx-ni){
      case 0:
	size=1;
	break;
	
      case 1:
	size=2;
	break;
	
      default:
	size=3;
	break;
      }
    }
    else{
      switch(nx-it-1){
      case 0:
	size=1;
	break;
	
      case 1:
	size=2;
	break;
	
      default:
	ni=it+1;
	size=3;
	break;
      }
    }
    
    if( size==1 )
      swtch1(&nx, &k, a, &sh);
    else if( size==2 )
      swtch2(&nx, &k, a, b, &sh, small);
    if(  size<=2 ){
      if( nx<ni ){
	if( ni==1 )
	  return 0;
	ni=1;
      }
      itloop=0;
    }
    else{
      if( swtchg(&nx, &ni, &k, a, b, &sh, small, m, ntop, &it) )
	return -1;
      itloop=1;
    }
  }
}

void size_2(double *a, double *b, double *v[], int *pn, double *psh,
	    double *pst, int nnm1, int njx, double small)
{
  double
    al,
    am,
    an,
    c,
    s,
    ta,
    tb,
    tc,
    cx,
    sx,
    cs;
  
  b += *pn - 1;
  a += *pn;
  
  al = -*b;
  am = (*(a-1) - *a)/2.0;
  if(fabs(al) > small * fabs(am)){
    an = hypot(al, am);
    c = sqrt((an + fabs(am))/2./an);
    s = ((am >= 0) ? 0.5 : -0.5)*al/an/c;
    ta = *(a-1);
    tb = *a;
    tc = *b;
    cx = c * c;
    sx = s * s;
    cs = c * s;
  
    *(a-1) = ta *cx + tb * sx - 2.* tc * cs;
    *a = ta * sx + tb * cx + 2. * tc * cs + *psh;
    *b = 2. * am * cs + tc * (cx - sx);
    s = -s;
  }
  else{				/* following detour avoids division by zero */
				/* when *b is very small */
				/* introduced on 4 June '93 */
    c = 1.0;
    s = 0.0;
    *a += *psh;
    *b = 0.0;
  }
  
  rotate(&v[--*pn], c, s, njx, nnm1);
  return;
}




void size_g(int k, int n, int ni, double *a, double *b, double *psh,
	    double *pst, double *v[], int nnm1, int njx)
{
  int i, ntop;
  double at, disc, r, s, cs, c, u, g, q;
  double *limit;
  
  a += n;
  b += n - 1;
  
  if(k != 1){
    at = *a + *(a - 1);
    *pst = at/2.0;
    disc = at * at -4. * (*a * *(a-1) - SQRE(*b));
    if(disc > 0.)
      *pst -= ( (*pst >=0) ? sqrt(disc) : -sqrt(disc) )/2.0;
  }
  *psh += *pst;
  limit = a - n;
  while( a > limit)
    *(a--) -= *pst;
  a += ni;
  b += ni - n +1;
  r = hypot(*a, *b);
  s = *b/r;
  cs = s;
  c = *a/r;
  u = SQRE(s) * (*a + *(a+1));
  *a += u;
  rotate(&v[ni], c, s, njx, nnm1);
  ntop = n -2;
  for(i = ni; i <= ntop; i++){
    g = *(++a) - u;
    q = c * *a - cs * *b;
    r = hypot(q, *(++b));
    *(b-1) = s * r;
    s = *b/r;
    cs = c * s;
    c = q / r;
    u = SQRE(s) * (g + *(a+1));
    *a = g + u;
    rotate(&v[i+1], c, s, njx, nnm1);
  }
  *b = s * (c * *(a += -ntop - 1 + n) - cs * *b);
  *a -= u;
}

int gage_b(int n, int ni, double *b, int *ptiny, double tol)
{	  	
  while(--n > ni){
    if(fabs(b[n]) <= tol){
      *ptiny = 1;
      return n;
    }
  }
  *ptiny = 0;
  return n;
}

int qrtn(double *a, double *b, double *v[], double *eig, int n, int m,
     double tol, double small, int njx)
     
{	
  double
    sh,
    st;
  
  int
    ni,
    nnm1,
    k,
    it,
    tiny_b;
  
  nnm1 = n - 1;
  ni = 1;
  sh = 0;
  it = tiny_b = 0;
  k = 0;
  
  for(;;){
    switch(n - it){
    case 2:
      size_2(a, b, v, &n, &sh, &st, nnm1, njx, small);
      
    case 1:
      a[n] += sh;
      n--;
      k=0;
      if(n < ni){
	if(ni == 1)
	  return 0;
	ni = 1;
      }
      it = ni - 1;
      tiny_b = 0;
      break;
      
    default:
      ni = it + 1;
      if(!tiny_b)
	st = eig[n] - sh;
      while(++k < m){
	if(!tiny_b){
	  it = gage_b(n, ni, b, &tiny_b, tol);
	  if(tiny_b)
	    break;
	}
	size_g(k, n, ni, a, b, &sh, &st, v, nnm1, njx);
	tiny_b = 0;
      }
      if(!tiny_b){
	printf("\nNo convergence in %d iterations", k);
	return -1;
      }
    }
  }
}

void simvec(int nsize, int n, int nroot, double *a, double *b[],
		 int *b6, double *vect[])
{	
  int jump, limit, j1;
  register int k, im;
  
  jump = nsize - n - 1;
  for(im = n - 1; im > 1; im--){
    limit = b6[im-1];
    j1=jump;
    for(k = im; k <= limit; k++){
      b[2][k]=a[j1];
      j1 += k;
    }
    
    idif=limit-im;
    for(k=1; k<=nroot; k++){
      fact = -2. * dot(b[2]+im, vect[k]+im);
      vecsum(vect[k]+im, b[2]+im);
    }
    
    jump-=im;
  }
}


void tridia(int n, int nsize, double *a, double *b[], int *b6,
			double del1)
{	
  int j, ia = 1, ic, id = 0, jp2, j1, limit, limles, limlo, ii, ip1;
  register int i, jj;
  double dtemp, temp, sum, dak;
  
  for( j=1; j <= n-2; j++){
    ia += j+2;
    id += j+1;
    jp2 = j+2;
    j1 = j + 1;
    
    limit = j1;
    ii = ia;
    for(i = jp2; i <= n; i++){
      b[0][i] = a[ii];
      if(fabs(a[ii] ) > del1)
	limit = i;
      ii += i;
    }
    dtemp = a[id];
    
    if(limit <= j1){
      b[1][j] = dtemp;
      a[id] = 0.0;
      continue;
    }
    
    idif = limit - jp2;
    sum = dot(b[0] + jp2, b[0] + jp2 );
    sum = sqrt(sum+SQRE(dtemp));
    b[1][j] = -SIGN( sum, dtemp );
    b[2][j+1] = sqrt( (1. + fabs(dtemp)/sum) / 2. );
    temp = SIGN( 0.5/b[2][j+1]/sum, dtemp);
    ii = ia;
    
    for(i = jp2; i <= limit; i++){
      b[2][i] = a[ii] * temp;
      ii+=i;
    }
    
    dak=0.0;
    ic = id + 1;
    limles = limit-1;
    for(i = j1; i <= limles; i++){
      idif = i - j1;
      dtemp = dot( b[2]+j1, a+ic );
      ic += i;
      ip1 = i + 1;
      jj = ic + idif;
      dtemp += dsum( b[2]+1, a+jj, ip1, limit ); 
      dak += dtemp * b[2][i];
      b[1][i] = dtemp;
    }
    
    idif = limit - j1;
    dtemp = dot( b[2] + j1, a + ic );
    dak += dtemp * b[2][limit];
    b[1][limit] = dtemp;
    idif = limit - j1;
    if( limit != n ){
      ic += limit;
      limlo = limit + 1;
      for(i = limlo; i <= n; i++){
	b[1][i] = dot(b[2] + j1, a + ic);
	b[2][i] = 0.0;
	ic += i;
      }
    }
    fact = -dak;
    vecsum( b[1] + j1, b[2] + j1 );
    jj=id;
    
    for(i = j1; i <= n; i++){
      a[jj] = b[2][i];
      if( i <= limit ){
	fact = -2. * b[2][i];
	idif = i - j1;
	vecsum(a + jj + 1, b[1] + j1);
      }
      fact = -2. * b[1][i];
      idif = MIN(i, limit) - j1;
      vecsum(a + jj + 1, b[2] + j1);
      jj += i;
    }
    b6[j] = limit;
  }
}

int givens(int n, int nrootx, int njx, double *a, double *root,
	    double *vect[])
				/* nrootx is the number of eigenvectors */
				/* wanted. */
				/* should be negative if no eigenvector is */
				/* wanted and zero if a cut_off is stored in */
				/* root[1], in which case givens returns the */
				/* number of eigenvalues below the cut_off */
				/* and computes only those many eigenvectors */
{
  int
    nroot,
    nsize,
    jump=1,
    nm1,
    nm2,
    i,
    j,
    ntop,
    *b6,
    ret = 0;
  
  char degenerate = 0;
  
  double
    trace = 0.,
    anorm = 0.,
    tolsc,
    *b[5],
    cut_off = THETA;		/* default cut_off energy (alternatively, */
				/* cut_off may be stored in root[1] and */
				/* signaled by argument nrootx == 0 */

  const double
    del1 = ETA / 100.,
    delta = ETA * ETA * 100.,
    small = ETA * ETA / 100.,
    delbig = THETA * ETA * ETA / 10.,
    theta1 = 1000. / THETA,
    toler = 100. * sqrt(ETA),
    emag= ETA;


  nsize = (n * (n + 1) ) / 2;
  nm1 = n - 1;
  nm2 = n - 2;
  nroot = abs(nrootx);
  
  if(n < 1)
    return -1;

  if(!nroot){			/* possibility nroot == 0 accepted as of */
    cut_off = root[1];		/* August 12, 1993*/
  }
  
  if(n == 1){
    root[1] = a[1];
    if(nrootx > 0){		
      vect[1][1] = 1.0;
      return 0;
    }
    else if(!nrootx && root[1] < cut_off){
      vect[1][1] = 1.0;
      return 1;
    }
    else{
      return 0;
    }
  }
  
  norm_a(n, nsize, nrootx, &anorm, &trace, a, root, vect);
  if( allocg(b, 4, n + 1) )
    return(-1);
  
  /*****************************TRIDIAGONALIZATION SECTION****************/
  
  if( alloci(&b6, 0, n + 1) ){
    printf("\nInsufficient memory to alloc b6\n");
    return -1 ;
  }
  
  if(nm2)
    tridia(n, nsize, a, b, b6, del1);

  b[1][n - 1] = a[nsize - 1];	/* this line and next were erroneously */
				/* placed inside tridia()*/
  a[nsize - 1] = 0.0;		/* bug caught 4 June '93  */
  
  jump = 0;
  for(j = 1; j<=nm1; j++){
    jump += j;
    root[j] = a[jump];
    b[3][j] = SQRE( b[1][j] );
  }
  root[n] = a[nsize];
  if( evqr(root, b[3], n, MAXITER, small) )
    return -1;
  for(j = 1; j<=n; j++)
    b[2][j] = root[j];
  
  
  sort(0, n, root, vect);	/* second argument changed from nroot to n */
				/* on 30 March '93 */
  cut_off -= trace;		// shift cut_off
  cut_off /= anorm;		// and scale it, to make it comparable to
				// the shifted and scaled energies in root.
  if(!nroot){			//  cut_off introduced August 12, 1993
	 ret = nroot = number_below(cut_off, root, n);
    nrootx = (nroot) ? nroot: -1;
  }

  if(nrootx < 0){
    rstroot(anorm, n, trace,root);
    freeg(b, 4);		// 2 (minor) leaks plugged Jul 9 2006
    freei(&b6, 0);
    return 0;
  }
  else if(nrootx > 1){
    ntop = nroot-1;
    tolsc =  toler / anorm;
    for(i = 1; i<=ntop; i++){
      if( fabs( root[i+1]-root[i] )<= tolsc ){
	degenerate = 1;
	break;
      }
      
    }
  }
  
  if(degenerate){
    jump = 0;
    for(j = 1; j<=nm1; j++){
      jump += j;
      root[j] = a[jump];
      b[3][j] = b[1][j];
    }
    
    root[n] = a[nsize];
    
    for(i = 1; i<=n; i++)
      for(j  =  1; j <= n; j++)
	vect[i][j] = (i == j? 1. : 0.); /* eigenvectors zeroed */
				/* bug caught */
				/* 2 April '93*/
    
    if( qrtn(root, b[3], vect, b[2], n, 25, emag, small, njx) )
      return -1;
    
    sort(n, n, root, vect);	/* second argument changed from */
				/* nroot to n  on 4 / 2 / 93 */ 
    
  }
  
  else{				/* !degenerate */
    bkiter(n, nroot, a, b, root, vect, delta,delbig,del1, theta1);
  }
  
				/* simvec section */
  
  if(nm2)
    simvec(nsize, n, nroot, a, b, b6, vect);
  
  freeg(b, 4);
  freei(&b6, 0);

  rstroot(anorm, n, trace, root);
  
  return ret;
}

#define END -7.7e77
#define TOLER 1.e-4

void cpy(int dim, double *c, double *a)
				// copies matrix a 
				// (lower triangular, dimension dim)	       
				// into c (previously allocated).
				// a and c are defined from 0,
				// but start at 1
{
  int size = dim * (dim + 1) / 2;
#ifdef _INCLUDE_HPUX_SOURCE
  a++;
  c++;
  vec_$dcopy(a, c, &size);
#elif defined _LINUX__
  a++;
  c++;
  int step = 1;
  dcopy_(&size, a, &step, c, &step);
#else
  while(size--)
    *(++c) = *(++a);
#endif // _INCLUDE_HPUX_SOURCE
}

int check(int n, double *a, double *c, double *root, double *vect[],
      int nroot, int chk)
     /****************************************************************
      *	function check						*
      *		called with chk=0, copies matrix a to matrix c.	*
      *		called with chk=1, checks diagonalization.	*
      *								*
      *		upon entry:	n contain matrix dimension	*
      *				a contains matrix to be 	*
      *				  diagonalized	(chk=0)		*
      *				  (stored in triangular form).	*
      *				c contains copy of a (chk=1).	*
      *				root contains eigenvalues.	*
      *				vect contains eigenvectors.	*
      *				nroot is number of roots for	*
      *				      which there are		*
      *				      eigenvectors in vect.	*
      *								*
      *		returns:	0 if check OK.			*
      *				-1 otherwise			*
      *								*
      *		uses:		lin2sq()			*
      * ***************************************************************/
     
{
  int i;			/* row index */
  int j;			/* column index */
  register int k;
  
  double sum;			/* accumulates sum of squares of differences */
				/* between c X vect and root X vect */ 
  double *limit;		/* limits growth of matrix pointers */
  double *pvecti;		/* pointers to vect matrix */
  double *pvectj;
  double norm;			/* accumulates norm of eigenvector */
  double dot_ij;		/* accumulates scalar product of two */
				/* eigenvectors */
  double c_times_vect;		/* accumulates product matrix X eigenvector */
  double c_max;			/* maximum element in row i of matrix c */
  double lin2sq();		/* finds element (i,j) in linear array */
  double c_ik;			/* element (i,k) of matrix c */
  
//  limit= a + n * (n+1)/2;

  if( !chk){			/* copy a to c */
    cpy(n, c, a);
//    while(a < limit)
//      *(++c) = *(++a);
  }
  
  else{
    for(i = 1; i <= nroot; i++){
      norm = 0.0;
      sum = 0.0;
      c_max = 0.0;
      
      for(j = 1; j <= n; j++){
	norm += SQRE(vect[i][j]);
	c_times_vect = 0.0;
	dot_ij=0.0;
	
	pvectj = vect[j];
	limit = pvectj + n;
	
	for(k = 1; k <= n; k++){
	  if(fabs(c_ik = c[ (i > k) ? i*(i-1)/2 + k : k*(k-1)/2 + i]) > c_max)
	    c_max = fabs(c_ik);
	  
	  c_times_vect += *(++pvectj) * c_ik;
	}
	
	if( i < j && j <= nroot){
	  pvecti = vect[i];
	  pvectj = vect[j];
#ifdef _INCLUDE_HPUX_SOURCE
	  pvecti++;
	  pvectj++;
	  dot_ij = vec_$ddot(pvecti, pvectj, &n);
#else
	  while(pvectj < limit)
	    dot_ij += *(++pvecti) * *(++pvectj);
#endif
	  if( SQRE(dot_ij) > TOLER ){
	    printf("EIGENVECTORS %d & %d NONORTHOGONAL:"
		   " %15.7g", i, j, SQRE(dot_ij));
	    return -1;
	  }
	  
	}
	if(c_max > DBL_MIN & j <= nroot)
	  sum += fabs(vect[j][i] * root[j] - c_times_vect)/c_max;
      }
      
      if(SQRE(norm - 1.) > TOLER ){
	printf(" EIGENVECTOR %d NOT NORMALIZED", i);
	return -1;
      }
      
      if( sum > TOLER){
	printf(" VECTOR %d NO EIGENVECTOR", i);
	return -1;
      }
    }
  }
  return 0;
}

