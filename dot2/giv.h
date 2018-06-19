#ifndef _GIV_
#define _GIV_
				/* n = dimensao da matriz */
				/* nrootx = no. de autovetores desejados, */
				/* em ordem crescente de autovalores */
				/* njx = dimensao da matriz */
				/* a[] = matriz a ser diagonalizada, */
				/* guardada na forma triangular inferior: */
				/* a[1] */
				/* a[2] a[3] */
				/* a[4] a[5] a[6] ... ; a diagonalizacao */
				/* destroi a matriz a */
				/* root[] armazena autovalores em ordem */
				/* crescente root[1] root[2] */
				/* vect[][] armazena autovetores */
				/* vect[1][1] vect[1][2] ... vect[1][n] */
				/* vect[2][1] vect[2][2] ... */
				/* ... */
				/* vect[nroot][1] vect[nroot][2] ...*/
extern int givens(int n, int nrootx, int njx, double *a, double *root,
	   double *vect[]);

				/* n = dimensao da matriz */
				/* a[] = mesma de acima */
				/* c[] = ponteiro para matriz onde sera' */
				/* guardada copia de a[], pre-alocada */
				/* root[] = mesma de acima */
				/* vect[][] = mesma de acima */
				/* nrott = mesmo de acima */
				/* chk = 0 para copiar a em c */
				/* chk = 1 para checar  */
extern int check(int n, double *a, double *c, double *root,
		       double* vect[], int nroot, int chk);
#endif
