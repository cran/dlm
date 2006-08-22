# include <stdio.h>
# include <R.h>
# include <Rinternals.h>

static double sqrarg;
# define SQR(x) ((sqrarg = (x)) == 0.0 ? 0.0 : sqrarg*sqrarg)

SEXP dlmLL();
SEXP dlmLL0(); /* for time-invariant models */
SEXP dlmFilter();
SEXP dlmFilter0(); /* for time-invariant models */
SEXP dlmSmooth();
SEXP dlmSmooth0(); /* for time-invariant models */
SEXP dlmForecast();
SEXP ARtranspar();

void F77_NAME(dgesdd)(const char *jobz,
                      const int *m, const int *n,
                      double *a, const int *lda, double *s,
                      double *u, const int *ldu,
                      double *vt, const int *ldvt,
                      double *work, const int *lwork, int *iwork, int *info);

void pmatrix();
void pIntMatrix();

void pmatrix(char *txt, double *x, /* print matrix, for debugging */
             int ldx, int m, int n)
{
    int i, j;
    printf("%s\n",txt);
    for (i=0; i<m; i++)
    {
        for (j=0; j<n; j++)
            printf("%6.5lg ",x[i+j*ldx]);
        printf("\n");
    }
    printf("\n");
}

void pIntMatrix(char *txt, int *x, /* print int matrix, for debugging */
             int ldx, int m, int n)
{
    int i, j;
    printf("%s\n",txt);
    for (i=0; i<m; i++)
    {
        for (j=0; j<n; j++)
            printf("%d ",x[i+j*ldx]);
        printf("\n");
    }
    printf("\n");
}

SEXP dlmLL(SEXP y, SEXP mod, SEXP tvFF, SEXP tvV, SEXP tvGG, SEXP tvW)
{
/***** Warning: the function relies on the order of the  *****/
/***** components of the list 'mod', not on their names. *****/     

    SEXP val;
    int i, j, k, l, p, m, n, t, max_m_p, la_m, la_n, la_info=0, la_lwork, *la_iwork;
    int stvFF=INTEGER(tvFF)[0], stvV=INTEGER(tvV)[0], stvGG=INTEGER(tvGG)[0], 
	stvW=INTEGER(tvW)[0], stvFV, *sJFF, *sJV, *sJGG, *sJW, nrJFF, nrJV, nrJGG, nrJW;
    double *sy=REAL(y), *sm0, *sFF, *sV, *sGG, *sW, *sX, *Ux,  
        *Dx, *sqrtV, *sqrtW, 
        *sqrtVinv, *a, *Ux_prior, *Dx_prior, *f, *Uy, *Dy,
        *e, *tF_Vinv, ll=0.0;
    double tmp, tmp1, *tmpMat, *tmpMat2, *la_s, *la_u, *la_vt, *la_work;
    char la_jobz='S';

    PROTECT(val = allocVector(REALSXP, 1));
    m = INTEGER(getAttrib( VECTOR_ELT(mod,2), R_DimSymbol ))[0];
    p = INTEGER(getAttrib( VECTOR_ELT(mod,2), R_DimSymbol ))[1];
    max_m_p = m > p ? m : p;
    la_n = max_m_p;
    la_m = la_n + p;
    n = length(y) / m; 
    sm0 = (double *) R_alloc( p, sizeof(double) );
    for (i = 0; i < p; i++)
        sm0[i] = REAL(VECTOR_ELT(mod,0))[i];
    sFF = REAL(VECTOR_ELT(mod,2));
    sV = REAL(VECTOR_ELT(mod,3));
    sGG = REAL(VECTOR_ELT(mod,4));
    sW = REAL(VECTOR_ELT(mod,5));
    sqrtV = (double *) R_alloc( m * m, sizeof(double) );
    sqrtVinv = (double *) R_alloc( m * m, sizeof(double) );
    sqrtW = (double *) R_alloc( p * p, sizeof(double) );
    tF_Vinv = (double *) R_alloc( p * m, sizeof(double) );
    Ux = (double *) R_alloc( p * p, sizeof(double) ); 
    Dx = (double *) R_alloc( p, sizeof(double) ); 
    Ux_prior = (double *) R_alloc( p * p, sizeof(double) ); 
    Dx_prior = (double *) R_alloc( p, sizeof(double) ); 
    a = (double *) R_alloc( p, sizeof(double) ); 
    f = (double *) R_alloc( m, sizeof(double) ); 
    Uy = (double *) R_alloc( m * m, sizeof(double) );
    Dy = (double *) R_alloc( m, sizeof(double) ); 
    e = (double *) R_alloc( m, sizeof(double) );
    tmpMat2 = (double *) R_alloc( m * p, sizeof(double) );

    /* allocate space for a la_m by la_n matrix */
    tmpMat = (double *) R_alloc( la_m * la_n, sizeof(double) ); 
    /* space for singular values */
    la_s = (double *) R_alloc( max_m_p, sizeof(double) );
    /* space for U matrix and Vt matrix (singular vectors) */
    la_u = (double *) R_alloc( la_m * max_m_p, sizeof(double) );
    la_vt = (double *) R_alloc( max_m_p * max_m_p, sizeof(double) );
    /* space for la_iwork */
    la_iwork = (int *) R_alloc( 8 * p, sizeof(int) );

    /* ask for optimal size of work array */
    la_lwork = -1;
    F77_CALL(dgesdd)(&la_jobz,
                     &la_m, &la_n, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     &tmp, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    la_lwork = (int) tmp;
    la_work = (double *) R_alloc( la_lwork, sizeof(double) );

    /** preliminaries: compute svd of C0, etc... **/
    for (i = 0; i < p; i++) {
        for (j = 0; j<i; j++) 
            tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
                REAL(VECTOR_ELT(mod,1))[i + p * j]; /* C0 */
        tmpMat[i + la_m * j] = REAL(VECTOR_ELT(mod,1))[i + p * j];
    }
    F77_CALL(dgesdd)(&la_jobz,
                     &p, &p, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     la_work, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    for (i = 0; i < p; i++) {
        for (j = 0; j < p; j++)
            Ux[i + j * p] = la_vt[j + i * max_m_p];
        Dx[i] = sqrt( la_s[i] );
    }
    /* time-varying matrices stuff */
    stvFV = stvFF || stvV;
    sX = REAL(VECTOR_ELT(mod,10));
    if (stvFF) {
	sJFF=INTEGER(VECTOR_ELT(mod,6));
	nrJFF=INTEGER(getAttrib( VECTOR_ELT(mod,6), R_DimSymbol ))[0];
    }
    if (stvV) {
	sJV=INTEGER(VECTOR_ELT(mod,7));
	nrJV=INTEGER(getAttrib( VECTOR_ELT(mod,7), R_DimSymbol ))[0];
    } else {
	/* compute sqrtV, time-invariant*/
	for (i = 0; i < m; i++) {
	    for (j = 0; j<i; j++) 
		tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
		    sV[i + m * j];
	    tmpMat[i + la_m * j] = sV[i + m * j];
	}
	F77_CALL(dgesdd)(&la_jobz,
			 &m, &m, tmpMat, &la_m, la_s,
			 la_u, &la_m,
			 la_vt, &max_m_p,
			 la_work, &la_lwork, la_iwork, &la_info);
	if (la_info != 0)
	    error("error code %d from Lapack routine dgesdd", la_info);
	for (i = 0; i < m; i++) {
	    tmp = sqrt( la_s[i] );
	    tmp1 = 1 / tmp;
	    tmp1 = R_FINITE(tmp1) ? tmp1 : 0.0;
	    for (j = 0; j<m; j++) {
		sqrtV[i + j * m] = tmp * la_vt[i + j * max_m_p];
		sqrtVinv[i + j * m] = tmp1 * la_vt[i + j * max_m_p];
	    }
	}
	if (!stvFF) {
	    /* compute also tF_Vinv, time-invariant */
	    for (i = 0; i < m; i++) {
		for (j = 0; j < i; j++) {
		    tmp = 0.0;
		    for (k = 0; k < m; k++)
			tmp += sqrtVinv[k + i * m] * sqrtVinv[k + j * m];
		    tmpMat[i + j * la_m] = tmpMat[j + i * la_m] = tmp;
		}
		tmp = 0.0;
		for (k = 0; k < m; k++)
		    tmp += SQR( sqrtVinv[k + i * m] );
		tmpMat[i + i * la_m] = tmp;
	    }
	    for (i = 0; i < p ; i++) 
		for (j = 0; j < m; j++) {
		    tmp = 0.0;
		    for (k = 0; k < m; k++)
			tmp += sFF[k + i * m] * tmpMat[k + j * la_m];
		    tF_Vinv[i + j * p] = tmp;
		}
	}
    }
    if (stvGG) {
	sJGG=INTEGER(VECTOR_ELT(mod,8));
	nrJGG=INTEGER(getAttrib( VECTOR_ELT(mod,8), R_DimSymbol ))[0];
    }
    if (stvW) {
	sJW=INTEGER(VECTOR_ELT(mod,9));
	nrJW=INTEGER(getAttrib( VECTOR_ELT(mod,9), R_DimSymbol ))[0];
    } else {
	/* compute sqrtW, time-invariant */
	for (i = 0; i < p; i++) {
	    for (j = 0; j < i; j++) 
		tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
		    sW[i + p * j];
	    tmpMat[i + la_m * j] = sW[i + p * j];
	}
	F77_CALL(dgesdd)(&la_jobz,
			 &p, &p, tmpMat, &la_m, la_s,
			 la_u, &la_m,
			 la_vt, &max_m_p,
			 la_work, &la_lwork, la_iwork, &la_info);
	if (la_info != 0)
	    error("error code %d from Lapack routine dgesdd", la_info);
	for (i = 0; i < p; i++) {
	    tmp = sqrt( la_s[i] );
	    for (j = 0; j < p; j++) 
		sqrtW[i + j * p] = tmp * la_vt[i + j * max_m_p];
	}
    }
    
/*      pmatrix("sqrtW",sqrtW,p,p,p); pmatrix("sqrtV",sqrtV,m,m,m); */
/*      pmatrix("sqrtVinv",sqrtVinv,m,m,m);pmatrix("tF_Vinv",tF_Vinv,p,p,m);error(""); */

    /** loop over observations **/
    for (t = 0; t < n; t++) { 
	/** set time-varying matrices **/
	if (stvFF) 
	    for (i = 0; i < nrJFF; i++)
		sFF[ sJFF[i] + m * sJFF[i + nrJFF] ] = sX[ t + n * sJFF[i + 2 * nrJFF] ];
	if (stvV) {
	    for (i = 0; i < nrJV; i++)
		sV[ sJV[i] + m * sJV[i + nrJV] ] = sX[ t + n * sJV[i + 2 * nrJV] ];
	    for (i = 0; i < m; i++) {
		for (j = 0; j<i; j++) 
		    tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
			sV[i + m * j];
		tmpMat[i + la_m * j] = sV[i + m * j];
	    }
	    F77_CALL(dgesdd)(&la_jobz,
			     &m, &m, tmpMat, &la_m, la_s,
			     la_u, &la_m,
			     la_vt, &max_m_p,
			     la_work, &la_lwork, la_iwork, &la_info);
	    if (la_info != 0)
		error("error code %d from Lapack routine dgesdd", la_info);
	    for (i = 0; i < m; i++) {
		tmp = sqrt( la_s[i] );
		tmp1 = 1 / tmp;
		tmp1 = R_FINITE(tmp1) ? tmp1 : 0.0;
		for (j = 0; j<m; j++) {
		    sqrtV[i + j * m] = tmp * la_vt[i + j * max_m_p];
		    sqrtVinv[i + j * m] = tmp1 * la_vt[i + j * max_m_p];
		}
	    }
	}
	if (stvGG) 
	    for (i = 0; i < nrJGG; i++)
		sGG[ sJGG[i] + p * sJGG[i + nrJGG] ] = sX[ t + n * sJGG[i + 2 * nrJGG] ];
	if (stvW) {
	    for (i = 0; i < nrJW; i++)
		sW[ sJW[i] + p * sJW[i + nrJW] ] = sX[ t + n * sJW[i + 2 * nrJW] ];
	    for (i = 0; i < p; i++) {
		for (j = 0; j < i; j++) 
		    tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
			sW[i + p * j];
		tmpMat[i + la_m * j] = sW[i + p * j];
	    }
	    F77_CALL(dgesdd)(&la_jobz,
			     &p, &p, tmpMat, &la_m, la_s,
			     la_u, &la_m,
			     la_vt, &max_m_p,
			     la_work, &la_lwork, la_iwork, &la_info);
	    if (la_info != 0)
		error("error code %d from Lapack routine dgesdd", la_info);
	    for (i = 0; i < p; i++) {
		tmp = sqrt( la_s[i] );
		for (j = 0; j < p; j++) 
		    sqrtW[i + j * p] = tmp * la_vt[i + j * max_m_p];
	    }
	}
	if (stvFV) {
	    for (i = 0; i < m; i++) {
		for (j = 0; j < i; j++) {
		    tmp = 0.0;
		    for (k = 0; k < m; k++)
			tmp += sqrtVinv[k + i * m] * sqrtVinv[k + j * m];
		    tmpMat[i + j * la_m] = tmpMat[j + i * la_m] = tmp;
		}
		tmp = 0.0;
		for (k = 0; k < m; k++)
		    tmp += SQR( sqrtVinv[k + i * m] );
		tmpMat[i + i * la_m] = tmp;
	    }
	    for (i = 0; i < p ; i++) 
		for (j = 0; j < m; j++) {
		    tmp = 0.0;
		    for (k = 0; k < m; k++)
			tmp += sFF[k + i * m] * tmpMat[k + j * la_m];
		    tF_Vinv[i + j * p] = tmp;
		}
	}
        /** Prior **/
        for (i = 0; i < p; i++) {
            tmp = 0.0;
            for (k = 0; k < p; k++)
                tmp += sGG[i + p * k] * sm0[k];
            a[i] = tmp;
        }
        
        for (i = 0; i < p; i++) {
            tmp1 = Dx[i];
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (l = 0; l < p; l++)
                    tmp += sGG[j + l * p] * Ux[i * p + l];
                tmpMat[i + j * la_m] = tmp * tmp1;
            }
            
        }
        for (i = 0; i < p; i++) 
            for (j = 0; j < p; j++) 
                tmpMat[i + p + j * la_m] = sqrtW[i + j * p];
        l = 2 * p;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &p, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < p; i++) {
            for (j = 0; j < p; j++)
                Ux_prior[i + j * p] = la_vt[j + i * max_m_p];
            Dx_prior[i] = la_s[i];
        }
/*	pmatrix("Ux_prior",Ux_prior,p,p,p); pmatrix("Dx_prior",Dx_prior,1,1,p); */
/*          pmatrix("a",a,1,1,p); */

        /** One-step forecast **/
        for (i = 0; i < m; i++) {
            tmp = 0.0;
            for (j = 0; j < p; j++)
                tmp += sFF[i + j * m] * a[j];
            f[i] = tmp;
        }
        for (i = 0; i < p; i++) {
            tmp1 = Dx_prior[i];
            for (j = 0; j < m; j++) {
                tmp = 0.0;
                for (l = 0; l < p; l++)
                    tmp += sFF[j + l * m] * Ux_prior[l + i * p];
                tmpMat[i + j * la_m] = tmp * tmp1;
            }
            
        }
        for (i = 0; i < m; i++) 
            for (j = 0; j < m; j++) 
                tmpMat[i + p + j * la_m] = sqrtV[i + j * m];
        
        l = p + m;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &m, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++)
                Uy[i + j * m] = la_vt[j + i * max_m_p];
            Dy[i] = la_s[i];
        }
/*          pmatrix("Uy",Uy,m,m,m); pmatrix("Dy",Dy,1,1,m);  */
/*          pmatrix("f",f,1,1,m); error("");  */
        
        /** Posterior **/
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (k = 0; k < m; k++)
                    tmp += sqrtVinv[i + k * m] * sFF[k + j * m];
                tmpMat2[i + j * m] = tmp;
            }
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += tmpMat2[i + k * m] * Ux_prior[k + j * p];
                tmpMat[i + j * la_m] = tmp;
            }
        for (i = 0; i < p; i++) {
            tmp = 1 / Dx_prior[i];
            tmpMat[i + m + i * la_m] = R_FINITE(tmp) ? tmp : 0.0;
            for (j = i + 1; j < p; j++) 
                tmpMat[j + m + i * la_m] = tmpMat[i + m + j * la_m] = 0.0;
        }
        l = p + m;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &p, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < p; i++) {
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += Ux_prior[i + k * p] * la_vt[j + k * max_m_p];
                Ux[i + j * p] = tmp;
            }
            tmp = 1 / la_s[i];
            Dx[i] = R_FINITE(tmp) ? tmp : 0.0;
        }
        for (i = 0; i < m; i++)
            e[i] = sy[t + i * n] - f[i];
	for (i = 0; i < p; i++)
            for (j = 0; j < m; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += SQR(Dx[i]) * Ux[k + i * p] * tF_Vinv[k + j * p];
                tmpMat[i + j * la_m] = tmp;
            }
        for (i = 0; i < p; i++)
            for (j = 0; j < m; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += Ux[i + k * p] * tmpMat[k + j * la_m];
                tmpMat2[j + i * m] = tmp;
            }
        for (i = 0; i < p; i++) {
            tmp = a[i];
            for (j = 0; j < m; j++)
                tmp += tmpMat2[j + i * m] * e[j];
            sm0[i] = tmp;
        }
/*          pmatrix("Ux",Ux,p,p,p); pmatrix("Dx",Dx, 1,1,p);  */
/*          pmatrix("m0",sm0,1,1,p); error(""); */

        /** Update negative log likelihood **/
        for (i = 0; i < m; i++) {
            tmp =0.0;
            for (j = 0; j < m; j++)
                tmp += Uy[j + i * m] * e[j];
            ll += SQR( tmp / Dy[i] ) + 2.0 * log( Dy[i] );
        }
    }

    REAL(val)[0] = 0.5 * ll;
    UNPROTECT(1);
    return(val);

}


SEXP dlmLL0(SEXP y, SEXP mod) 
{
/***** Warning: the function relies on the order of the  *****/
/***** components of the list 'mod', not on their names. *****/     

    SEXP val;
    int i, j, k, l, p, m, n, t, max_m_p, la_m, la_n, la_info=0, la_lwork, *la_iwork;
    double *sy=REAL(y), *sm0, *sFF, *sGG, *Ux,  
        *Dx, *sqrtV, *sqrtW, 
        *sqrtVinv, *a, *Ux_prior, *Dx_prior, *f, *Uy, *Dy,
        *e, *tF_Vinv, ll=0.0;
    double tmp, tmp1, *tmpMat, *tmpMat2, *la_s, *la_u, *la_vt, *la_work;
    char la_jobz='S';

    PROTECT(val = allocVector(REALSXP, 1));
    m = INTEGER(getAttrib( VECTOR_ELT(mod,2), R_DimSymbol ))[0];
    p = INTEGER(getAttrib( VECTOR_ELT(mod,2), R_DimSymbol ))[1];
    max_m_p = m > p ? m : p;
    la_n = max_m_p;
    la_m = la_n + p;
    n = length(y) / m; 
    sm0 = (double *) R_alloc( p, sizeof(double) );
    for (i = 0; i < p; i++)
        sm0[i] = REAL(VECTOR_ELT(mod,0))[i];
    sFF = REAL(VECTOR_ELT(mod,2));
    sGG = REAL(VECTOR_ELT(mod,4));
    sqrtV = (double *) R_alloc( m * m, sizeof(double) );
    sqrtVinv = (double *) R_alloc( m * m, sizeof(double) );
    sqrtW = (double *) R_alloc( p * p, sizeof(double) );
    tF_Vinv = (double *) R_alloc( p * m, sizeof(double) );
    Ux = (double *) R_alloc( p * p, sizeof(double) ); 
    Dx = (double *) R_alloc( p, sizeof(double) ); 
    Ux_prior = (double *) R_alloc( p * p, sizeof(double) ); 
    Dx_prior = (double *) R_alloc( p, sizeof(double) ); 
    a = (double *) R_alloc( p, sizeof(double) ); 
    f = (double *) R_alloc( m, sizeof(double) ); 
    Uy = (double *) R_alloc( m * m, sizeof(double) );
    Dy = (double *) R_alloc( m, sizeof(double) ); 
    e = (double *) R_alloc( m, sizeof(double) );
    tmpMat2 = (double *) R_alloc( m * p, sizeof(double) );

    /* allocate space for a la_m by la_n matrix */
    tmpMat = (double *) R_alloc( la_m * la_n, sizeof(double) ); 
    /* space for singular values */
    la_s = (double *) R_alloc( max_m_p, sizeof(double) );
    /* space for U matrix and Vt matrix (singular vectors) */
    la_u = (double *) R_alloc( la_m * max_m_p, sizeof(double) );
    la_vt = (double *) R_alloc( max_m_p * max_m_p, sizeof(double) );
    /* space for la_iwork */
    la_iwork = (int *) R_alloc( 8 * p, sizeof(int) );

    /* ask for optimal size of work array */
    la_lwork = -1;
    F77_CALL(dgesdd)(&la_jobz,
                     &la_m, &la_n, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     &tmp, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    la_lwork = (int) tmp;
    la_work = (double *) R_alloc( la_lwork, sizeof(double) );

    /** preliminaries: compute svd of C0, sqrt(V), sqrt(W), etc... **/
    for (i = 0; i < p; i++) {
        for (j = 0; j<i; j++) 
            tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
                REAL(VECTOR_ELT(mod,1))[i + p * j]; /* C0 */
        tmpMat[i + la_m * j] = REAL(VECTOR_ELT(mod,1))[i + p * j];
    }
    F77_CALL(dgesdd)(&la_jobz,
                     &p, &p, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     la_work, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    for (i = 0; i < p; i++) {
        for (j = 0; j < p; j++)
            Ux[i + j * p] = la_vt[j + i * max_m_p];
        Dx[i] = sqrt( la_s[i] );
    }
    for (i = 0; i < p; i++) {
        for (j = 0; j < i; j++) 
            tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
                REAL(VECTOR_ELT(mod,5))[i + p * j]; /* W */
        tmpMat[i + la_m * j] = REAL(VECTOR_ELT(mod,5))[i + p * j];
    }
    F77_CALL(dgesdd)(&la_jobz,
                     &p, &p, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     la_work, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    for (i = 0; i < p; i++) {
        tmp = sqrt( la_s[i] );
        for (j = 0; j < p; j++) 
            sqrtW[i + j * p] = tmp * la_vt[i + j * max_m_p];
    }
    for (i = 0; i < m; i++) {
        for (j = 0; j<i; j++) 
            tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
                REAL(VECTOR_ELT(mod,3))[i + m * j]; /* V */
        tmpMat[i + la_m * j] = REAL(VECTOR_ELT(mod,3))[i + m * j];
    }
    F77_CALL(dgesdd)(&la_jobz,
                     &m, &m, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     la_work, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    for (i = 0; i < m; i++) {
        tmp = sqrt( la_s[i] );
        tmp1 = 1 / tmp;
        tmp1 = R_FINITE(tmp1) ? tmp1 : 0.0;
        for (j = 0; j<m; j++) {
            sqrtV[i + j * m] = tmp * la_vt[i + j * max_m_p];
            sqrtVinv[i + j * m] = tmp1 * la_vt[i + j * max_m_p];
        }
    }
    for (i = 0; i < m; i++) {
        for (j = 0; j < i; j++) {
            tmp = 0.0;
            for (k = 0; k < m; k++)
                tmp += sqrtVinv[k + i * m] * sqrtVinv[k + j * m];
            tmpMat[i + j * la_m] = tmpMat[j + i * la_m] = tmp;
        }
        tmp = 0.0;
        for (k = 0; k < m; k++)
            tmp += SQR( sqrtVinv[k + i * m] );
        tmpMat[i + i * la_m] = tmp;
    }
    for (i = 0; i < p ; i++) 
        for (j = 0; j < m; j++) {
            tmp = 0.0;
            for (k = 0; k < m; k++)
                tmp += sFF[k + i * m] * tmpMat[k + j * la_m];
            tF_Vinv[i + j * p] = tmp;
        }

/*      pmatrix("sqrtW",sqrtW,p,p,p); pmatrix("sqrtV",sqrtV,m,m,m); */
/*      pmatrix("sqrtVinv",sqrtVinv,m,m,m);pmatrix("tF_Vinv",tF_Vinv,p,p,m);error(""); */

    /** loop over observations **/
    for (t = 0; t < n; t++) { 
        /** Prior **/
        for (i = 0; i < p; i++) {
            tmp = 0.0;
            for (k = 0; k < p; k++)
                tmp += sGG[i + p * k] * sm0[k];
            a[i] = tmp;
        }
        
        for (i = 0; i < p; i++) {
            tmp1 = Dx[i];
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (l = 0; l < p; l++)
                    tmp += sGG[j + l * p] * Ux[i * p + l];
                tmpMat[i + j * la_m] = tmp * tmp1;
            }
            
        }
        for (i = 0; i < p; i++) 
            for (j = 0; j < p; j++) 
                tmpMat[i + p + j * la_m] = sqrtW[i + j * p];
        l = 2 * p;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &p, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < p; i++) {
            for (j = 0; j < p; j++)
                Ux_prior[i + j * p] = la_vt[j + i * max_m_p];
            Dx_prior[i] = la_s[i];
        }
/*	pmatrix("Ux_prior",Ux_prior,p,p,p); pmatrix("Dx_prior",Dx_prior,1,1,p); */
/*          pmatrix("a",a,1,1,p); */

        /** One-step forecast **/
        for (i = 0; i < m; i++) {
            tmp = 0.0;
            for (j = 0; j < p; j++)
                tmp += sFF[i + j * m] * a[j];
            f[i] = tmp;
        }
        for (i = 0; i < p; i++) {
            tmp1 = Dx_prior[i];
            for (j = 0; j < m; j++) {
                tmp = 0.0;
                for (l = 0; l < p; l++)
                    tmp += sFF[j + l * m] * Ux_prior[l + i * p];
                tmpMat[i + j * la_m] = tmp * tmp1;
            }
            
        }
        for (i = 0; i < m; i++) 
            for (j = 0; j < m; j++) 
                tmpMat[i + p + j * la_m] = sqrtV[i + j * m];
        
        l = p + m;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &m, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++)
                Uy[i + j * m] = la_vt[j + i * max_m_p];
            Dy[i] = la_s[i];
        }
/*          pmatrix("Uy",Uy,m,m,m); pmatrix("Dy",Dy,1,1,m);  */
/*          pmatrix("f",f,1,1,m); error("");  */
        
        /** Posterior **/
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (k = 0; k < m; k++)
                    tmp += sqrtVinv[i + k * m] * sFF[k + j * m];
                tmpMat2[i + j * m] = tmp;
            }
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += tmpMat2[i + k * m] * Ux_prior[k + j * p];
                tmpMat[i + j * la_m] = tmp;
            }
        for (i = 0; i < p; i++) {
            tmp = 1 / Dx_prior[i];
            tmpMat[i + m + i * la_m] = R_FINITE(tmp) ? tmp : 0.0;
            for (j = i + 1; j < p; j++) 
                tmpMat[j + m + i * la_m] = tmpMat[i + m + j * la_m] = 0.0;
        }
        l = p + m;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &p, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < p; i++) {
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += Ux_prior[i + k * p] * la_vt[j + k * max_m_p];
                Ux[i + j * p] = tmp;
            }
            tmp = 1 / la_s[i];
            Dx[i] = R_FINITE(tmp) ? tmp : 0.0;
        }
        for (i = 0; i < m; i++)
            e[i] = sy[t + i * n] - f[i];
	for (i = 0; i < p; i++)
            for (j = 0; j < m; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += SQR(Dx[i]) * Ux[k + i * p] * tF_Vinv[k + j * p];
                tmpMat[i + j * la_m] = tmp;
            }
        for (i = 0; i < p; i++)
            for (j = 0; j < m; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += Ux[i + k * p] * tmpMat[k + j * la_m];
                tmpMat2[j + i * m] = tmp;
            }
        for (i = 0; i < p; i++) {
            tmp = a[i];
            for (j = 0; j < m; j++)
                tmp += tmpMat2[j + i * m] * e[j];
            sm0[i] = tmp;
        }
/*          pmatrix("Ux",Ux,p,p,p); pmatrix("Dx",Dx, 1,1,p);  */
/*          pmatrix("m0",sm0,1,1,p); error(""); */

        /** Update negative log likelihood **/
        for (i = 0; i < m; i++) {
            tmp =0.0;
            for (j = 0; j < m; j++)
                tmp += Uy[j + i * m] * e[j];
            ll += SQR( tmp / Dy[i] ) + 2.0 * log( Dy[i] );
        }
    }

    REAL(val)[0] = 0.5 * ll;
    UNPROTECT(1);
    return(val);

}


SEXP dlmFilter(SEXP y, SEXP mod, SEXP tvFF, SEXP tvV, SEXP tvGG, SEXP tvW)
{
/***** Warning: the function relies on the order of the  *****/
/***** components of the list 'mod', not on their names. *****/     

    SEXP val, mR, UxR, DxR, aR, Ux_priorR, Dx_priorR, fR;
    int i, j, k, l, p, m, n, nPlus, t, max_m_p, la_m, la_n, la_info=0, la_lwork, *la_iwork;
    int stvFF=INTEGER(tvFF)[0], stvV=INTEGER(tvV)[0], stvGG=INTEGER(tvGG)[0], 
	stvW=INTEGER(tvW)[0], stvFV, *sJFF, *sJV, *sJGG, *sJW, nrJFF, nrJV, nrJGG, nrJW;
    double *sy=REAL(y), *sm0, *sFF, *sV, *sGG, *sW, *sX, *Ux,  
        *Dx, *sqrtV, *sqrtW, 
        *sqrtVinv, *a, *Ux_prior, *Dx_prior, *f, *Uy, *Dy,
        *e, *tF_Vinv;
    double tmp, tmp1, *tmpMat, *tmpMat2, *la_s, *la_u, *la_vt, *la_work;
    char la_jobz='S';

    m = INTEGER(getAttrib( VECTOR_ELT(mod,2), R_DimSymbol ))[0];
    p = INTEGER(getAttrib( VECTOR_ELT(mod,2), R_DimSymbol ))[1];
    max_m_p = m > p ? m : p;
    la_n = max_m_p;
    la_m = la_n + p;
    n = length(y) / m; 
    nPlus = n + 1;
    PROTECT(mR = allocMatrix(REALSXP, nPlus, p));
    PROTECT(UxR = allocVector(VECSXP, nPlus));
    PROTECT(DxR = allocMatrix(REALSXP, nPlus, p));
    PROTECT(aR = allocMatrix(REALSXP, n, p));
    PROTECT(Ux_priorR = allocVector(VECSXP, n));
    PROTECT(Dx_priorR = allocMatrix(REALSXP, n, p));
    PROTECT(fR = allocMatrix(REALSXP, n, m));
    for (i = 0; i < p; i++)
	REAL(mR)[i * nPlus] = REAL(VECTOR_ELT(mod,0))[i];
    sm0 = REAL(mR);
    SET_VECTOR_ELT(UxR, 0, allocMatrix(REALSXP, p, p));
    Ux = REAL(VECTOR_ELT(UxR, 0));
    Dx = REAL(DxR);
    a = REAL(aR);
    f = REAL(fR);
    Dx_prior = REAL(Dx_priorR);
    sFF = REAL(VECTOR_ELT(mod,2));
    sV = REAL(VECTOR_ELT(mod,3));
    sGG = REAL(VECTOR_ELT(mod,4));
    sW = REAL(VECTOR_ELT(mod,5));
    sqrtV = (double *) R_alloc( m * m, sizeof(double) );
    sqrtVinv = (double *) R_alloc( m * m, sizeof(double) );
    sqrtW = (double *) R_alloc( p * p, sizeof(double) );
    tF_Vinv = (double *) R_alloc( p * m, sizeof(double) );
    Uy = (double *) R_alloc( m * m, sizeof(double) );
    Dy = (double *) R_alloc( m, sizeof(double) ); 
    e = (double *) R_alloc( m, sizeof(double) );
    tmpMat2 = (double *) R_alloc( m * p, sizeof(double) );
    /* allocate space for a la_m by la_n matrix */
    tmpMat = (double *) R_alloc( la_m * la_n, sizeof(double) ); 
    /* space for singular values */
    la_s = (double *) R_alloc( max_m_p, sizeof(double) );
    /* space for U matrix and Vt matrix (singular vectors) */
    la_u = (double *) R_alloc( la_m * max_m_p, sizeof(double) );
    la_vt = (double *) R_alloc( max_m_p * max_m_p, sizeof(double) );
    /* space for la_iwork */
    la_iwork = (int *) R_alloc( 8 * p, sizeof(int) );

    /* ask for optimal size of work array */
    la_lwork = -1;
    F77_CALL(dgesdd)(&la_jobz,
                     &la_m, &la_n, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     &tmp, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    la_lwork = (int) tmp;
    la_work = (double *) R_alloc( la_lwork, sizeof(double) );

    /** preliminaries: compute svd of C0, etc... **/
    for (i = 0; i < p; i++) {
        for (j = 0; j<i; j++) 
            tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
                REAL(VECTOR_ELT(mod,1))[i + p * j]; /* C0 */
        tmpMat[i + la_m * j] = REAL(VECTOR_ELT(mod,1))[i + p * j];
    }
    F77_CALL(dgesdd)(&la_jobz,
                     &p, &p, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     la_work, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    for (i = 0; i < p; i++) {
        for (j = 0; j < p; j++)
            Ux[i + j * p] = la_vt[j + i * max_m_p];
        Dx[i * nPlus] = sqrt( la_s[i] );
    }
    /* time-varying matrices stuff */
    stvFV = stvFF || stvV;
    sX = REAL(VECTOR_ELT(mod,10));
    if (stvFF) {
	sJFF=INTEGER(VECTOR_ELT(mod,6));
	nrJFF=INTEGER(getAttrib( VECTOR_ELT(mod,6), R_DimSymbol ))[0];
    }
    if (stvV) {
	sJV=INTEGER(VECTOR_ELT(mod,7));
	nrJV=INTEGER(getAttrib( VECTOR_ELT(mod,7), R_DimSymbol ))[0];
    } else {
	/* compute sqrtV, time-invariant*/
	for (i = 0; i < m; i++) {
	    for (j = 0; j<i; j++) 
		tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
		    sV[i + m * j];
	    tmpMat[i + la_m * j] = sV[i + m * j];
	}
	F77_CALL(dgesdd)(&la_jobz,
			 &m, &m, tmpMat, &la_m, la_s,
			 la_u, &la_m,
			 la_vt, &max_m_p,
			 la_work, &la_lwork, la_iwork, &la_info);
	if (la_info != 0)
	    error("error code %d from Lapack routine dgesdd", la_info);
	for (i = 0; i < m; i++) {
	    tmp = sqrt( la_s[i] );
	    tmp1 = 1 / tmp;
	    tmp1 = R_FINITE(tmp1) ? tmp1 : 0.0;
	    for (j = 0; j<m; j++) {
		sqrtV[i + j * m] = tmp * la_vt[i + j * max_m_p];
		sqrtVinv[i + j * m] = tmp1 * la_vt[i + j * max_m_p];
	    }
	}
	if (!stvFF) {
	    /* compute also tF_Vinv, time-invariant */
	    for (i = 0; i < m; i++) {
		for (j = 0; j < i; j++) {
		    tmp = 0.0;
		    for (k = 0; k < m; k++)
			tmp += sqrtVinv[k + i * m] * sqrtVinv[k + j * m];
		    tmpMat[i + j * la_m] = tmpMat[j + i * la_m] = tmp;
		}
		tmp = 0.0;
		for (k = 0; k < m; k++)
		    tmp += SQR( sqrtVinv[k + i * m] );
		tmpMat[i + i * la_m] = tmp;
	    }
	    for (i = 0; i < p ; i++) 
		for (j = 0; j < m; j++) {
		    tmp = 0.0;
		    for (k = 0; k < m; k++)
			tmp += sFF[k + i * m] * tmpMat[k + j * la_m];
		    tF_Vinv[i + j * p] = tmp;
		}
	}
    }
    if (stvGG) {
	sJGG=INTEGER(VECTOR_ELT(mod,8));
	nrJGG=INTEGER(getAttrib( VECTOR_ELT(mod,8), R_DimSymbol ))[0];
    }
    if (stvW) {
	sJW=INTEGER(VECTOR_ELT(mod,9));
	nrJW=INTEGER(getAttrib( VECTOR_ELT(mod,9), R_DimSymbol ))[0];
    } else {
	/* compute sqrtW, time-invariant */
	for (i = 0; i < p; i++) {
	    for (j = 0; j < i; j++) 
		tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
		    sW[i + p * j];
	    tmpMat[i + la_m * j] = sW[i + p * j];
	}
	F77_CALL(dgesdd)(&la_jobz,
			 &p, &p, tmpMat, &la_m, la_s,
			 la_u, &la_m,
			 la_vt, &max_m_p,
			 la_work, &la_lwork, la_iwork, &la_info);
	if (la_info != 0)
	    error("error code %d from Lapack routine dgesdd", la_info);
	for (i = 0; i < p; i++) {
	    tmp = sqrt( la_s[i] );
	    for (j = 0; j < p; j++) 
		sqrtW[i + j * p] = tmp * la_vt[i + j * max_m_p];
	}
    }

    /** loop over observations **/
    for (t = 0; t < n; t++) { 
	/** allocate matrices and make pointers point to 'current' **/
	SET_VECTOR_ELT(Ux_priorR, t, allocMatrix(REALSXP, p, p));
	Ux_prior = REAL(VECTOR_ELT(Ux_priorR, t));
	SET_VECTOR_ELT(UxR, t+1, allocMatrix(REALSXP, p, p));
	/** set time-varying matrices **/
	if (stvFF) 
	    for (i = 0; i < nrJFF; i++)
		sFF[ sJFF[i] + m * sJFF[i + nrJFF] ] = sX[ t + n * sJFF[i + 2 * nrJFF] ];
	if (stvV) {
	    for (i = 0; i < nrJV; i++)
		sV[ sJV[i] + m * sJV[i + nrJV] ] = sX[ t + n * sJV[i + 2 * nrJV] ];
	    for (i = 0; i < m; i++) {
		for (j = 0; j<i; j++) 
		    tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
			sV[i + m * j];
		tmpMat[i + la_m * j] = sV[i + m * j];
	    }
	    F77_CALL(dgesdd)(&la_jobz,
			     &m, &m, tmpMat, &la_m, la_s,
			     la_u, &la_m,
			     la_vt, &max_m_p,
			     la_work, &la_lwork, la_iwork, &la_info);
	    if (la_info != 0)
		error("error code %d from Lapack routine dgesdd", la_info);
	    for (i = 0; i < m; i++) {
		tmp = sqrt( la_s[i] );
		tmp1 = 1 / tmp;
		tmp1 = R_FINITE(tmp1) ? tmp1 : 0.0;
		for (j = 0; j<m; j++) {
		    sqrtV[i + j * m] = tmp * la_vt[i + j * max_m_p];
		    sqrtVinv[i + j * m] = tmp1 * la_vt[i + j * max_m_p];
		}
	    }
	}
	if (stvGG) 
	    for (i = 0; i < nrJGG; i++)
		sGG[ sJGG[i] + p * sJGG[i + nrJGG] ] = sX[ t + n * sJGG[i + 2 * nrJGG] ];
	if (stvW) {
	    for (i = 0; i < nrJW; i++)
		sW[ sJW[i] + p * sJW[i + nrJW] ] = sX[ t + n * sJW[i + 2 * nrJW] ];
	    for (i = 0; i < p; i++) {
		for (j = 0; j < i; j++) 
		    tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
			sW[i + p * j];
		tmpMat[i + la_m * j] = sW[i + p * j];
	    }
	    F77_CALL(dgesdd)(&la_jobz,
			     &p, &p, tmpMat, &la_m, la_s,
			     la_u, &la_m,
			     la_vt, &max_m_p,
			     la_work, &la_lwork, la_iwork, &la_info);
	    if (la_info != 0)
		error("error code %d from Lapack routine dgesdd", la_info);
	    for (i = 0; i < p; i++) {
		tmp = sqrt( la_s[i] );
		for (j = 0; j < p; j++) 
		    sqrtW[i + j * p] = tmp * la_vt[i + j * max_m_p];
	    }
	}
	if (stvFV) {
	    for (i = 0; i < m; i++) {
		for (j = 0; j < i; j++) {
		    tmp = 0.0;
		    for (k = 0; k < m; k++)
			tmp += sqrtVinv[k + i * m] * sqrtVinv[k + j * m];
		    tmpMat[i + j * la_m] = tmpMat[j + i * la_m] = tmp;
		}
		tmp = 0.0;
		for (k = 0; k < m; k++)
		    tmp += SQR( sqrtVinv[k + i * m] );
		tmpMat[i + i * la_m] = tmp;
	    }
	    for (i = 0; i < p ; i++) 
		for (j = 0; j < m; j++) {
		    tmp = 0.0;
		    for (k = 0; k < m; k++)
			tmp += sFF[k + i * m] * tmpMat[k + j * la_m];
		    tF_Vinv[i + j * p] = tmp;
		}
	}
        /** Prior **/
        for (i = 0; i < p; i++) {
            tmp = 0.0;
            for (k = 0; k < p; k++)
                tmp += sGG[i + p * k] * sm0[k * nPlus];
            a[i * n] = tmp;
        }
        
        for (i = 0; i < p; i++) {
            tmp1 = Dx[i * nPlus];
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (l = 0; l < p; l++)
                    tmp += sGG[j + l * p] * Ux[i * p + l];
                tmpMat[i + j * la_m] = tmp * tmp1;
            }
        }
        for (i = 0; i < p; i++) 
            for (j = 0; j < p; j++) 
                tmpMat[i + p + j * la_m] = sqrtW[i + j * p];
        l = 2 * p;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &p, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < p; i++) {
            for (j = 0; j < p; j++)
                Ux_prior[i + j * p] = la_vt[j + i * max_m_p];
            Dx_prior[i * n] = la_s[i];
        }

        /** One-step forecast **/
        for (i = 0; i < m; i++) {
            tmp = 0.0;
            for (j = 0; j < p; j++)
                tmp += sFF[i + j * m] * a[j * n];
            f[i * n] = tmp;
        }
        for (i = 0; i < p; i++) {
            tmp1 = Dx_prior[i * n];
            for (j = 0; j < m; j++) {
                tmp = 0.0;
                for (l = 0; l < p; l++)
                    tmp += sFF[j + l * m] * Ux_prior[l + i * p];
                tmpMat[i + j * la_m] = tmp * tmp1;
            }
        }
        for (i = 0; i < m; i++) 
            for (j = 0; j < m; j++) 
                tmpMat[i + p + j * la_m] = sqrtV[i + j * m];
        
        l = p + m;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &m, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++)
                Uy[i + j * m] = la_vt[j + i * max_m_p];
            Dy[i] = la_s[i];
        }
        
        /** Posterior **/
        sm0++; Dx++;
	Ux = REAL(VECTOR_ELT(UxR, t+1));
	for (i = 0; i < m; i++)
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (k = 0; k < m; k++)
                    tmp += sqrtVinv[i + k * m] * sFF[k + j * m];
                tmpMat2[i + j * m] = tmp;
            }
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += tmpMat2[i + k * m] * Ux_prior[k + j * p];
                tmpMat[i + j * la_m] = tmp;
            }
        for (i = 0; i < p; i++) {
            tmp = 1 / Dx_prior[i * n];
            tmpMat[i + m + i * la_m] = R_FINITE(tmp) ? tmp : 0.0;
            for (j = i + 1; j < p; j++) 
                tmpMat[j + m + i * la_m] = tmpMat[i + m + j * la_m] = 0.0;
        }
        l = p + m;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &p, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < p; i++) {
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += Ux_prior[i + k * p] * la_vt[j + k * max_m_p];
                Ux[i + j * p] = tmp;
            }
            tmp = 1 / la_s[i];
            Dx[i * nPlus] = R_FINITE(tmp) ? tmp : 0.0;
        }
        for (i = 0; i < m; i++)
            e[i] = sy[i * n] - f[i * n];
	for (i = 0; i < p; i++)
            for (j = 0; j < m; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += SQR(Dx[i * nPlus]) * Ux[k + i * p] * tF_Vinv[k + j * p];
                tmpMat[i + j * la_m] = tmp;
            }
        for (i = 0; i < p; i++)
            for (j = 0; j < m; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += Ux[i + k * p] * tmpMat[k + j * la_m];
                tmpMat2[j + i * m] = tmp;
            }
        for (i = 0; i < p; i++) {
            tmp = a[i * n];
            for (j = 0; j < m; j++)
                tmp += tmpMat2[j + i * m] * e[j];
            sm0[i * nPlus] = tmp;
        }

	/** increment pointers **/
	sy++; a++; Dx_prior++, f++;
    }

    PROTECT(val = allocVector(VECSXP, 7));
    SET_VECTOR_ELT(val, 0, mR);
    SET_VECTOR_ELT(val, 1, UxR);
    SET_VECTOR_ELT(val, 2, DxR);
    SET_VECTOR_ELT(val, 3, aR);
    SET_VECTOR_ELT(val, 4, Ux_priorR);
    SET_VECTOR_ELT(val, 5, Dx_priorR);
    SET_VECTOR_ELT(val, 6, fR);
    UNPROTECT(8);
    
    return(val);
}

SEXP dlmFilter0(SEXP y, SEXP mod)
{
/***** Warning: the function relies on the order of the  *****/
/***** components of the list 'mod', not on their names. *****/     

    SEXP val, mR, UxR, DxR, aR, Ux_priorR, Dx_priorR, fR;
    int i, j, k, l, p, m, n, nPlus, t, max_m_p, la_m, la_n, la_info=0, la_lwork, *la_iwork;
    double *sy=REAL(y), *sm0, *sFF, *sGG, *Ux,  
        *Dx, *sqrtV, *sqrtW, 
        *sqrtVinv, *a, *Ux_prior, *Dx_prior, *f, *Uy, *Dy,
        *e, *tF_Vinv;
    double tmp, tmp1, *tmpMat, *tmpMat2, *la_s, *la_u, *la_vt, *la_work;
    char la_jobz='S';

    m = INTEGER(getAttrib( VECTOR_ELT(mod,2), R_DimSymbol ))[0];
    p = INTEGER(getAttrib( VECTOR_ELT(mod,2), R_DimSymbol ))[1];
    max_m_p = m > p ? m : p;
    la_n = max_m_p;
    la_m = la_n + p;
    n = length(y) / m; 
    nPlus = n + 1;
    PROTECT(mR = allocMatrix(REALSXP, nPlus, p));
    PROTECT(UxR = allocVector(VECSXP, nPlus));
    PROTECT(DxR = allocMatrix(REALSXP, nPlus, p));
    PROTECT(aR = allocMatrix(REALSXP, n, p));
    PROTECT(Ux_priorR = allocVector(VECSXP, n));
    PROTECT(Dx_priorR = allocMatrix(REALSXP, n, p));
    PROTECT(fR = allocMatrix(REALSXP, n, m));
    for (i = 0; i < p; i++)
	REAL(mR)[i * nPlus] = REAL(VECTOR_ELT(mod,0))[i];
    sm0 = REAL(mR);
    SET_VECTOR_ELT(UxR, 0, allocMatrix(REALSXP, p, p));
    Ux = REAL(VECTOR_ELT(UxR, 0));
    Dx = REAL(DxR);
    a = REAL(aR);
    f = REAL(fR);
    Dx_prior = REAL(Dx_priorR);
    sFF = REAL(VECTOR_ELT(mod,2));
    sGG = REAL(VECTOR_ELT(mod,4));
    sqrtV = (double *) R_alloc( m * m, sizeof(double) );
    sqrtVinv = (double *) R_alloc( m * m, sizeof(double) );
    sqrtW = (double *) R_alloc( p * p, sizeof(double) );
    tF_Vinv = (double *) R_alloc( p * m, sizeof(double) );
    Uy = (double *) R_alloc( m * m, sizeof(double) );
    Dy = (double *) R_alloc( m, sizeof(double) ); 
    e = (double *) R_alloc( m, sizeof(double) );
    tmpMat2 = (double *) R_alloc( m * p, sizeof(double) );

    /* allocate space for a la_m by la_n matrix */
    tmpMat = (double *) R_alloc( la_m * la_n, sizeof(double) ); 
    /* space for singular values */
    la_s = (double *) R_alloc( max_m_p, sizeof(double) );
    /* space for U matrix and Vt matrix (singular vectors) */
    la_u = (double *) R_alloc( la_m * max_m_p, sizeof(double) );
    la_vt = (double *) R_alloc( max_m_p * max_m_p, sizeof(double) );
    /* space for la_iwork */
    la_iwork = (int *) R_alloc( 8 * p, sizeof(int) );

    /* ask for optimal size of work array */
    la_lwork = -1;
    F77_CALL(dgesdd)(&la_jobz,
                     &la_m, &la_n, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     &tmp, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    la_lwork = (int) tmp;
    la_work = (double *) R_alloc( la_lwork, sizeof(double) );

    /** preliminaries: compute svd of C0, sqrt(V), sqrt(W), etc... **/
    for (i = 0; i < p; i++) {
	for (j = 0; j<i; j++) 
            tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
                REAL(VECTOR_ELT(mod,1))[i + p * j]; /* C0 */
        tmpMat[i + la_m * j] = REAL(VECTOR_ELT(mod,1))[i + p * j];
    }
    F77_CALL(dgesdd)(&la_jobz,
                     &p, &p, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     la_work, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    for (i = 0; i < p; i++) {
        for (j = 0; j < p; j++)
            Ux[i + j * p] = la_vt[j + i * max_m_p];
        Dx[i * nPlus] = sqrt( la_s[i] );
    }
    for (i = 0; i < p; i++) {
        for (j = 0; j < i; j++) 
            tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
                REAL(VECTOR_ELT(mod,5))[i + p * j]; /* W */
        tmpMat[i + la_m * j] = REAL(VECTOR_ELT(mod,5))[i + p * j];
    }
    F77_CALL(dgesdd)(&la_jobz,
                     &p, &p, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     la_work, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    for (i = 0; i < p; i++) {
        tmp = sqrt( la_s[i] );
        for (j = 0; j < p; j++) 
            sqrtW[i + j * p] = tmp * la_vt[i + j * max_m_p];
    }
    for (i = 0; i < m; i++) {
        for (j = 0; j<i; j++) 
            tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
                REAL(VECTOR_ELT(mod,3))[i + m * j]; /* V */
        tmpMat[i + la_m * j] = REAL(VECTOR_ELT(mod,3))[i + m * j];
    }
    F77_CALL(dgesdd)(&la_jobz,
                     &m, &m, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     la_work, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    for (i = 0; i < m; i++) {
        tmp = sqrt( la_s[i] );
        tmp1 = 1 / tmp;
        tmp1 = R_FINITE(tmp1) ? tmp1 : 0.0;
        for (j = 0; j<m; j++) {
            sqrtV[i + j * m] = tmp * la_vt[i + j * max_m_p];
            sqrtVinv[i + j * m] = tmp1 * la_vt[i + j * max_m_p];
        }
    }
    for (i = 0; i < m; i++) {
        for (j = 0; j < i; j++) {
            tmp = 0.0;
            for (k = 0; k < m; k++)
                tmp += sqrtVinv[k + i * m] * sqrtVinv[k + j * m];
            tmpMat[i + j * la_m] = tmpMat[j + i * la_m] = tmp;
        }
        tmp = 0.0;
        for (k = 0; k < m; k++)
            tmp += SQR( sqrtVinv[k + i * m] );
        tmpMat[i + i * la_m] = tmp;
    }
    for (i = 0; i < p ; i++) 
        for (j = 0; j < m; j++) {
            tmp = 0.0;
            for (k = 0; k < m; k++)
                tmp += sFF[k + i * m] * tmpMat[k + j * la_m];
            tF_Vinv[i + j * p] = tmp;
        }

    /** loop over observations **/
    for (t = 0; t < n; t++) { 
	/** allocate matrices and make pointers point to 'current' **/
	SET_VECTOR_ELT(Ux_priorR, t, allocMatrix(REALSXP, p, p));
	Ux_prior = REAL(VECTOR_ELT(Ux_priorR, t));
	SET_VECTOR_ELT(UxR, t+1, allocMatrix(REALSXP, p, p));
        /** Prior **/
        for (i = 0; i < p; i++) {
            tmp = 0.0;
            for (k = 0; k < p; k++)
                tmp += sGG[i + p * k] * sm0[k * nPlus];
            a[i * n] = tmp;
        }
        
        for (i = 0; i < p; i++) {
            tmp1 = Dx[i * nPlus];
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (l = 0; l < p; l++)
                    tmp += sGG[j + l * p] * Ux[i * p + l];
                tmpMat[i + j * la_m] = tmp * tmp1;
            }
            
        }
        for (i = 0; i < p; i++) 
            for (j = 0; j < p; j++) 
                tmpMat[i + p + j * la_m] = sqrtW[i + j * p];
        l = 2 * p;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &p, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < p; i++) {
            for (j = 0; j < p; j++)
                Ux_prior[i + j * p] = la_vt[j + i * max_m_p];
            Dx_prior[i * n] = la_s[i];
        }

        /** One-step forecast **/
        for (i = 0; i < m; i++) {
            tmp = 0.0;
            for (j = 0; j < p; j++)
                tmp += sFF[i + j * m] * a[j * n];
            f[i * n] = tmp;
        }
        for (i = 0; i < p; i++) {
            tmp1 = Dx_prior[i * n];
            for (j = 0; j < m; j++) {
                tmp = 0.0;
                for (l = 0; l < p; l++)
                    tmp += sFF[j + l * m] * Ux_prior[l + i * p];
                tmpMat[i + j * la_m] = tmp * tmp1;
            }
            
        }
        for (i = 0; i < m; i++) 
            for (j = 0; j < m; j++) 
                tmpMat[i + p + j * la_m] = sqrtV[i + j * m];
        
        l = p + m;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &m, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++)
                Uy[i + j * m] = la_vt[j + i * max_m_p];
            Dy[i] = la_s[i];
        }
        
        /** Posterior **/
        sm0++; Dx++;
	Ux = REAL(VECTOR_ELT(UxR, t+1));
	for (i = 0; i < m; i++)
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (k = 0; k < m; k++)
                    tmp += sqrtVinv[i + k * m] * sFF[k + j * m];
                tmpMat2[i + j * m] = tmp;
            }
        for (i = 0; i < m; i++)
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += tmpMat2[i + k * m] * Ux_prior[k + j * p];
                tmpMat[i + j * la_m] = tmp;
            }
        for (i = 0; i < p; i++) {
            tmp = 1 / Dx_prior[i * n];
            tmpMat[i + m + i * la_m] = R_FINITE(tmp) ? tmp : 0.0;
            for (j = i + 1; j < p; j++) 
                tmpMat[j + m + i * la_m] = tmpMat[i + m + j * la_m] = 0.0;
        }
        l = p + m;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &p, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < p; i++) {
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += Ux_prior[i + k * p] * la_vt[j + k * max_m_p];
                Ux[i + j * p] = tmp;
            }
            tmp = 1 / la_s[i];
            Dx[i * nPlus] = R_FINITE(tmp) ? tmp : 0.0;
        }
        for (i = 0; i < m; i++)
            e[i] = sy[i * n] - f[i * n];
	for (i = 0; i < p; i++)
            for (j = 0; j < m; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += SQR(Dx[i * nPlus]) * Ux[k + i * p] * tF_Vinv[k + j * p];
                tmpMat[i + j * la_m] = tmp;
            }
        for (i = 0; i < p; i++)
            for (j = 0; j < m; j++) {
                tmp = 0.0;
                for (k = 0; k < p; k++)
                    tmp += Ux[i + k * p] * tmpMat[k + j * la_m];
                tmpMat2[j + i * m] = tmp;
            }
        for (i = 0; i < p; i++) {
            tmp = a[i * n];
            for (j = 0; j < m; j++)
                tmp += tmpMat2[j + i * m] * e[j];
            sm0[i * nPlus] = tmp;
        }
/*          pmatrix("Ux",Ux,p,p,p); pmatrix("Dx",Dx, 1,1,p);  */
/*          pmatrix("m0",sm0,1,1,p); error(""); */

	/** increment pointers **/
	sy++; a++; Dx_prior++, f++;
    }

    PROTECT(val = allocVector(VECSXP, 7));
    SET_VECTOR_ELT(val, 0, mR);
    SET_VECTOR_ELT(val, 1, UxR);
    SET_VECTOR_ELT(val, 2, DxR);
    SET_VECTOR_ELT(val, 3, aR);
    SET_VECTOR_ELT(val, 4, Ux_priorR);
    SET_VECTOR_ELT(val, 5, Dx_priorR);
    SET_VECTOR_ELT(val, 6, fR);
    UNPROTECT(8);
    
    return(val);
}


SEXP dlmSmooth0(SEXP mod, SEXP big)
{
/***** Warning: the function relies on the order of the  *****/
/***** components of the list 'mod', not on their names. *****/     

    SEXP val, smoothR, US_R, DS_R;
    int t, i, j, k, p, n, nPlus, la_m, la_n, *la_iwork, la_lwork, la_info=0;
    double *sm, *smooth, *DS, *US, *sGG, *DC, *UC, *DR, *UR, *Ht, *C, *Rinv, 
	*sqrtWinv, tmp, tmp1,
	*tmpMat, *dPointer, *dptr, *dPointer1, *dptr1, *dptr2,
	*la_s, *la_u, *la_vt, *la_work, dBig = REAL(big)[0], eps;
    char la_jobz='S';

    eps = 1.0/dBig;
    p = INTEGER(getAttrib( VECTOR_ELT(mod,7), R_DimSymbol ))[0];
    sGG = REAL(VECTOR_ELT(mod, 6));
    nPlus = length( VECTOR_ELT(mod,1) ); 
    n = nPlus - 1;
    la_m = 2 * p; 
    la_n = p;

    PROTECT(smoothR = allocMatrix(REALSXP, nPlus, p));
    PROTECT(US_R = allocVector(VECSXP, nPlus));
    PROTECT(DS_R = allocMatrix(REALSXP, nPlus, p));
    /* set last value */
    SET_VECTOR_ELT(US_R, n, allocMatrix(REALSXP, p, p));
    US = REAL(VECTOR_ELT(US_R, n));
    DS = REAL(DS_R) + n;
    smooth = REAL(smoothR) + n;
    dptr = REAL(VECTOR_ELT(VECTOR_ELT(mod, 1), n));
    for (i = 0, k = p * p; i < k; i++)
	US[i] = dptr[i];
    dPointer = DS;
    dPointer1 = smooth;
    dptr = REAL(VECTOR_ELT(mod, 2))+ n;
    dptr1 = REAL(VECTOR_ELT(mod, 0)) + n;
    for (i = 0; i < p; i++) {
	*dPointer1 = *dptr1;
	*dPointer = *dptr;
	dPointer += nPlus;
	dPointer1 += nPlus;
	dptr += nPlus;
	dptr1 += nPlus;
    }
    sm = REAL(VECTOR_ELT(mod, 0)) + n - 1; /* m[n-2,] */

    /* allocate space for a la_m by la_n matrix */
    tmpMat = (double *) R_alloc( la_m * la_n, sizeof(double) ); 
    /* space for singular values */
    la_s = (double *) R_alloc( la_n, sizeof(double) );
    /* space for U matrix and Vt matrix (singular vectors) */
    la_u = (double *) R_alloc( la_m * la_n, sizeof(double) );
    la_vt = (double *) R_alloc( la_m * la_n, sizeof(double) );
    /* space for la_iwork */
    la_iwork = (int *) R_alloc( 8 * la_n, sizeof(int) );
    /* space for the square root of W^{-1} */
    sqrtWinv = (double *) R_alloc( p * p, sizeof(double) );
    C = (double *) R_alloc( p * p, sizeof(double) );
    Rinv = (double *) R_alloc( p * p, sizeof(double) );
    Ht = (double *) R_alloc( p * p, sizeof(double) );

    /* ask for optimal size of work array */
    la_lwork = -1;
    F77_CALL(dgesdd)(&la_jobz,
		     &la_m, &la_n, tmpMat, &la_m, la_s,
		     la_u, &la_m,
		     la_vt, &la_n,
		     &tmp, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
	error("error code %d from Lapack routine dgesdd", la_info);
    la_lwork = (int) tmp;
    la_work = (double *) R_alloc( la_lwork, sizeof(double) );
	
    /** preliminaries: compute sqrt(W^(-1)) **/
    dPointer = REAL(VECTOR_ELT(mod,7)); /* W */
    for (i = 0; i < p; i++) {
	for (j = 0; j < i; j++)
	    tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
		dPointer[i + p * j];
	tmpMat[i + la_m * i] = dPointer[i + p * i];
    }
    F77_CALL(dgesdd)(&la_jobz,
		     &p, &p, tmpMat, &la_m, la_s,
		     la_u, &la_m,
		     la_vt, &la_n,
		     la_work, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
	error("error code %d from Lapack routine dgesdd", la_info);
    for (i = 0; i < p; i++) {
	tmp = sqrt(la_s[i]);
 	tmp = ( tmp < eps ? dBig : 1.0 / tmp ); 
	for (j = 0; j < p; j++)
	    sqrtWinv[i + j * p] = tmp * la_vt[i + j * la_n];
    }
    
    /** loop over observations **/
    for (t = n; t > 0; t--) { 
	/* compute H */
	dPointer = tmpMat;
	dptr = UC = REAL(VECTOR_ELT(VECTOR_ELT(mod, 1), t-1)); /* U.C[[t]] */
	dptr1 = DC = REAL(VECTOR_ELT(mod, 2)) + t - 1; /* D.C[t,] */
	for (i = 0; i < p; i++) {
	    tmp = *dptr1;
	    for (j = 0; j < p; j++)
		*dPointer++ = *dptr++ * tmp;
	    dptr1 += nPlus;
	    dPointer += p;
	}
	for (i = 0; i< p; i++) {
	    for (j = 0; j < i; j++) {
		dptr = tmpMat + i;
		dptr1 = tmpMat + j;
		tmp = 0.0;
		for (k = 0; k < p; k++) {
		    tmp += *dptr * *dptr1;
		    dptr += la_m;
		    dptr1 += la_m;
		}
		C[i + j * p] = tmp;
	    }
	    dptr = tmpMat + i;
	    tmp = 0.0;
	    for (k = 0; k < p; k++) {
		tmp += SQR(*dptr);
		dptr += la_m;
	    }
	    C[i + i * p] = tmp;
	} /* lower-tri of C */
	dPointer = tmpMat;
	dptr = REAL(VECTOR_ELT(VECTOR_ELT(mod, 4), t-1)); /* U.R[[t]] */
	dptr1 = REAL(VECTOR_ELT(mod, 5)) + t - 1; /* D.R[t,] */
	for (i = 0; i < p; i++) {
	    tmp = 1.0 / *dptr1;
	    if ( !R_FINITE(tmp) ) tmp = 0.0;
	    for (j = 0; j < p; j++)
		*dPointer++ = *dptr++ * tmp;
	    dptr1 += n;
	    dPointer += p;
	}
	for (i = 0; i< p; i++) {
	    for (j = 0; j < i; j++) {
		dptr = tmpMat + i;
		dptr1 = tmpMat + j;
		tmp = 0.0;
		for (k = 0; k < p; k++) {
		    tmp += *dptr * *dptr1;
		    dptr += la_m;
		    dptr1 += la_m;
		}
		Rinv[i + j * p] = tmp;
	    }
	    dptr = tmpMat + i;
	    tmp = 0.0;
	    for (k = 0; k < p; k++) {
		tmp += SQR(*dptr);
		dptr += la_m;
	    }
	    Rinv[i + i * p] = tmp;
	} /* lower-tri of Rinv */
	for (i = 0; i < p; i++) {
	    dPointer = tmpMat + i;
	    for (j = 0; j < p; j++) {
		dptr = Rinv + i;
		dptr1 = sGG + p * j;
		tmp = 0.0;
		for (k = 0; k < i; k++) {
		    tmp += *dptr * *dptr1++;
		    dptr += p;
		}
		tmp += *dptr++ * *dptr1++;
		for (k = i + 1; k < p; k++) 
		    tmp += *dptr++ * *dptr1++;
		*dPointer = tmp;
		dPointer += la_m;
	    }
	}
	for (i = 0; i < p; i++) {
	    dPointer = Ht + i + p * (p - 1);
	    for (j = p - 1; j >= 0; j--) {
		dptr = tmpMat + i + la_m * (p - 1);
		dptr1 = C + p * (1 + j) - 1;
		tmp = 0.0;
		for (k = p - 1; k > j; k--) { 
		    tmp += *dptr * *dptr1--;
		    dptr -= la_m;
		}
		tmp += *dptr * *dptr1;
		for (k = j - 1; k >= 0; k--) {
		    dptr1 -= p;
		    dptr -= la_m;
		    tmp += *dptr * *dptr1;
		}
		*dPointer = tmp;
		dPointer -= p;
	    }
	} /* t(H) */
	for (j = 0; j < p; j++) {
	    dPointer = tmpMat + j * la_m;
	    for (i = 0; i < p; i++) {
		dptr = sqrtWinv + i;
		dptr1 = sGG + j * p;
		tmp = 0.0;
		for (k = 0; k < p; k++) {
		    tmp += *dptr * *dptr1++;
		    dptr += p;
		}
		*dPointer++ = tmp;
	    }
	}
	for (i = 0; i < p; i++) {
	    dPointer = tmpMat + p + i;
	    tmp = 1.0 / *(DC + i * nPlus);
	    if ( !R_FINITE(tmp) ) tmp = 0.0;
	    dptr = UC + i * p;
	    for (j = 0; j < p; j++) {
		*dPointer = *dptr++ * tmp;
		dPointer += la_m;
	    }
	}
	F77_CALL(dgesdd)(&la_jobz,
			 &la_m, &la_n, tmpMat, &la_m, la_s,
			 la_u, &la_m,
			 la_vt, &la_n,
			 la_work, &la_lwork, la_iwork, &la_info);
	if (la_info != 0)
	    error("error code %d from Lapack routine dgesdd", la_info);
	for (i = 0; i < p; i++) {
	    dPointer = tmpMat + i;
	    tmp = 1.0 / la_s[i];
	    if ( !R_FINITE(tmp) ) tmp = 0.0;
	    dptr = la_vt + i;
	    for (j = 0; j < p; j++) {
		*dPointer = *dptr * tmp;
		dPointer += la_m;
		dptr += la_n;
	    }
	}
	for (i = 0; i < p; i++) {
	    tmp1 = *(DS + i * nPlus);
	    dPointer = tmpMat + p + i;
	    for (j = 0; j < p; j++) {
		dptr = US + i * p;
		dptr1 = Ht + j * p;
		tmp = 0.0;
		for (k = 0; k < p; k++)
		    tmp += *dptr++ * *dptr1++;
		*dPointer = tmp * tmp1;
		dPointer += la_m;
	    }
	}
	F77_CALL(dgesdd)(&la_jobz,
			 &la_m, &la_n, tmpMat, &la_m, la_s,
			 la_u, &la_m,
			 la_vt, &la_n,
			 la_work, &la_lwork, la_iwork, &la_info);
	if (la_info != 0)
	    error("error code %d from Lapack routine dgesdd", la_info);
	/* store results */
	SET_VECTOR_ELT(US_R, t-1, allocMatrix(REALSXP, p, p));
	US = REAL(VECTOR_ELT(US_R, t-1));
  	DS--;
	dPointer1 = DS;
	dptr = la_vt;
	dptr1 = la_s;
	for (i = 0; i < p; i++) {
	    *dPointer1 = *dptr1++;
	    dPointer1 += nPlus;
	    dPointer = US + i;
	    for (j = 0; j < p; j++) {
		*dPointer = *dptr++;
		dPointer += p;
	    }
	}
	dPointer = tmpMat;
	dptr = smooth;
	dptr1 = REAL(VECTOR_ELT(mod, 3)) + t - 1; /* a[t,] */
	for (i = 0; i < p; i++) {
	    *dPointer++ = *dptr - *dptr1;
	    dptr += nPlus;
	    dptr1 += n;
	}
	dPointer = --smooth;
	dptr = Ht;
	dptr2 = sm--;
	for (i = 0; i < p; i++) {
	    dptr1 = tmpMat;
	    *dPointer = 0.0;
	    for (j = 0; j < p; j++) 
		*dPointer += *dptr++ * *dptr1++;
	    *dPointer += *dptr2;
	    dPointer += nPlus;
	    dptr2 += nPlus;
	}
    }
	
    /* set up return value */
    PROTECT(val = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(val, 0, smoothR);
    SET_VECTOR_ELT(val, 1, US_R);
    SET_VECTOR_ELT(val, 2, DS_R);
    UNPROTECT(4);
    return(val);

}

SEXP dlmSmooth(SEXP mod, SEXP tvGG, SEXP tvW, SEXP big)
{
/***** Warning: the function relies on the order of the  *****/
/***** components of the list 'mod', not on their names. *****/     

    SEXP val, smoothR, US_R, DS_R;
    int t, i, j, k, p, n, nPlus, la_m, la_n, *la_iwork, la_lwork, la_info=0,
	stvGG=INTEGER(tvGG)[0], stvW=INTEGER(tvW)[0], *sJGG, *sJW, nrJGG, nrJW,
	nr;
    double *sm, *smooth, *DS, *US, *sGG, *DC, *UC, *DR, *UR, *Ht, *C, *Rinv, 
	*sW, *sqrtWinv, *sX, tmp, tmp1,
	*tmpMat, *dPointer, *dptr, *dPointer1, *dptr1, *dptr2,
	*la_s, *la_u, *la_vt, *la_work, dBig = REAL(big)[0], eps;
    char la_jobz='S';

    eps = 1.0/dBig;
    p = INTEGER(getAttrib( VECTOR_ELT(mod,7), R_DimSymbol ))[0];
    sGG = REAL(VECTOR_ELT(mod, 6));
    sX = REAL(VECTOR_ELT(mod, 10));
    nPlus = length( VECTOR_ELT(mod,1) ); 
    n = nPlus - 1;
    la_m = 2 * p; 
    la_n = p;

    PROTECT(smoothR = allocMatrix(REALSXP, nPlus, p));
    PROTECT(US_R = allocVector(VECSXP, nPlus));
    PROTECT(DS_R = allocMatrix(REALSXP, nPlus, p));
    /* set last value */
    SET_VECTOR_ELT(US_R, n, allocMatrix(REALSXP, p, p));
    US = REAL(VECTOR_ELT(US_R, n));
    DS = REAL(DS_R) + n;
    smooth = REAL(smoothR) + n;
    dptr = REAL(VECTOR_ELT(VECTOR_ELT(mod, 1), n));
    for (i = 0, k = p * p; i < k; i++)
	US[i] = dptr[i];
    dPointer = DS;
    dPointer1 = smooth;
    dptr = REAL(VECTOR_ELT(mod, 2))+ n;
    dptr1 = REAL(VECTOR_ELT(mod, 0)) + n;
    for (i = 0; i < p; i++) {
	*dPointer1 = *dptr1;
	*dPointer = *dptr;
	dPointer += nPlus;
	dPointer1 += nPlus;
	dptr += nPlus;
	dptr1 += nPlus;
    }
    sm = REAL(VECTOR_ELT(mod, 0)) + n - 1; /* m[n-2,] */

    /* allocate space for a la_m by la_n matrix */
    tmpMat = (double *) R_alloc( la_m * la_n, sizeof(double) ); 
    /* space for singular values */
    la_s = (double *) R_alloc( la_n, sizeof(double) );
    /* space for U matrix and Vt matrix (singular vectors) */
    la_u = (double *) R_alloc( la_m * la_n, sizeof(double) );
    la_vt = (double *) R_alloc( la_m * la_n, sizeof(double) );
    /* space for la_iwork */
    la_iwork = (int *) R_alloc( 8 * la_n, sizeof(int) );
    /* space for the square root of W^{-1} */
    sqrtWinv = (double *) R_alloc( p * p, sizeof(double) );
    C = (double *) R_alloc( p * p, sizeof(double) );
    Rinv = (double *) R_alloc( p * p, sizeof(double) );
    Ht = (double *) R_alloc( p * p, sizeof(double) );

    /* ask for optimal size of work array */
    la_lwork = -1;
    F77_CALL(dgesdd)(&la_jobz,
		     &la_m, &la_n, tmpMat, &la_m, la_s,
		     la_u, &la_m,
		     la_vt, &la_n,
		     &tmp, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
	error("error code %d from Lapack routine dgesdd", la_info);
    la_lwork = (int) tmp;
    la_work = (double *) R_alloc( la_lwork, sizeof(double) );
	
    /** preliminaries: compute sqrt(W^(-1)) if time-invariant **/
    if (stvW) {
	sJW = INTEGER(VECTOR_ELT(mod, 9));
	nrJW = INTEGER(getAttrib(VECTOR_ELT(mod, 9), R_DimSymbol))[0];
	sW = REAL(VECTOR_ELT(mod, 7));
    }
    else {
	dPointer = REAL(VECTOR_ELT(mod,7)); /* W */
	for (i = 0; i < p; i++) {
	    for (j = 0; j < i; j++)
		tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
		    dPointer[i + p * j];
	    tmpMat[i + la_m * i] = dPointer[i + p * i];
	}
	F77_CALL(dgesdd)(&la_jobz,
			 &p, &p, tmpMat, &la_m, la_s,
			 la_u, &la_m,
			 la_vt, &la_n,
			 la_work, &la_lwork, la_iwork, &la_info);
	if (la_info != 0)
	    error("error code %d from Lapack routine dgesdd", la_info);
	for (i = 0; i < p; i++) {
	    tmp = sqrt(la_s[i]);
	    tmp = ( tmp < eps ? dBig : 1.0 / tmp ); 
	    for (j = 0; j < p; j++)
		sqrtWinv[i + j * p] = tmp * la_vt[i + j * la_n];
	}
    }
    if (stvGG) {
	sJGG = INTEGER(VECTOR_ELT(mod, 8));
	nrJGG = INTEGER(getAttrib(VECTOR_ELT(mod, 8), R_DimSymbol))[0];
    }

    
    /** loop over observations **/
    for (t = n; t > 0; t--) { 
	/* set time-varying matrices */
	if (stvW) {
	    for (i = 0; i < nrJW; i++)
		sW[ sJW[i] + p * sJW[i + nrJW] ] = sX[ t - 1 + n * sJW[i + 2 * nrJW] ];
	    dPointer = sW;
	    for (i = 0; i < p; i++) {
		for (j = 0; j < i; j++)
		    tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
			dPointer[i + p * j];
		tmpMat[i + la_m * i] = dPointer[i + p * i];
	    }
	    F77_CALL(dgesdd)(&la_jobz,
			     &p, &p, tmpMat, &la_m, la_s,
			     la_u, &la_m,
			     la_vt, &la_n,
			     la_work, &la_lwork, la_iwork, &la_info);
	    if (la_info != 0)
		error("error code %d from Lapack routine dgesdd", la_info);
	    for (i = 0; i < p; i++) {
		tmp = sqrt(la_s[i]);
		tmp = ( tmp < eps ? dBig : 1.0 / tmp ); 
		for (j = 0; j < p; j++)
		    sqrtWinv[i + j * p] = tmp * la_vt[i + j * la_n];
	    }
	}
	if (stvGG) {
	    for (i = 0; i < nrJGG; i++)
		sGG[ sJGG[i] + p * sJGG[i + nrJGG] ] = sX[ t - 1 + n * sJGG[i + 2 * nrJGG] ];
	}
	/* compute H */
	dPointer = tmpMat;
	dptr = UC = REAL(VECTOR_ELT(VECTOR_ELT(mod, 1), t-1)); /* U.C[[t]] */
	dptr1 = DC = REAL(VECTOR_ELT(mod, 2)) + t - 1; /* D.C[t,] */
	for (i = 0; i < p; i++) {
	    tmp = *dptr1;
	    for (j = 0; j < p; j++)
		*dPointer++ = *dptr++ * tmp;
	    dptr1 += nPlus;
	    dPointer += p;
	}
	for (i = 0; i< p; i++) {
	    for (j = 0; j < i; j++) {
		dptr = tmpMat + i;
		dptr1 = tmpMat + j;
		tmp = 0.0;
		for (k = 0; k < p; k++) {
		    tmp += *dptr * *dptr1;
		    dptr += la_m;
		    dptr1 += la_m;
		}
		C[i + j * p] = tmp;
	    }
	    dptr = tmpMat + i;
	    tmp = 0.0;
	    for (k = 0; k < p; k++) {
		tmp += SQR(*dptr);
		dptr += la_m;
	    }
	    C[i + i * p] = tmp;
	} /* lower-tri of C */
	dPointer = tmpMat;
	dptr = REAL(VECTOR_ELT(VECTOR_ELT(mod, 4), t-1)); /* U.R[[t]] */
	dptr1 = REAL(VECTOR_ELT(mod, 5)) + t - 1; /* D.R[t,] */
	for (i = 0; i < p; i++) {
	    tmp = 1.0 / *dptr1;
	    if ( !R_FINITE(tmp) ) tmp = 0.0;
	    for (j = 0; j < p; j++)
		*dPointer++ = *dptr++ * tmp;
	    dptr1 += n;
	    dPointer += p;
	}
	for (i = 0; i< p; i++) {
	    for (j = 0; j < i; j++) {
		dptr = tmpMat + i;
		dptr1 = tmpMat + j;
		tmp = 0.0;
		for (k = 0; k < p; k++) {
		    tmp += *dptr * *dptr1;
		    dptr += la_m;
		    dptr1 += la_m;
		}
		Rinv[i + j * p] = tmp;
	    }
	    dptr = tmpMat + i;
	    tmp = 0.0;
	    for (k = 0; k < p; k++) {
		tmp += SQR(*dptr);
		dptr += la_m;
	    }
	    Rinv[i + i * p] = tmp;
	} /* lower-tri of Rinv */
	for (i = 0; i < p; i++) {
	    dPointer = tmpMat + i;
	    for (j = 0; j < p; j++) {
		dptr = Rinv + i;
		dptr1 = sGG + p * j;
		tmp = 0.0;
		for (k = 0; k < i; k++) {
		    tmp += *dptr * *dptr1++;
		    dptr += p;
		}
		tmp += *dptr++ * *dptr1++;
		for (k = i + 1; k < p; k++) 
		    tmp += *dptr++ * *dptr1++;
		*dPointer = tmp;
		dPointer += la_m;
	    }
	}
	for (i = 0; i < p; i++) {
	    dPointer = Ht + i + p * (p - 1);
	    for (j = p - 1; j >= 0; j--) {
		dptr = tmpMat + i + la_m * (p - 1);
		dptr1 = C + p * (1 + j) - 1;
		tmp = 0.0;
		for (k = p - 1; k > j; k--) { 
		    tmp += *dptr * *dptr1--;
		    dptr -= la_m;
		}
		tmp += *dptr * *dptr1;
		for (k = j - 1; k >= 0; k--) {
		    dptr1 -= p;
		    dptr -= la_m;
		    tmp += *dptr * *dptr1;
		}
		*dPointer = tmp;
		dPointer -= p;
	    }
	} /* t(H) */
	for (j = 0; j < p; j++) {
	    dPointer = tmpMat + j * la_m;
	    for (i = 0; i < p; i++) {
		dptr = sqrtWinv + i;
		dptr1 = sGG + j * p;
		tmp = 0.0;
		for (k = 0; k < p; k++) {
		    tmp += *dptr * *dptr1++;
		    dptr += p;
		}
		*dPointer++ = tmp;
	    }
	}
	for (i = 0; i < p; i++) {
	    dPointer = tmpMat + p + i;
	    tmp = 1.0 / *(DC + i * nPlus);
	    if ( !R_FINITE(tmp) ) tmp = 0.0;
	    dptr = UC + i * p;
	    for (j = 0; j < p; j++) {
		*dPointer = *dptr++ * tmp;
		dPointer += la_m;
	    }
	}
	F77_CALL(dgesdd)(&la_jobz,
			 &la_m, &la_n, tmpMat, &la_m, la_s,
			 la_u, &la_m,
			 la_vt, &la_n,
			 la_work, &la_lwork, la_iwork, &la_info);
	if (la_info != 0)
	    error("error code %d from Lapack routine dgesdd", la_info);
	for (i = 0; i < p; i++) {
	    dPointer = tmpMat + i;
	    tmp = 1.0 / la_s[i];
	    if ( !R_FINITE(tmp) ) tmp = 0.0;
	    dptr = la_vt + i;
	    for (j = 0; j < p; j++) {
		*dPointer = *dptr * tmp;
		dPointer += la_m;
		dptr += la_n;
	    }
	}
	for (i = 0; i < p; i++) {
	    tmp1 = *(DS + i * nPlus);
	    dPointer = tmpMat + p + i;
	    for (j = 0; j < p; j++) {
		dptr = US + i * p;
		dptr1 = Ht + j * p;
		tmp = 0.0;
		for (k = 0; k < p; k++)
		    tmp += *dptr++ * *dptr1++;
		*dPointer = tmp * tmp1;
		dPointer += la_m;
	    }
	}
	F77_CALL(dgesdd)(&la_jobz,
			 &la_m, &la_n, tmpMat, &la_m, la_s,
			 la_u, &la_m,
			 la_vt, &la_n,
			 la_work, &la_lwork, la_iwork, &la_info);
	if (la_info != 0)
	    error("error code %d from Lapack routine dgesdd", la_info);
	/* store results */
	SET_VECTOR_ELT(US_R, t-1, allocMatrix(REALSXP, p, p));
	US = REAL(VECTOR_ELT(US_R, t-1));
  	DS--;
	dPointer1 = DS;
	dptr = la_vt;
	dptr1 = la_s;
	for (i = 0; i < p; i++) {
	    *dPointer1 = *dptr1++;
	    dPointer1 += nPlus;
	    dPointer = US + i;
	    for (j = 0; j < p; j++) {
		*dPointer = *dptr++;
		dPointer += p;
	    }
	}
	dPointer = tmpMat;
	dptr = smooth;
	dptr1 = REAL(VECTOR_ELT(mod, 3)) + t - 1; /* a[t,] */
	for (i = 0; i < p; i++) {
	    *dPointer++ = *dptr - *dptr1;
	    dptr += nPlus;
	    dptr1 += n;
	}
	dPointer = --smooth;
	dptr = Ht;
	dptr2 = sm--;
	for (i = 0; i < p; i++) {
	    dptr1 = tmpMat;
	    *dPointer = 0.0;
	    for (j = 0; j < p; j++) 
		*dPointer += *dptr++ * *dptr1++;
	    *dPointer += *dptr2;
	    dPointer += nPlus;
	    dptr2 += nPlus;
	}
    }
	
    /* set up return value */
    PROTECT(val = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(val, 0, smoothR);
    SET_VECTOR_ELT(val, 1, US_R);
    SET_VECTOR_ELT(val, 2, DS_R);
    UNPROTECT(4);
    return(val);

}

SEXP dlmForecast(SEXP mod, SEXP nAhead)
{
/***** Warning: the function relies on the order of the  *****/
/***** components of the list 'mod', not on their names. *****/     

    SEXP val, UxR, DxR, aR, fR, UyR, DyR;
    int i, j, k, l, p, m, n, nPlus, t, max_m_p, 
	la_m, la_n, la_info=0, la_lwork, *la_iwork;
    double *sFF, *sGG, *Ux, *Ux_next,  
        *Dx, *sqrtV, *sqrtW, *a, *f, *Uy, *Dy;
    double tmp, tmp1, *tmpMat, *la_s, *la_u, *la_vt, *la_work;
    char la_jobz='S';

    m = INTEGER(getAttrib( VECTOR_ELT(mod,2), R_DimSymbol ))[0];
    p = INTEGER(getAttrib( VECTOR_ELT(mod,2), R_DimSymbol ))[1];
    max_m_p = m > p ? m : p;
    la_n = max_m_p;
    la_m = 2 * la_n;
    n = INTEGER(nAhead)[0]; 
    nPlus = n + 1;
    PROTECT(aR = allocMatrix(REALSXP, nPlus, p));
    PROTECT(UxR = allocVector(VECSXP, nPlus));
    PROTECT(DxR = allocMatrix(REALSXP, nPlus, p));
    PROTECT(fR = allocMatrix(REALSXP, n, m));
    PROTECT(UyR = allocVector(VECSXP, n));
    PROTECT(DyR = allocMatrix(REALSXP, n, m));
    for (i = 0; i < p; i++)
	REAL(aR)[i * nPlus] = REAL(VECTOR_ELT(mod,0))[i];
    a = REAL(aR);
    Dx = REAL(DxR);
    SET_VECTOR_ELT(UxR, 0, allocMatrix(REALSXP, p, p));
    Ux = REAL(VECTOR_ELT(UxR, 0));
    f = REAL(fR);
    Dy = REAL(DyR);
    SET_VECTOR_ELT(UyR, 0, allocMatrix(REALSXP, m, m));
    Uy = REAL(VECTOR_ELT(UyR, 0));
    sFF = REAL(VECTOR_ELT(mod,2));
    sGG = REAL(VECTOR_ELT(mod,4));
    sqrtV = (double *) R_alloc( m * m, sizeof(double) );
    sqrtW = (double *) R_alloc( p * p, sizeof(double) );

    /* allocate space for a la_m by la_n matrix */
    tmpMat = (double *) R_alloc( la_m * la_n, sizeof(double) ); 
    /* space for singular values */
    la_s = (double *) R_alloc( max_m_p, sizeof(double) );
    /* space for U matrix and Vt matrix (singular vectors) */
    la_u = (double *) R_alloc( la_m * max_m_p, sizeof(double) );
    la_vt = (double *) R_alloc( max_m_p * max_m_p, sizeof(double) );
    /* space for la_iwork */
    la_iwork = (int *) R_alloc( 8 * p, sizeof(int) );

    /* ask for optimal size of work array */
    la_lwork = -1;
    F77_CALL(dgesdd)(&la_jobz,
                     &la_m, &la_n, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     &tmp, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    la_lwork = (int) tmp;
    la_work = (double *) R_alloc( la_lwork, sizeof(double) );

    /** preliminaries: compute svd of C0, sqrt(V), sqrt(W), etc... **/
    for (i = 0; i < p; i++) {
        for (j = 0; j<i; j++) 
            tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
                REAL(VECTOR_ELT(mod,1))[i + p * j]; /* C0 */
        tmpMat[i + la_m * j] = REAL(VECTOR_ELT(mod,1))[i + p * j];
    }
    F77_CALL(dgesdd)(&la_jobz,
                     &p, &p, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     la_work, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    for (i = 0; i < p; i++) {
        for (j = 0; j < p; j++)
            Ux[i + j * p] = la_vt[j + i * max_m_p];
        Dx[i * nPlus] = sqrt( la_s[i] );
    }
    for (i = 0; i < p; i++) {
        for (j = 0; j < i; j++) 
            tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
                REAL(VECTOR_ELT(mod,5))[i + p * j]; /* W */
        tmpMat[i + la_m * j] = REAL(VECTOR_ELT(mod,5))[i + p * j];
    }
    F77_CALL(dgesdd)(&la_jobz,
                     &p, &p, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     la_work, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    for (i = 0; i < p; i++) {
        tmp = sqrt( la_s[i] );
        for (j = 0; j < p; j++) 
            sqrtW[i + j * p] = tmp * la_vt[i + j * max_m_p];
    }
    for (i = 0; i < m; i++) {
        for (j = 0; j<i; j++) 
            tmpMat[i + la_m * j] = tmpMat[j + la_m * i] = 
                REAL(VECTOR_ELT(mod,3))[i + m * j]; /* V */
        tmpMat[i + la_m * j] = REAL(VECTOR_ELT(mod,3))[i + m * j];
    }
    F77_CALL(dgesdd)(&la_jobz,
                     &m, &m, tmpMat, &la_m, la_s,
                     la_u, &la_m,
                     la_vt, &max_m_p,
                     la_work, &la_lwork, la_iwork, &la_info);
    if (la_info != 0)
        error("error code %d from Lapack routine dgesdd", la_info);
    for (i = 0; i < m; i++) {
        tmp = sqrt( la_s[i] );
        for (j = 0; j<m; j++) 
            sqrtV[i + j * m] = tmp * la_vt[i + j * max_m_p];
    }
    
    /** loop over future times **/
    for (t = 0; t < n; t++) { 
	/** allocate matrices and make pointers point to 'current' **/
	SET_VECTOR_ELT(UxR, t+1, allocMatrix(REALSXP, p, p));
	Ux_next = REAL(VECTOR_ELT(UxR, t+1));
	SET_VECTOR_ELT(UyR, t, allocMatrix(REALSXP, m, m));
	Uy = REAL(VECTOR_ELT(UyR, t));
        /** Prior **/
        for (i = 0; i < p; i++) {
            tmp = 0.0;
            for (k = 0; k < p; k++)
                tmp += sGG[i + p * k] * a[k * nPlus];
            a[i * nPlus + 1] = tmp;
        }
        
        for (i = 0; i < p; i++) {
            tmp1 = Dx[i * nPlus];
            for (j = 0; j < p; j++) {
                tmp = 0.0;
                for (l = 0; l < p; l++)
                    tmp += sGG[j + l * p] * Ux[i * p + l];
                tmpMat[i + j * la_m] = tmp * tmp1;
            }
            
        }
        for (i = 0; i < p; i++) 
            for (j = 0; j < p; j++) 
                tmpMat[i + p + j * la_m] = sqrtW[i + j * p];
        l = 2 * p;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &p, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < p; i++) {
            for (j = 0; j < p; j++)
                Ux_next[i + j * p] = la_vt[j + i * max_m_p];
            Dx[i * nPlus + 1] = la_s[i];
        }

	/** increment pointers (states) **/
	a++; Dx++;
	Ux = Ux_next;

        /** One-step forecast **/
        for (i = 0; i < m; i++) {
            tmp = 0.0;
            for (j = 0; j < p; j++)
                tmp += sFF[i + j * m] * a[j * nPlus];
            f[i * n] = tmp;
        }
        for (i = 0; i < p; i++) {
            tmp1 = Dx[i * nPlus];
            for (j = 0; j < m; j++) {
                tmp = 0.0;
                for (l = 0; l < p; l++)
                    tmp += sFF[j + l * m] * Ux[l + i * p];
                tmpMat[i + j * la_m] = tmp * tmp1;
            }
            
        }
        for (i = 0; i < m; i++) 
            for (j = 0; j < m; j++) 
                tmpMat[i + p + j * la_m] = sqrtV[i + j * m];
        
        l = p + m;
        F77_CALL(dgesdd)(&la_jobz,
                         &l, &m, tmpMat, &la_m, la_s,
                         la_u, &la_m,
                         la_vt, &max_m_p,
                         la_work, &la_lwork, la_iwork, &la_info);
        if (la_info != 0)
            error("error code %d from Lapack routine dgesdd", la_info);
        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++)
                Uy[i + j * m] = la_vt[j + i * max_m_p];
            Dy[i * n] = la_s[i];
        }

	/** increment pointers (observables) **/
	f++; Dy++;
    }

    PROTECT(val = allocVector(VECSXP, 6));
    SET_VECTOR_ELT(val, 0, aR);
    SET_VECTOR_ELT(val, 1, UxR);
    SET_VECTOR_ELT(val, 2, DxR);
    SET_VECTOR_ELT(val, 3, fR);
    SET_VECTOR_ELT(val, 4, UyR);
    SET_VECTOR_ELT(val, 5, DyR);
    UNPROTECT(7);
    
    return(val);
}


SEXP ARtranspar(SEXP Rp, SEXP Rraw)
{
/*** Adapted from 'partrans' in arima.c ***/ 
    int j, k, p;
    double a, work[50], *new, *raw=REAL(Rraw);
    SEXP Rnew;

    p = INTEGER(Rp)[0];
    if(p > 50) error("can only transform 50 pars in ARtranspar");
    PROTECT( Rnew = allocVector(REALSXP, p));
    new = REAL(Rnew);

    /* Step one: map (-Inf, Inf) to (-1, 1) via tanh
       The parameters are now the pacf phi_{kk} */
    for(j = 0; j < p; j++) work[j] = new[j] = tanh(raw[j]);
    /* Step two: run the Durbin-Levinson recursions to find phi_{j.},
       j = 2, ..., p and phi_{p.} are the autoregression coefficients */
    for(j = 1; j < p; j++) {
        a = new[j];
        for(k = 0; k < j; k++)
            work[k] -= a * new[j - k - 1];
        for(k = 0; k < j; k++) new[k] = work[k];
    }
    UNPROTECT(1);
    return( Rnew );
}

