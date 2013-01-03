#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#define STR(SE) CHAR(STRING_ELT(SE,0))

SEXP gcrma_getSeq(SEXP);


SEXP gcrma_getSeq(SEXP psequence) {
    SEXP outMatrix;
    const char *pseq;
    int i;

    pseq = STR(psequence);
    
    PROTECT(outMatrix = allocMatrix(INTSXP, 4, strlen(pseq)));
    
    for (i = 0; i < strlen(pseq); i++) {
	if (pseq[i] == 'A')
	    INTEGER(outMatrix)[(i * 4)] = 1;
	else
	    INTEGER(outMatrix)[(i * 4)] = 0;
	
	if (pseq[i] == 'C')
	    INTEGER(outMatrix)[(i * 4)+1] = 1;
	else
	    INTEGER(outMatrix)[(i * 4)+1] = 0;

	if (pseq[i] == 'G')
	    INTEGER(outMatrix)[(i * 4)+2] = 1;
	else
	    INTEGER(outMatrix)[(i * 4)+2] = 0;

	if (pseq[i] == 'T')
	    INTEGER(outMatrix)[(i * 4)+3] = 1;
	else
	    INTEGER(outMatrix)[(i * 4)+3] = 0;
    }

    UNPROTECT(1);
    return(outMatrix);
}
