/*
Authors: Sarah Al-Assam, Stephen Clark and Dieter Jaksch
$LastChangedDate: 2013-10-09 10:02:48 +0100 (Wed, 09 Oct 2013) $
(c) University of Oxford 2014
*/

/*! \file tntMpsExpecOutput.c
 *  \brief This contains the routine for calculating, displaying then saving the expectation value of single site operators for a 1D MPS with open boundary conditions
 *
 */

/* Include the header for the TNT MPS library */
#include "tntMpsInternal.h"

/* Internal function to avoid repeated code */
void _tnt_MpsExpecOutputOneOp(tntComplexArray *Exval, int printOutput, int saveOutput, char *savename, char *label, unsigned counter);


/*!
 * \ingroup mps
 * 
 * Calculates various types of one-site and two-site expectation values, the operators for which are passed in the ::tntMpsExOp structure.
 *
 * The function can output the expectation values to the screen as well as saving them to the named output file. 
 * Use the appropriate flags to specify whether results should be output to the screen or to an output file. 
 * The function will have no effect (and nothing will be calculated) if both the screen output flag is \c 0 and the printOutput flag is \c 0 - 
 * in this case a warning will be printed to screen but execution will continue.
 * 
 * Note the expecation value of operator \f$ \hat{o} \f$ is given by \f$ \langle\psi|\hat{o}|\psi\rangle/\langle\psi|\psi\rangle \f$ i.e. it uses the normalised value of the MPS wave function.
 * However the function does not renormalise the MPS wave function to do this, it simply divides by the norm-squared \f$ \langle\psi|\psi\rangle \f$.
 * The value of the norm squared will be printed to the screen and saved to the output file if the appropriate flags are given.
 *
 * \return No return value - the expectation values are output to the screen and/or to an output file.
 *
 * \see tntMpsProcessExpecOptions() for a function which converts command line arguments to the arrays of operators that are required by this function.
 */
void tntMpsExpecOutput(tntNetwork mps,      /*!< The network representing the MPS. Unchanged by the function. */
                       tntMpsExOp *Op,      /*!< The operators for calculating expectation values */
                       int orthCentre,      /*!< Orthogonality centre of the MPS. */
                       unsigned printOutput,/*!< Set to 1 to print expectation values to screen, 0 otherwise */
                       unsigned saveOutput, /*!< Set to 1 to save values to output file, 0 otherwise */
                       char *savename,      /*!< Path to the output file. If the extension is missing ".mat" is added. Only used if saveOutput is 1. */
                       unsigned counter)    /*!< can be used if saving multiple batches of expectation values e.g. for different timesteps. 
                                                 It appends the values to a pre-exisiting array in the position given by counter.
                                                 Pass zero if a counter is not required. */
{

    tntComplexArray Exval;  /* array for holding expectation values */
    double normsq;          /* Norm squared of the wave function. Always calculated */
    unsigned loop;          /* Used to loop through the different expectation value types */
    tntNodeArray op;        /* Node array for holding operators for expectation value currently being calculated */
    tntIntArray sitenum;    /* Array for holding site number for expectation value currently being calculated */
    int L;                  /* Length of the MPS */    
    tntNodeArray betap, gammap, *betas = &betap, *gammas = &gammap;    /* Node arrays for precontracted parts of the network */
    
    if (!printOutput && !saveOutput) {
        tntWarningPrint("No expectation values are being calculated as both output flags are set to off"); /* NO_COVERAGE */
        return; /* NO_COVERAGE */
    }  /* NO_COVERAGE */
    
    tntSysQNClearWarnOff();
    
    /* Find the length of the MPS */
    L = (int) tntMpsLength(mps);
    
    /* Calculate the norm squared */
    normsq = tntMpsSelfProduct(mps,orthCentre);
    
    if (printOutput) {
        tntPrintf("\n~~~~~~~             Norm squared is % #6.6g             ~~~~~~~~\n\n",normsq);
    }
    if (saveOutput) {
        if (counter) {
            tntDoubleParamsUpdate(savename, counter, normsq);
        } else {
            tntDoubleParamsSave(savename, normsq);
        }
        
    }
    
    /* Calculate precontracted nodes since multiple types of the same contraction will be performed */
    if (0 == orthCentre) {
        tntMpsSelfInit(mps, betas, NULL);
        gammas = NULL;
    } else if (L-1 == orthCentre) {
        tntMpsSelfInit(mps, NULL, gammas);
        betas = NULL;
    } else {
        tntMpsSelfInit(mps, betas, gammas);
    }
    
    
    /* ----- Onsite expectation values ----- */
    
    if (Op->ExOpOs.sz > 0) {
    
        /* Allocate arrays of size one for onsite operators */
        op = tntNodeArrayAlloc(1);
        sitenum = tntIntArrayAlloc(1);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~          Onsite expectation values          ~~~~~~~~\n");
        }
    
    }
    /* Find the single site expectation values of the calculated ground state */
    for (loop = 0; loop < Op->ExOpOs.sz; loop++) {
        
        /* Allocate memory for storing single-site expectation values */
        Exval = tntComplexArrayAlloc(L);
        
        /* Copy over the required operator */
        op.vals[0] = tntNodeCopy(Op->ExOpOs.vals[loop],0);
        
        for (sitenum.vals[0] = 0; sitenum.vals[0] < L; (sitenum.vals[0])++) {
            /* Find the expectation value on each site using a product MPO. */
            Exval.vals[sitenum.vals[0]] = tntMpsPmpoMpsProduct(mps, &op, &sitenum, orthCentre, betas, gammas);
            Exval.vals[sitenum.vals[0]].re /= normsq;
            Exval.vals[sitenum.vals[0]].im /= normsq;
        }
        
        /* output and free values */
        _tnt_MpsExpecOutputOneOp(&Exval, printOutput, saveOutput, savename, Op->LbOpOs.vals[loop], counter);
        
        /* Free the copy of the operator */
        tntNodeFree(&(op.vals[0]));
    }
    if (Op->ExOpOs.sz > 0) {
        
        /* free the arrays */
        tntIntArrayFree(&sitenum);
        tntNodeArrayFree(&op);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        }
    }
    
    
    /* ---------- Nearest-neighbour expectation values ----------- */
    
    if (Op->ExOp2nn.sz > 0) {
        
        /* Allocate arrays of size two for two-site operator site numbers */
        sitenum = tntIntArrayAlloc(2);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~    Nearest-neighbour expectation values     ~~~~~~~~\n");
        }
        
    }
    for (loop = 0; loop < Op->ExOp2nn.sz; loop+=2) {
        
        /* Allocate memory for storing two-site expectation values */
        Exval = tntComplexArrayAlloc(L-1);
        
        /* Allocate arrays of size two for two-site operators */
        op = tntNodeArrayAlloc(2);
        
        /* Copy over the required operator */
        op.vals[0] = tntNodeCopy(Op->ExOp2nn.vals[loop],0);
        op.vals[1] = tntNodeCopy(Op->ExOp2nn.vals[loop+1],0);
        
        for (sitenum.vals[0] = 0; sitenum.vals[0] < L-1; (sitenum.vals[0])++) {
            /* assign value for second site */
            sitenum.vals[1] = sitenum.vals[0]+1;
            /* Find the expectation value on each site using a product MPO. */
            Exval.vals[sitenum.vals[0]] = tntMpsPmpoMpsProduct(mps, &op, &sitenum, orthCentre, betas, gammas);
            Exval.vals[sitenum.vals[0]].re /= normsq;
            Exval.vals[sitenum.vals[0]].im /= normsq;
        }
        
        /* Free the operators */
        tntNodeArrayFree(&op);
        
        /* output and free values */
        _tnt_MpsExpecOutputOneOp(&Exval, printOutput, saveOutput, savename, Op->LbOp2nn.vals[loop/2], counter);
        
    }
    if (Op->ExOp2nn.sz > 0) {
        
        /* free the arrays */
        tntIntArrayFree(&sitenum);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        }
    }
    
    /* ---------- Centre-site expectation values ----------- */
    
    if (Op->ExOp2cs.sz > 0) {
        
        /* Allocate arrays of size two for two-site operator site numbers */
        sitenum = tntIntArrayAlloc(2);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~       Centre-site expectation values        ~~~~~~~~\n");
        }
    }
    for (loop = 0; loop < Op->ExOp2cs.sz; loop+=2) {
        
        /* Allocate memory for storing two-site expectation values */
        Exval = tntComplexArrayAlloc(L);

        /* Allocate first site that expectation value will be found on */
        sitenum.vals[0] = L/2;
        
        /* Allocate arrays of size two for two-site operators */
        op = tntNodeArrayAlloc(2);
        
        /* Copy over the required operator */
        op.vals[0] = tntNodeCopy(Op->ExOp2cs.vals[loop]);
        op.vals[1] = tntNodeCopy(Op->ExOp2cs.vals[loop+1]);
        
        for (sitenum.vals[1] = 0; sitenum.vals[1] < L; (sitenum.vals[1])++) {
            /* Find the expectation value on each site using a product MPO. */
            Exval.vals[sitenum.vals[1]] = tntMpsPmpoMpsProduct(mps, &op, &sitenum, orthCentre, betas, gammas);
            Exval.vals[sitenum.vals[1]].re /= normsq;
            Exval.vals[sitenum.vals[1]].im /= normsq;
        }
        
        /* Free the operators */
        tntNodeArrayFree(&op);
        
        /* output and free values */
        _tnt_MpsExpecOutputOneOp(&Exval, printOutput, saveOutput, savename, Op->LbOp2cs.vals[loop/2], counter);
        
    }
    if (Op->ExOp2cs.sz > 0) {
        
        /* free the arrays */
        tntIntArrayFree(&sitenum);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        }
    }
    
    
    /* ---------- All pairs of two-site expectation values ----------- */
    
    if (Op->ExOp2ap.sz > 0) {
        
        /* Allocate arrays of size two for sitenumbers for two-site operators */
        sitenum = tntIntArrayAlloc(2);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~~~ All pairs of two-site expectation values ~~~~~~~~~\n");
        }
        
    }
    for (loop = 0; loop < Op->ExOp2ap.sz; loop+=2) {
        
        int ind; /* total matrix index */
        
        /* Allocate memory for storing two-site expectation values */
        Exval = tntComplexArrayAlloc(L,L);
        
        /* Allocate arrays of size two for two-site operators */
        op = tntNodeArrayAlloc(2);
        
        /* Copy over the required operator */
        op.vals[0] = tntNodeCopy(Op->ExOp2ap.vals[loop]);
        op.vals[1] = tntNodeCopy(Op->ExOp2ap.vals[loop+1]);
        
        /* Find the expectation value on each pair of sites using a product MPO. */
        for (sitenum.vals[0] = 0; sitenum.vals[0] < L; (sitenum.vals[0])++) {
            for (sitenum.vals[1] = 0; sitenum.vals[1] < L; (sitenum.vals[1])++) {
                ind = sitenum.vals[0] + L*sitenum.vals[1];
                Exval.vals[ind] = tntMpsPmpoMpsProduct(mps, &op, &sitenum, orthCentre, betas, gammas);
                Exval.vals[ind].re /= normsq;
                Exval.vals[ind].im /= normsq;
            }
        }
        
        /* Free the operators */
        tntNodeArrayFree(&op);
        
        /* output and free values */
        _tnt_MpsExpecOutputOneOp(&Exval, printOutput, saveOutput, savename, Op->LbOp2ap.vals[loop/2], counter);

    }
    if (Op->ExOp2ap.sz > 0) {
        
        /* free the arrays */
        tntIntArrayFree(&sitenum);
        
        if (printOutput) {
            tntPrintf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        }
        
    }
    
    if (0 == orthCentre) {
        tntNodeArrayFree(betas);
    } else if (L-1 == orthCentre) {
        tntNodeArrayFree(gammas);
    } else {
        tntNodeArrayFree(betas);
        tntNodeArrayFree(gammas);
    }

    tntSysQNClearWarnOn();
}

/* prints, saves and frees the complex array, converting to double if necessary, and taking account of flags and counters */
void _tnt_MpsExpecOutputOneOp(tntComplexArray *Exval, int printOutput, int saveOutput, char *savename, char *label, unsigned counter) {
    
    tntDoubleArray Exvalr; /* array for holding real part of expectation values */
    tntComplexArray Exvalc; /* array for holding real part of expectation values */
    
    /* Now take different action depending on whether the expectation values are real or complex */
    if (tntComplexArrayIsReal(Exval)) {
        /* change complex array to real array */
        Exvalr = tntComplexArrayToReal(Exval);
        /* Print real output */
        if (printOutput) {
            tntNamedDoubleArrayPrint(label,&Exvalr);
        }
        /* Save this expectation value */
        if (saveOutput) {
            if (counter) tntDoubleArraysNamedUpdate(savename, counter,Exvalr, label);
            else         tntDoubleArraysNamedSave(savename, Exvalr, label);

        }
        
        /* free real array */
        tntDoubleArrayFree(&Exvalr);
    } else {
        Exvalc = *Exval;
        /* Print complex output */
        if (printOutput) {
            tntNamedComplexArrayPrint(label,Exval);
        }
        /* Save complex output */
        if (saveOutput) {
            if (counter) tntComplexArraysNamedUpdate(savename, counter, Exvalc, label);
            else         tntComplexArraysNamedSave(savename, Exvalc, label);
        }
        /* Free complex array */
        tntComplexArrayFree(Exval);
    }
}

