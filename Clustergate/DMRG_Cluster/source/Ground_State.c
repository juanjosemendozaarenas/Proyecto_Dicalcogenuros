/*
Authors: Fernando G—mez-Ruiz (1,2), Adolfo del Campo (2).
        (1) Universidad de Los Andes.
        (2) University of Massachusetts Boston

Date: Febrary 2018

(c) University of Massachusetts Boston

This code calculates the ground state of a 1D lattice, using a matlab initialisation file

*/

#include "tntMps.h"

int main(int argc, char **argv)
{

    char initfile[TNT_STRLEN];          /* Path to initialisation data */
    char saveprefix[TNT_STRLEN];        /* Path to output file */
    char savefile[TNT_STRLEN];          /* Name of output file */

    tntNetwork wf, H, prop;             /* The MPS wave function, the MPO Hamiltonian and the staircase propagator */
    int chi, L, rand_wf;                /* Maximum internal dimension for the MPS, length of the system, and start from rand wf */
    tntNodeArray nnl, nnr, os;          /* operators for building the hamiltonian */
    tntComplexArray nnparam, osparam;   /* parameters for building the hamiltonian */
    double prec, E, Eprev;              /* Precision for the ground state, energy  */
    tntIntArray qn_tot;                 /* Quantum numbers for the ground state */
    int i, i_max;                       /* Iteration counter, maximum number of iterations */
    tntNodeArray HeffL, HeffR;          /* Precontracted nodes for the right and left sides of the network */
    unsigned use_symm;               /* Set to 1 to use symmetries, 0 otherwise */
    int chi_ini_rand, chi_max, delta_chi; /* Initial chi for random initial state (no symmetries), maximal chi, and change of chi when error is too large */
    double err_max;                     /* Maximal error of DMRG sweep, beyond which chi is increased */
    tntMpsExOp ExOp;                    /* Defines all the operators for calculating expectation values */
    double err = 0.0, err_previous, err_step = 0.0;      /* The accumulated truncation error, previous accumulated error, and error of current DMRG step */
    int intermediate;                   /* Set to 1 to continue DMRG simulation when previous one was not complete */

    /* Initialize the TNT library */
    tntInitialize();

    /* Process the command line options */
    tntProcessCLOptions(argc, argv, NULL, initfile, saveprefix, NULL, NULL);

    /* Load all the system parameters */
    tntNodeArraysNamedLoad(initfile,nnl,"nnlg",nnr,"nnrg",os,"osg");
    tntComplexArraysNamedLoad(initfile,nnparam,"nnparamg",osparam,"osparamg");
    tntIntParamsLoad(initfile, chi, L, i_max, rand_wf);
    tntIntParamsLoad(initfile, chi_ini_rand, chi_max, delta_chi);
    tntIntParamsLoad(initfile, use_symm, intermediate);
    tntDoubleParamsLoad(initfile, prec, err_max);
    tntStringsLoad(initfile, savefile);

    /* Read operators for expectation values */
    tntExOpsLoad(initfile, ExOp);

    /* Complete name of saving file */
    strcat(saveprefix, savefile);

    /* Start with a random state with the quantum numbers given if flag given, otherwise load wf, and orthonormalise */
    if ((rand_wf==1) && (use_symm==1)) {
        printf("Creating random initial state with symmetries.\n");
        tntIntArraysLoad(initfile, qn_tot);
        wf = tntMpsCreateSymmRandom(L, qn_tot.vals);
        tntMpsOrthNorm(wf,0);
    } else if ((rand_wf==1) && (use_symm==0)){
        printf("Creating random initial state with no symmetries.\n");
        wf = tntMpsCreateRandom(L,chi_ini_rand);
        tntMpsOrthNorm(wf,0);
        tntNodePrintAll(tntNodeFindFirst(wf));
    } else if (intermediate==1){
        printf("Loading starting state from output file of previous simulation with same parameters.\n");
        tntNetworksNamedLoad(saveprefix, wf, "wf_g");
        tntMpsOrthNorm(wf,0);
    } else {
        printf("Loading starting state from initialisation file.\n");
        tntNetworksLoad(initfile, wf);
        tntMpsOrthNorm(wf,0);
    }
    /* Create the MPO from the system parameters */
    H = tntMpsCreateMpo(L,&nnl,&nnr,&nnparam,&os,&osparam);
    tntNodeArrayFree(&nnl);
    tntNodeArrayFree(&nnr);
    tntNodeArrayFree(&os);
    tntComplexArrayFree(&nnparam);
    tntComplexArrayFree(&osparam);

    /* Initialise the nodes that will be needed for the DMRG sweep */
    tntMpsVarMinMpoInit(wf, H, &HeffL, &HeffR);

    /* Determine the energy of the start state (the norm of the state is assumed to be 1) */
    E = tntComplexToReal(tntMpsMpoMpsProduct(wf, H));
    tntDoubleParamsSave(saveprefix,E);
    tntDoubleParamsSave(saveprefix,err_step);

    printf("--------------------------------------------------------------------------------\n");
    printf("                Starting minimisation sweeps to find ground state               \n");
    printf("--------------------------------------------------------------------------------\n");
    printf("Starting energy is %g. \n", E);

    err_previous = 0;

    /* Perform minimization sweeps */
    for (i = 0; i < i_max; i++) {

        Eprev = E;
        E = tntMpsVarMinMpo2sStep(wf, chi, H, &HeffL, &HeffR, &err);

        err_step = err - err_previous; /* Define error of the current DMRG step */
        err_previous = err; /* Update */

        /* Update energy in output file */
        tntDoubleParamsUpdate(saveprefix,i+1,E);
        tntDoubleParamsUpdate(saveprefix,i+1,err_step);
        printf("Energy is %g. Difference is %g\n", E, Eprev-E);
        printf("For chi = %d, accumulated DMRG error is %g, and error of this step is %g.\n", chi, err, err_step);

        /* Increase truncation parameter if error is too large and chi_max has not been reached */
        if((err_step>err_max) && (chi<chi_max)) {
            chi += delta_chi;
            printf("Increasing truncation parameter to chi = %d.\n", chi);
        }

        /* Check if convergence in energy has been achieved */
        //if ((Eprev - E < prec) && (Eprev - E > 0.0)) break;
        if ((fabs(Eprev - E)) < prec) break;

        /* Save network at every intermediate step, for safety */
        tntNetworksNamedSave(saveprefix, wf, "wf_g");
    }

    /* Save network at the end of the DMRG sweeps */
    tntNetworksNamedSave(saveprefix, wf, "wf_g");

    tntNodeArrayFree(&HeffL);
    tntNodeArrayFree(&HeffR);

    /* Calculate expectation values */
    tntMpsExpecOutput(wf, &ExOp, 0, 1, 1, saveprefix, 0);

    /*  Free all the dynamically allocated nodes and associated dynamically allocated arrays. */
    tntNetworkFree(&wf);
    tntMpsExOpFree(&ExOp);

    /* Finish with the TNT library */
    tntFinalize();

    return 0;
}
