/*
Authors: Sarah Al-Assam, Juan Jose Mendoza-Arenas, Stephen Clark and Dieter Jaksch
Date:    December 2015
(c) University of Oxford 2015
*/

/*! This header file contains prototypes for all the network functions in
 *  the djscripts directory. In particular they contain the prototypes for
 *  functions written by Juan Jose Mendoza-Arenas for time evolution of a
 *  mixed state using TEBD.
 */

/* List of function prototypes */
/* Listed in alphabetical order */

#define tntMpsExOp tntExOp

tntComplex network_expec_single_mixed_obc1d(tntNetwork wf,  /* The current network */
                                        tntNodeArray *op, /* the operator to find the expectation value of */
                                        int site, /* Site at which the reduced density matrix is calculated */
                                        tntNode discard_legs, /* tnode helping to discard legs */
                                        tntNode basis); /* tnode for the local basis */

tntComplex network_expec_two_mixed_obc1d(tntNetwork wf,  /* The current network */
                                     tntNodeArray *op, /* the array of operators to find the expectation value of */
                                     int site1, /* First site */
                                     int site2, /* Second site */
                                     tntNode discard_legs, /* tnode helping to discard legs */
                                     tntNode basis); /* tnode for the local basis */

tntComplex network_expec_three_mixed_obc1d(tntNetwork wf,  /* The current network */
                                       tntNodeArray *op, /* the array of operators to find the expectation value of */
                                       int site, /* First site */
                                       tntNode discard_legs, /* node helping to discard legs */
                                       tntNode basis); /* node for the local basis */

tntNode network_rhoone_mixed_obc1d(tntNetwork wf,  /* The current network */
        int site, /* Site at which the reduced density matrix is calculated */
        tntNode discard_legs, /* tnode helping to discard legs */
        tntNode basis); /* tnode for the local basis */

tntNode network_rhotwo_mixed_obc1d(tntNetwork wf,  /* The current network */
        int site1, /* First site of the reduced density matrix */
        int site2, /* Second site of the reduced density matrix */
        struct tnode *discard_legs, /* tnode helping to discard legs */
        struct tnode *basis); /* tnode for the local basis */

tntNode network_rhothree_mixed_obc1d(tntNetwork wf,  /* The current network */
                                     int site, /* First site of the reduced density matrix */
                                     tntNode discard_legs, /* node helping to discard legs */
                                     tntNode basis); /* node for the local basis */

tntComplex network_trace_mixed_obc1d(tntNetwork wf,  /* The current network */
                                  tntNode discard_legs); /* tnode helping to discard legs */

tntNodeArray tntMpsCreatePropArray_nn_ready(unsigned L, /* Length of system. */
                                   tntComplex h, /* Uniform scale factor to apply to all terms. See the main description for information on how real and imaginary parts are applied */
                                   tntNodeArray *nn, /* Array of nearest-neighbour operators. Send NULL if there are no nearest neighbour terms. Unchanged by function - copies are used. */
                                   tntComplexArray *nnparam, /* Array of parameters for nearest-neighbour operators. Send NULL if there are no nearest neighbour terms. */
                                   tntNodeArray *os,  /* Array of onsite operators. Send NULL if there are no on-site operators. Unchanged by function - copies are used. */
                                   tntComplexArray *osparam);  /* Parameters for the on-site operators. Send NULL if there are no on-site operators. */

tntNetwork tntMpsCreatePropST2sc_nn_ready(unsigned L,  /* Length of system. */
                                 tntComplex dtc,            /* Size of the time step. See the main description for information on how real and imaginary parts are applied */
                                 tntNodeArray *nn,         /* Array of nearest-neighbour operators. Send NULL if there are no nearest neighbour terms. */
                                 tntComplexArray *nnparam,  /* Array of parameters for nearest-neighbour operators. Send NULL if there are no nearest neighbour terms. */
                                 tntNodeArray *os,          /* Array of onsite operators. Send NULL if there are no on-site operators. */
                                 tntComplexArray *osparam);  /* Parameters for the on-site operators. Send NULL if there are no on-site operators. */

tntNode tntMpsCreateTwoSiteOp_nn_ready(tntNodeArray *Lnodes, /* Array of operators for the left legs */
                              tntNodeArray *Rnodes, /* Array of operators for the right legs that should correspond to those for the left operators */
                              tntComplexArray *params_os, /* Parameters for each two-site term from single-site couplings */
                              tntNodeArray *NNnodes, /* Array of operators for the nearest-neighbour couplings */
                              tntComplexArray *params_nn); /* Parameters for each two-site nearest-neighbour coupling */

void tntMpsExpecOutput_open_MPI(tntNetwork wf,  /* The current network */
                                tntMpsExOp *Op,  /* The operators for calculating expectation values */
                                unsigned printOutput,/* Set to 1 to print expectation values to screen, 0 otherwise */
                                unsigned saveOutput, /* Set to 1 to save values to output file, 0 otherwise */
                                char *savename,      /* Path to the output file. If the extension is missing ".mat" is added. Only used if saveOutput is 1. */
                                tntNode discard_legs, /* Node helping to discard legs */
                                tntNode basis, /* Node for the local basis */
                                int counter,  /* Counter for saving on different timesteps */
                                int comm_sz,  /* Total number of processes */
                                int my_rank);   /* Current process */

void tntMpsExpecThreeSites_open_MPI(tntNetwork mps,      /* The network representing the MPS. Unchanged by the function. */
                                    tntNodeArray op,     /* The operator array for calculating expectation values */
                                    char *label,         /* Label of expectation values */
                                    unsigned printOutput,/* Set to 1 to print expectation values to screen, 0 otherwise */
                                    unsigned saveOutput, /* Set to 1 to save values to output file, 0 otherwise */
                                    char *savename,      /* Path to the output file. If the extension is missing ".mat" is added. Only used if saveOutput is 1. */
                                    tntNode discard_legs, /* Node helping to discard legs */
                                    tntNode basis, /* Node for the local basis */
                                    unsigned counter,    /* can be used if saving multiple batches of expectation values e.g. for different timesteps. It appends the values to a pre-exisiting array in the position given by counter. Pass zero if a counter is not required. */
                                    int comm_sz,         /* Total number of processes */
                                    int my_rank);        /* Current process */

void network_apply_oper_single_site_mixed_obc1d(tntNetwork wf,  /*!< The current network */
                                   int site, /*!< Site at which the operator (Pauli matrix) is applied */
                                   tntNode correl_oper_TNT) /*!< Operator to be applied */
