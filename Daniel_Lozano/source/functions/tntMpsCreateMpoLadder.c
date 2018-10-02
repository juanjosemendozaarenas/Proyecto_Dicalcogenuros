/*
Authors: Sarah Al-Assam, Stephen Clark and Dieter Jaksch
$LastChangedDate$
(c) University of Oxford 2014
*/
/* Include the internal header for the TNT library */
//#include "tntMpsInternal.h"
#include "tnt.h"
#include "tntMps.h"
#include "../tntLadders.h"


/*! \ingroup mps
 * Creates an MPO network represeting a site-wide operator \f$\hat{O}\f$ formed from a sum of nearest-neighbour and onsite terms.
 *
 * \f[
 * \hat{O} = \sum_{j=0}^{L-2}\sum_i^{n_n}\alpha_{i,j}\hat{o}^l_{i}\otimes\hat{o}^r_i + \sum_{j=0}^{L-1}\sum_i^{n_o}\beta_{i,j}\hat{o}^s_{i}
 * \f]
 *
 * Nearest-neighbour operators \f$\hat{o}^l_{i}\f$ and \f$\hat{o}^r_i\f$ should be provided in arrays \c nnl and \c nnr respectively both having length \f$n_n\f$.
 * Onsite operators \f$\hat{o}^s_{i}\f$ should be provided in array \c os having length \f$n_o\f$.
 * The operators should be single-site operators or product MPOs, i.e. no internal legs, and two physical legs with the legs labelled as follows:
 *
 * \image html single_site_op.png
 * <CENTER> \image latex single_site_op.png ""  </CENTER>
 *
 * All the operators should have the same dimension for the physical legs.
 *
 * The parameters \f$\alpha_{i,j}\f$ and \f$\beta_{i,j}\f$ for the nearest neighbour are onsite terms are supplied in matrices \c nnparam and \c osparam respectively.
 * The matrix must have a number of rows equal to the length \f$n_n, n_o\f$ of its respective operators array, but there are two options for the number of columns:
 *
 * 1. All the parameters are uniform accross the lattice.
 * In this case the parameters array should have one column (which is the default behaviour if it is created with only one dimension specified).
 * The parameter \f$\alpha_{i,j}\f$ or \f$\beta_{i,j}\f$ should be in position \c i,\c 1  in the matrix for all sites.
 *
 * 2. One or more of the parameters can vary accross the lattice.
 * In this case the parameters matrix should have \c L-1 columns for nearest neighbour operators and \c L columns for onsite operators.
 * The parameter \f$\alpha_{i,j}\f$ or \f$\beta_{i,j}\f$ for operator \c i and site \c j should be at position \c i,\c j in the matrix.
 * Any uniform operators should have identical entries for all the sites.
 *
 * A non-product site-wide MPO is then created which represents the sum of these operators, where now each operator in the network will have non-singleton dimension internal legs.
 * The physical legs will have the same dimension as the original single-site operators, the internal legs will have a dimension equal to the number of nearest neighbour terms + 2.
 * The legs are labelled as follows:
 *
 * \image html mpo_op_sing.png
 * <CENTER> \image latex mpo_op_sing.png ""  </CENTER>
 *
 * They are connected to form the complete network:
 *
 * \image html mpo_op_nw.png
 * <CENTER> \image latex mpo_op_nw.png ""  </CENTER>
 *
 * \return The network representing the matrix product operator.
 */
tntNetwork tntMpsCreateMpoLadder(unsigned L, /*!< Length of system. */
                            tntNodeArray *nnl, /*!< Array of nearest neighbour operators for left site. Send NULL if there are no nearest neighbour terms. */
                            tntNodeArray *nnr, /*!< Array of nearest neighbour operators for right site. Send NULL if there are no nearest neighbour terms. */
                            tntComplexArray *nnparam, /*!< Array of parameters for nearest neighbour operators. Send NULL if there are no nearest neighbour terms. */
                            tntNodeArray *nnnl, /*!< Array of next nearest neighbour operators for left site. Send NULL if there are no nearest neighbour terms. */
                            tntNodeArray *nnnr, /*!< Array of next nearest neighbour operators for right site. Send NULL if there are no nearest neighbour terms. */
                            tntComplexArray *nnnparam, /*!< Array of parameters for next nearest neighbour operators. Send NULL if there are no nearest neighbour terms. */
                            tntNodeArray *os, /*!< Array of on-site operators. Send NULL if there are no on-site operators. */
                            tntComplexArray *osparam, /*!< Parameters for the on-site operators. Send NULL if there are no on-site operators. */
                            char *saveprefix,
                            unsigned ignoreQN) /*!< Optional last argument. If it is set to 1, then even if symmetries are turned on for your system, they will not be applied to this MPO */
{
    tntNetwork mponw; /* Network representing the MPO */
    tntComplexArray mpogenvals,parityvals; /* Values for one of the 'building blocks' of the MPO */
    tntNode mpogen, mpobase, mpoc, eye, opss, opssterm, opnn, parity, mpoparity; /* Nodes for building the MPO */
    tntNode Lv, Rv; /* Left and right boundary vectors for the MPO */
    unsigned Dmpo; /* Internal dimension of the MPO */
    unsigned i, j; /* Used for looping over terms and sites respectively */
    unsigned numnodes; /* The number of different nodes that need to be created */
    tntComplex prm; /* The current paramter */
    tntIntArray qnums_phys, qnums_int, qnums_op; /* Quantum numbers for the physical legs and the internal legs, and for the current operator */

    /* Check that there are the same number of left and right nearest neighbour operators */
    if ((NULL == nnl && NULL != nnr) || (NULL != nnl && NULL == nnr) || ((NULL != nnl && NULL != nnr) && nnl->sz != nnr->sz))
        tntErrorPrint("Cannot create a matrix product operator|The number of left nearest terms is not equal to the number of right nearest neighbour terms");  /* NO_COVERAGE */

    /* Check that there is a parameter for each of the left and right operators, or for each of the operators and each of the sites */
    if ((nnl != NULL) && (NULL == nnparam || nnl->sz != nnparam->numrows))
        tntErrorPrint("Cannot create a matrix product operator|The number of rows for the nearest neighbour parameters is not equal to the number of nearest neighbour terms"); /* NO_COVERAGE */

    if ((NULL == nnnl && NULL != nnnr) || (NULL != nnnl && NULL == nnnr) || ((NULL != nnnl && NULL != nnnr) && nnnl->sz != nnnr->sz))
        tntErrorPrint("Cannot create a matrix product operator|The number of left next nearest terms is not equal to the number of right nearest neighbour terms");  /* NO_COVERAGE */

    /* Check that there is a parameter for each of the left and right operators, or for each of the operators and each of the sites */
    if ((nnnl != NULL) && (NULL == nnnparam || nnnl->sz != nnnparam->numrows))
        tntErrorPrint("Cannot create a matrix product operator|The number of rows for the next nearest neighbour parameters is not equal to the number of nearest neighbour terms"); /* NO_COVERAGE */

    /* Check that there is a parameter for each of the onsite operators */
    if ((os != NULL) && (NULL == osparam || os->sz != osparam->numrows))
        tntErrorPrint("Cannot create a matrix product operator|The number of rows for on-site parameters is not equal to the number of on-site terms"); /* NO_COVERAGE */

    /* Check that there is a non-zero number of operators */
    if ((NULL == nnl || 0 == nnl->sz) && (NULL == nnnl || 0 == nnnl->sz) && (NULL == os || 0 == os->sz)) tntErrorPrint("Cannot create a matrix product operator, as there are no terms");  /* NO_COVERAGE */


    /* Check that the length of the system is more than 1 */
    if (L < 2) tntErrorPrint("Function to create MPO is only valid for a system length of two or more"); /* NO_COVERAGE */

    /* Calulcate how many nodes to create: one node if the operators are uniform, L nodes if any of the operators vary */
    if ((NULL == os || 0 == os->sz || 1 == osparam->numcols) && (NULL == nnl || 0 == nnl->sz || 1 == nnparam->numcols)) numnodes = 1;
    else numnodes = L;

    /* Create the empty network  that will contain all the MPOs and that will be returned by the function. */
    mponw = tntNetworkCreate();

    /* Create the identity node basing it on the supplied single-site operators */
    if (NULL != nnl && 0 != nnl->sz) eye = tntNodeCreateEyeOp(nnl->vals[0]);
    else eye = tntNodeCreateEyeOp(os->vals[0]);

    /* The internal dimension of the MPO will be the number of nn terms + 2. */
    Dmpo = (NULL==nnl)? nnnl->sz + 2: nnnl->sz + nnl->sz + 2;//(NULL == nnl) ? 2 : nnl->sz + 2;

    printf("DMPO=%d\n",Dmpo);
    printf("tntSymmNumGet()=%d\n",tntSymmNumGet());
    /* Parity operator created */
    parityvals = tntComplexArrayAlloc(4*4);
    parityvals.vals[0].re = parityvals.vals[4*4 - 1].re = 1.0;
    parityvals.vals[5].re = parityvals.vals[10].re = -1.0;
    parity = tntNodeCreate(&parityvals, "UD", 4, 4); /* Node represeting the parity  */


    /* ----------------- SETTING QN INFO HERE ---------------- */
    /* If symmetries are being used, allocate memory for the lists of quantum numbers for the internal legs, and assign the values using the supplied operators */


    if (tntSymmTypeGet() && (0 == ignoreQN)) {

        /* Turn off warning since we know we are going to strip QN */
        tntSysQNClearWarnOff();

        /* Get the quantum numbers for the physical leg */
        qnums_phys = tntNodeGetQN(tntSysBasisOpGet(), "D");

        /* Allocate array of the correct size for the internal leg quantum numbers - all the values will be initialised to zero */
        qnums_int = tntIntArrayAlloc(Dmpo*tntSymmNumGet());

        /* Loop through the nearest neighbour and next nearest neighbours operators, finding the correct internal quantum number for each one */
        if( (nnl != NULL) || (nnnl != NULL) ){
            for (i = 0; i < ( (nnl->sz) + (nnnl->sz) ); i++) {
                /* Get quantum number for left operator since the function calculates QN for outgoing leg, and the left operator will be associated with the outgoing leg */

                if( i % 2 == 0 ){ /* i%2==0 i < nnl->sz */
                    qnums_op = tntMpsOpGetQN(nnl->vals[i/2]);
                }
                else{
                    qnums_op = tntMpsOpGetQN(nnnl->vals[(i-1)/2]);
                    }

                /* copy it to the relevant entry for the internal leg */
                for (j = 0; j < tntSymmNumGet(); j++) qnums_int.vals[(i+1)*tntSymmNumGet() + j] = qnums_op.vals[j];

                /* free the array containing the quantum number label */
                tntIntArrayFree(&qnums_op);
            }
        }


        /* Find the quantum number for the first index in the internal leg */
        /* If there is an onsite term it will be equal to the quantum number for this */
        /* Otherwise the resultant qn of the nearest neighbour operators needs to be found */
        if (os != NULL && 0 != os->sz) {
            qnums_op = tntMpsOpGetQN(os->vals[0]);
            /* Copy the quantum numbers to the first position in the array */
            for (j = 0; j < tntSymmNumGet(); j++)  qnums_int.vals[j] = qnums_op.vals[j];

            tntIntArrayFree(&qnums_op);
        } else {
            tntIntArray qnums_op_second;

            qnums_op = tntMpsOpGetQN(nnl->vals[0]);
            qnums_op_second = tntMpsOpGetQN(nnr->vals[0]);

            /* TODO: Need to make this more general when other symmetry types are added */
            /* Add the quantum numbers together */
            for (j = 0; j < tntSymmNumGet(); j++) qnums_int.vals[j] = qnums_op.vals[j] + qnums_op_second.vals[j];


            tntIntArrayFree(&qnums_op);
            tntIntArrayFree(&qnums_op_second);
        }


    }
    /* ----------------- END OF SETTING QN INFO ----------------- */


    /* Create the generator for the MPO that doesn't depend on operators i.e. identity in first and last element */
    mpogenvals = tntComplexArrayAlloc(Dmpo*Dmpo);
    mpogenvals.vals[0].re = mpogenvals.vals[Dmpo*Dmpo - 1].re = 1.0;

    /* Create extra terms in the matrix if nnn terms appear other than the hopping terms*/
    if(nnl->sz>4){
        for(i=4;i<=Dmpo/2-2;i++){
            mpogenvals.vals[2*(i+1)+Dmpo*(2*i+1)].re=1.0; //Starting from 4 since the first non parity term is in position (10,9)
        }
    }

    /* Create node from mpogenvals */
    mpogen = tntNodeCreate(&mpogenvals, "LR", Dmpo, Dmpo);

    /* Contract the gen with the identity operator to form the base mpo term */
    /* Note that the contraction step will free both mpogen and eye */
    mpobase = tntNodeContract(mpogen, eye);

    /* Reset the numbers for the generator values */
    mpogenvals.vals[0].re = mpogenvals.vals[Dmpo*Dmpo - 1].re = 0.0;
    if(nnl->sz>4){
        for(i=4;i<=Dmpo/2-2;i++){
            mpogenvals.vals[2*(i+1)+Dmpo*(2*i+1)].re = 0.0; //Starting from 4 since the first non parity term is in position (10,9)
        }
    }

    /* Adding the parity operators within the matrix for the nnn */
    for(i=0;i<4;i++){
        mpogenvals.vals[2*(i+1)+Dmpo*(2*i+1)].re = 1.0;
    }

    mpogen = tntNodeCreate(&mpogenvals, "LR", Dmpo, Dmpo);/* Parity matrix positions in MPO  */
    mpoparity = tntNodeContract(mpogen, parity);

    tntNodeAdd(mpobase, mpoparity);  /* Identities and Parities added to the total matrix */

    /* reset the numbers for the generator values */
    for(i=0;i<4;i++){
        mpogenvals.vals[2*(i+1)+Dmpo*(2*i+1)].re = 0.0;
    }

    /* Create 1 copy of the MPO base in the network if the operators are uniform, and L copies of the MPO base in the network
       if any one of the operators vary through the system */
    /* Insert first node in the network */
    tntNodeInsertAtEnd(mpobase, "L", "R", mponw);

    /* Insert any remaining nodes in the network */
    for (j = 1; j < numnodes; j++) tntNodeInsertAtEnd(tntNodeCopy(mpobase), "L", "R", mponw);

    /* Now create the on-site terms: Create a different term for each site if the parameters vary. */
    if (os != NULL && 0 != os->sz) {

        /* Populate the bottom-left element of the MPO generator and create node, then reset array */
        mpogenvals.vals[Dmpo-1].re = 1.0;
        mpogen = tntNodeCreate(&mpogenvals, "LR", Dmpo, Dmpo);
        mpogenvals.vals[Dmpo-1].re = 0.0;

        for (j = 0; j < numnodes; j++) {

            /* Let the onsite term be equal to the first operator */
            opss = tntNodeCopy(os->vals[0]); /* First term of the OS operators */

            /* Set the parameters depending on whether the on-site parameters also vary */
            if (1 == osparam->numcols) prm = osparam->vals[0];
            else prm = osparam->vals[j * os->sz];

            /* Scale the operator by the parameter */
            tntNodeScaleComplex(opss, prm);

            /* Create each addition term, and add it to the total on-site operator */
            for (i = 1; i < os->sz; i++) {
                /* Copy operator  */
                opssterm = tntNodeCopy(os->vals[i]); /* i-th term of the OS operators */

                /* Set the parameters depending on whether the on-site parameters also vary */
                if (1 == osparam->numcols) prm = osparam->vals[i];
                else prm = osparam->vals[i + j * os->sz];

                /* Determine the correct parameter and scale */
                tntNodeScaleComplex(opssterm, prm);

                /* Add to the total term */
                tntNodeAdd(opss,opssterm);

                /* Free as it will no longer be required */
                tntNodeFree(&opssterm);
            }

            /* Create a copy of the MPO generator */
            mpoc = tntNodeCopy(mpogen);

            /* Contract the MPO generator with the on-site operator */
            mpoc = tntNodeContract(mpoc, opss);

            /* Add the MPO term to the total MPO */
            tntNodeAdd(mpobase, mpoc);

            /* free this generating node as it is no longer required */
            tntNodeFree(&mpoc);

            /* Find the next mpobase term */
            mpobase = tntNodeFindConn(mpobase, "R");

        }
        /* Free the generating node */
        tntNodeFree(&mpogen);
    }


    /* Now create the nearest neighbour terms */
    if (nnl != NULL && 0 != nnl->sz) {

        /* Create each term, using the MPO generator to insert each nn term to a different position in the MPO tensor. */

        for (i = 0; i < nnl->sz; i++) {


            /* Make the MPO generator: position of left operator will be bottom row, and (i + 1)st column (if col numbering starts from zero) */
            mpogenvals.vals[Dmpo-1 + (2*i+1)*Dmpo].re = 1.0;
            mpogen = tntNodeCreate(&mpogenvals, "LR", Dmpo, Dmpo);

            /* Reset the array for the MPO generator */
            mpogenvals.vals[Dmpo-1 + (2*i+1)*Dmpo].re = 0.0;

            /* Find the first MPO in the network */
            mpobase = tntNodeFindFirst(mponw);

            /* Now loop through sites for the left operator of the pair.
               Note there are L MPOs, but only L-1 parameters for the case of
               varying parameters. Put the parameter on the L node, which although built for the last site, will be
               destroyed by the boundary node. Therefore for the last node do not scale */

            for (j = 0; j < numnodes; j++) {
                /* Copy the left operator and scale by parameter */
                opnn = tntNodeCopy(nnl->vals[i]);

                /* Determine the correct parameter */
                if (j == L-1) {
                    prm.re = 1.0;
                    prm.im = 0.0;
                } else {
                    /* Set the parameters depending on whether the nearest neighbour parameters also vary */
                    if (1 == nnparam->numcols) prm = nnparam->vals[i];
                    else prm = nnparam->vals[i + j*nnl->sz];
                }

                tntNodeScaleComplex(opnn, prm);

                /* Make a copy of the MPO generator */
                mpoc = tntNodeCopy(mpogen);

                /* Contract the MPO generator and the single site operator */
                mpoc = tntNodeContract(mpoc, opnn);


                /* Add the MPO term to the total MPO */
                tntNodeAdd(mpobase, mpoc);


                /* free this generating node as it is no longer required */
                tntNodeFree(&mpoc);

                /* Find the next MPO in the network */
                mpobase = tntNodeFindConn(mpobase,"R");
            }

            /* free this generating node as it is no longer required */
            tntNodeFree(&mpogen);

            /* Make the MPO generator: position of right operator will be first column, and (i + 1)st row (if row numbering starts from zero) */
            mpogenvals.vals[2*i+1].re = 1.0;
            mpogen = tntNodeCreate(&mpogenvals, "LR", Dmpo, Dmpo);

            /* Reset the array for the MPO generator */
            mpogenvals.vals[2*i+1].re = 0.0;

            /* Find the first MPO in the network */
            mpobase = tntNodeFindFirst(mponw);

            /* Now loop through the sites for the right node of the pair. These
               do not need to be scaled */
            for (j = 0; j < numnodes; j++) {

                /* Copy the right operator (this does not need to be scaled by the parameter) */
                opnn = tntNodeCopy(nnr->vals[i]);

                /* Make a copy of the MPO generator */
                mpoc = tntNodeCopy(mpogen);

                /* Contract the MPO generator and the single site operator */
                mpoc = tntNodeContract(mpoc, opnn);

                /* Add the MPO term to the total MPO */
                tntNodeAdd(mpobase, mpoc);

                /* free this generating node as it is no longer required */
                tntNodeFree(&mpoc);

                /* Find the next MPO in the network */
                mpobase = tntNodeFindConn(mpobase,"R");

            }

            /* free this generating node as it is no longer required */
            tntNodeFree(&mpogen);
        }
    }
//    tntNodeFree(&mpobase);//Verificar para descomentar despues

    /* Now create the next nearest neighbour terms */

    if (nnnl != NULL && 0 != nnnl->sz) {
        /* Create each term, using the MPO generator to insert each nnn term to a different position in the MPO tensor. */
        for (i = 0; i < nnnl->sz; i++) {

            /* Make the MPO generator: position of left operator will be bottom row, and 2*(i + 1)st column (if col numbering starts from zero) */
            mpogenvals.vals[Dmpo-1 + 2*(i+1)*Dmpo].re = 1.0;
            mpogen = tntNodeCreate(&mpogenvals, "LR", Dmpo, Dmpo);

            /* Reset the array for the MPO generator */
            mpogenvals.vals[Dmpo-1 + 2*(i+1)*Dmpo].re = 0.0;

            /* Find the first MPO in the network */
            mpobase = tntNodeFindFirst(mponw);

            /* Now loop through sites for the left operator of the pair.
             Note there are L MPOs, but only L-1 parameters for the case of
             varying parameters. Put the parameter on the L node, which although built for the last site, will be
             destroyed by the boundary node. Therefore for the last node do not scale */

            for (j = 0; j < numnodes; j++) {
                /* Copy the left operator and scale by parameter */
                opnn = tntNodeCopy(nnnl->vals[i]);

                /* Determine the correct parameter */
                if (j == L-1) {
                    prm.re = 1.0;
                    prm.im = 0.0;
                } else {
                    /* Set the parameters depending on whether the nearest neighbour parameters also vary */
                    if (1 == nnnparam->numcols) prm = nnnparam->vals[i];
                    else prm = nnnparam->vals[i + j*nnnl->sz];
                }

                tntNodeScaleComplex(opnn, prm);

                /* Make a copy of the MPO generator */
                mpoc = tntNodeCopy(mpogen);

                /* Contract the MPO generator and the single site operator */
                mpoc = tntNodeContract(mpoc, opnn);

                /* Add the MPO term to the total MPO */
                tntNodeAdd(mpobase, mpoc);

                /* free this generating node as it is no longer required */
                tntNodeFree(&mpoc);

                /* Find the next MPO in the network */
                mpobase = tntNodeFindConn(mpobase,"R");

            }

            /* free this generating node as it is no longer required */
            tntNodeFree(&mpogen);

            /* free this generating node as it is no longer required */
            tntNodeFree(&mpogen);
        }
    }

    /* ----------------- SETTING QN INFO HERE ---------------- */

    /* Now put all the operators in the network in blocks form */
    if (tntSymmTypeGet() && (0 == ignoreQN)) {

        mpoc = tntNodeFindFirst(mponw);

        /* Set the quantum numbers for the first MPO. */
        tntNodeSetQN(mpoc,"D",&qnums_phys,TNT_QN_IN);
        tntNodeSetQN(mpoc,"U",&qnums_phys,TNT_QN_OUT);
        tntNodeSetQN(mpoc,"L",&qnums_int,TNT_QN_IN);
        tntNodeSetQN(mpoc,"R",&qnums_int,TNT_QN_OUT);

        while (mpoc != tntNodeFindLast(mponw)) {

            /* Move on the the next node */
            mpoc = tntNodeFindConn(mpoc,"R");

            /* Set the quantum numbers for the MPO. */
            tntNodeSetQN(mpoc,"D",&qnums_phys,TNT_QN_IN);
            tntNodeSetQN(mpoc,"U",&qnums_phys,TNT_QN_OUT);
            tntNodeSetQN(mpoc,"L",&qnums_int,TNT_QN_IN);
            tntNodeSetQN(mpoc,"R",&qnums_int,TNT_QN_OUT);
        }

        /* free the quantum numbers */
        tntIntArrayFree(&qnums_phys);
    }
    /* ----------------- END OF SETTING QN INFO ----------------- */

    /* Now if only one node has been created, create additional L-1 identical copies to put in the network */
    if (1 == numnodes) {
        mpobase = tntNodeFindFirst(mponw);
        /* Keep inserting copies of the mpo nodes at the end. */
        for (j = 1; j < L; j++) tntNodeInsertAtEnd(tntNodeCopy(mpobase), "L", "R", mponw);
    }

    /* Create the leftmost and rightmost MPO nodes by contracting with the boundary nodes */
    /* Elements for the boundary nodes - reallocate array of values */
    tntComplexArrayFree(&mpogenvals);
    mpogenvals = tntComplexArrayAlloc(Dmpo);

    /* Left boundary vector: Want to pick out the final row, so choose row vector with last element 1 */
    mpogenvals.vals[Dmpo-1].re = 1.0;
    Lv = tntNodeCreate(&mpogenvals,"RL", Dmpo);

    /* Right boundary vector: Want to pick out the first row, so choose column vector with the first element 1 */
    mpogenvals.vals[Dmpo-1].re = 0.0;
    mpogenvals.vals[0].re = 1.0;
    Rv = tntNodeCreate(&mpogenvals,"LR", Dmpo);

    /* ----------------- SETTING QN INFO HERE ---------------- */

    /* Now put the boundary vectors in blocks form */
    if (tntSymmTypeGet() && (0 == ignoreQN)) {

        mpoc = tntNodeFindFirst(mponw);

        /* Set the quantum numbers for the right boundary vector. */
        tntNodeSetQN(Rv,"L",&qnums_int,TNT_QN_IN);
        tntNodeSetQN(Rv,"R",&qnums_int,TNT_QN_OUT);

        /* Set the outgoing quantum numbers for the left boundary vector. */
        tntNodeSetQN(Lv,"R",&qnums_int,TNT_QN_OUT);

        /* Incoming boundary vectors are always zero */
        for (j = 0; j < tntSymmNumGet(); j++) qnums_int.vals[j] = 0;
        tntNodeSetQN(Lv,"L",&qnums_int,TNT_QN_IN);

        /* Now quantum numbers for internal legs are finished with - free the array */
        tntIntArrayFree(&qnums_int);

        /* Turn the warning back on */
        tntSysQNClearWarnOn();
    }
    /* ----------------- END OF SETTING QN INFO ----------------- */

    /* Insert the left and right boundary vectors in the network */
    tntNodeInsertAtStart(Lv, "L", "R", mponw);
    tntNodeInsertAtEnd(Rv, "L", "R", mponw);

    /* Contract left and right boundary vectors their neighbouring MPOs */
    mpoc = tntNodeFindConn(Lv, "R");
    Lv = tntNodeContract(Lv, mpoc);

    mpoc = tntNodeFindConn(Rv, "L");
    Rv = tntNodeContract(Rv, mpoc);

    /* Free the arrays and nodes that are no longer required */
    tntComplexArrayFree(&mpogenvals);

//    tntNodePrintAsMatrix(tntNodeFindFirst(mponw),"UL","DR");
//    //tntNodePrintAsMatrix(tntNodeFindConn(tntNodeFindFirst(mponw),"R"),"UL","DR");
//    tntNodePrintAsMatrix(tntNodeFindLast(mponw),"UL","DR");
//    //mpobase = tntNodeFindConn(tntNodeFindConn(tntNodeFindFirst(mponw),"R"),"R");
//    mpobase = tntNodeFindConn(tntNodeFindFirst(mponw),"R");
//    tntNodePrintAsMatrix(mpobase,"UL","DR");
//    tntNodesSave(saveprefix,mpobase);

    /* Return the MPO network */
    return mponw;
}
