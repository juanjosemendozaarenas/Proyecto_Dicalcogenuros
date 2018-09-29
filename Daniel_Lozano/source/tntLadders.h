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
                            unsigned ignoreQN); /*!< Optional last argument. If it is set to 1, then even if symmetries are turned on for your system, they will not be applied to this MPO */

