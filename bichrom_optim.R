bichrom_optim <- 
function (phy.tree, tip.values, model = "bichrom", model.args = list(size = 28), 
          starting.val, root.dist, optim.options = list(algorithm = "NLOPT_LN_SBPLX", 
                                                        ftol_rel = 1e-08, print_level = 1, maxtime = 1.7e+08, 
                                                        maxeval = 1000)) 
{
  if (model == "bichrom") {
    results <- rep(0, 12)
    mle <- nloptr::nloptr(x0 = starting.val, eval_f = chromploid_nllike,
                          opts = optim.options, phy = phy.tree, tip.values = tip.values, 
                          pi.0 = root.dist, Q.FUN = Q_bichrom, Q.ARGS = model.args)
    print(mle)
    results[1:10] <- mle$solution
    results[11] <- mle$objective
    results[12] <- mle$status
    results <- as.data.frame(results)
    rownames(results) <- c("lambda0", "lambda1", "mu0", "mu1", "rho0", "rho1", "q01", "q10", "epsilon0", "epsilon1", "nloglike", "convergencestatus")
    
  }
  else {
    results <- rep(0, 11)
    mle <- nloptr::nloptr(x0 = starting.val, eval_f = chromploid_nllike,
                          opts = optim.options, phy = phy.tree, tip.values = tip.values, 
                          pi.0 = root.dist, Q.FUN = Q_reducedbichrom, Q.ARGS = model.args)
    print(mle)
    results[1:9] <- mle$solution
    results[10] <- mle$objective
    results[11] <- mle$status
    results <- as.data.frame(results)
    rownames(results) <- c("lambda0", "lambda1", "mu0", "mu1", "rho", "q01", "q10", "epsilon0", "epsilon1", "nloglike", "convergencestatus")
  }
    return(results)
}