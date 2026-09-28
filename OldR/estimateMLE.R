#' Estimate Functional ARCH Model
#'
#' @description est.fArch function estimates the Functional ARCH(q) model by using the Quasi-Maximum Likelihood Estimation method.
#'
#' @param fdata The functional data object with N paths.
#' @param basis The M-dimensional basis functions.
#' @param q The order of the depedence on past squared observations. If it is missing, q=1.
#'
#' @return List of model paramters:
#' @return d: d Parameter vector, for intercept function \eqn{\delta}.
#' @return As: A Matrices, for \eqn{\alpha} operators.
#'
#' @export
#'
#' @importFrom nloptr cobyla
#' @importFrom fda Data2fd create.bspline.basis eval.fd
#'
#' @details
#' This function estimates the Functional ARCH(q) model:\cr
#' \eqn{x_i(t)=\sigma_i(t)\varepsilon_i(t)}, for \eqn{t \in [0,1]} and \eqn{1\leq i \leq N},\cr
#' \eqn{\sigma_i^2(t)=\omega(t)+ \sum_{j=1}^q \int \alpha_j(t,s) x^2_{i-j}(s)ds}.
#'
#' @seealso \code{\link{est.fGarch}} \code{\link{est.fGarchx}} \code{\link{diagnostic.fGarch}}
#'
#' @examples
#' \dontrun{
#' # generate discrete evaluations of the FARCH process and smooth them into a functional data object.
#' yd = dgp.fgarch(grid_point=50, N=200, "arch")
#' yd = yd$garch_mat
#' fd = fda::Data2fd(argvals=seq(0,1,len=50),y=yd,fda::create.bspline.basis(nbasis=32))
#'
#' # extract data-driven basis functions through the truncated FPCA method.
#' basis_est = basis.est(yd, M=2, "tfpca")$basis
#'
#' # estimate an FARCH(1) model with basis when M=1.
#' arch1_est = est.fArch(fd, basis_est[,1])
#' }
#' @references
#' Aue, A., Horvath, L., F. Pellatt, D. (2017). Functional generalized autoregressive conditional heteroskedasticity. Journal of Time Series Analysis. 38(1), 3-21. <doi:10.1111/jtsa.12192>.\cr
#' Cerovecki, C., Francq, C., Hormann, S., Zakoian, J. M. (2019). Functional GARCH models: The quasi-likelihood approach and its applications. Journal of Econometrics. 209(2), 353-375. <doi:10.1016/j.jeconom.2019.01.006>.\cr
#' Hormann, S., Horvath, L., Reeder, R. (2013). A functional version of the ARCH model. Econometric Theory. 29(2), 267-288. <doi:10.1017/S0266466612000345>.\cr
est.fArch <- function(fdata, basis, q=1,NBASIS=20){
  max_eval=10000

  ## projecting the squared process onto given basis functions to get functional scores.
  basis=as.matrix(basis)
  N=dim(fdata$coefs)[2]
  grid_point=dim(basis)[1]
  M=dim(basis)[2]

  fd_eval=eval.fd(seq(0, 1, length.out=grid_point), fdata)
  y_vec=matrix(NA,N,M)
  for(i in 1:N){
    for(j in 1:M){
      y_vec[i,j]=dot_integrate(fd_eval[,i]^2*basis[,j])
    }
  }
  y_vec = as.matrix(y_vec)

  # Objective function 2
  obj_fn2_q=function(theta,ysq,M,q){

    params=get_dAq_from_theta(theta,M,q)
    d=params$d
    As=params$As

    #for cpp
    As=do.call(rbind,As)
    d=matrix(d,ncol = 1)

    return(obj_fn_q_cpp(d,As,ysq,q))
  }

  # Reformats theta (vector of parameters) into the matrices d, As and Bs.
  get_dAq_from_theta=function(theta,M,q){
    #first M entries are d
    d=theta[1:M]

    #Next p*M^2 entries are A, A_inds are indices of A matrix params
    num_mat_params=(M^2-M*(M-1)/2)
    A_inds=M+(1:(q*(num_mat_params)))

    theta_A=theta[A_inds]

    #This could be changed to search over M^2-M(M-1)/2 paramters instead
    As=list()
    for(i in 1:q){
      curr_A_vals=theta_A[1:num_mat_params+(num_mat_params)*(i-1)]
      A=matrix(0,ncol=M,nrow=M)
      diag(A)=curr_A_vals[1:M]
      A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]
      A=t(A)
      A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]

      As=append(As,list(A))
    }
    return(list(d=d,As=As))
  }

  #If univariate fix type
  if(is.null(dim(y_vec)))
    y_vec=matrix(y_vec,ncol=1)

  #num of basis functions
  M=dim(y_vec)[2]
  #num params in search space
  num_params=q*(M^2-M*(M-1)/2)+M #; num_params
  #num of params in B matrices
  num_mat_params=(M^2-M*(M-1)/2)

  #initial values for optimization

  # lower bounds
  d_lower=rep(10^-20,M)
  a_lower=rep(10^-10,q*num_mat_params)

  #constrained optimization

  #formualte objective in terms of theta only
  obj_fn3_q=function(theta){obj_fn2_q(theta=theta,ysq=y_vec,M=M,q=q)}

  #nonlinear inequality constraints on parameters
  eval_g0<-function(theta){
    M = M
    q = q
    h <- numeric(1)
    params=get_dAq_from_theta(theta, M, q)
    upb = 1-(10^-20)
    h[1] <- -(sum(Reduce("+",params$As))-upb)
    return(h)
  }

  #non linear local optimization
  #see nlopt documentation
  #run 10 times to reduce variability in optimization
  for(i in 1:10){
    init=c(runif(num_params,0,1))
    # https://github.com/astamm/nloptr/pull/157
    ress=nloptr::cobyla(x0 = init,fn = obj_fn3_q, lower=c(d_lower,a_lower),
                        upper=c(rep(1,num_params)),
                        hin = eval_g0, nl.info = FALSE,
                        control = list(xtol_rel = 1e-18, maxeval=max_eval))
    # print("1")
    if(i==1){
      curr_val=ress$value
      curr_res=ress
      theta=ress$par #; print(res)
      params=get_dAq_from_theta(theta,M,q)

      d=params$d # ; d
      if(sum(d==0)==M)
        flag=T
      else
        flag=F
    }
    else{
      theta=ress$par #; print(res)
      params=get_dAq_from_theta(theta,M,q)

      d=params$d #; d

      if(((ress$value<curr_val)&(sum(d==0)!=M))||flag){
        # print("3")
        flag=FALSE
        curr_val=ress$value
        curr_res=ress
        # print(curr_res)
      }
    }
  }
  # print("4")
  res=curr_res

  #printed in case there was any problems with the optimization
  theta=res$par #; print(res)
  params=get_dAq_from_theta(theta,M,q)

  d=params$d #; d
  As=params$As #; As

  return(params)
}


#' Estimate Functional GARCH Model
#'
#' @description est.fGarch function estimates the Functional GARCH(p,q) model by using the Quasi-Maximum Likelihood Estimation method.
#'
#' @param fdata The functional data object with N paths.
#' @param basis The M-dimensional basis functions.
#' @param p order of the depedence on past volatilities.
#' @param q order of the depedence on past squared observations.
#'
#' @return List of model paramters:
#' @return d: d Parameter vector, for intercept function \eqn{\delta}.
#' @return As: A Matrices, for \eqn{\alpha} operators.
#' @return Bs: B Matrices, for \eqn{\beta} operators.
#'
#' @export
#'
#' @import stats
#' @import sfsmisc
#' @importFrom nloptr cobyla
#' @importFrom fda Data2fd create.bspline.basis eval.fd
#'
#' @seealso \code{\link{est.fArch}} \code{\link{est.fGarchx}} \code{\link{diagnostic.fGarch}}
#'
#' @examples
#' \dontrun{
#' # generate discrete evaluations of the FGARCH process and smooth them into a functional data object.
#' yd = dgp.fgarch(grid_point=50, N=200, "garch")
#' yd = yd$garch_mat
#' fd = fda::Data2fd(argvals=seq(0,1,len=50),y=yd,fda::create.bspline.basis(nbasis=32))
#'
#' # extract data-driven basis functions through the truncated FPCA method.
#' basis_est = basis.est(yd, M=2, "tfpca")$basis
#'
#' # estimate an FGARCH(1,1) model with basis when M=1.
#' garch11_est = est.fGarch(fd, basis_est[,1])
#' }
#'
#' @references
#' Aue, A., Horvath, L., F. Pellatt, D. (2017). Functional generalized autoregressive conditional heteroskedasticity. Journal of Time Series Analysis. 38(1), 3-21. <doi:10.1111/jtsa.12192>.\cr
#' Cerovecki, C., Francq, C., Hormann, S., Zakoian, J. M. (2019). Functional GARCH models: The quasi-likelihood approach and its applications. Journal of Econometrics. 209(2), 353-375. <doi:10.1016/j.jeconom.2019.01.006>.\cr
est.fGarch <- function(fdata, basis, p=1, q=1){
  max_eval=10000

  ## projecting the squared process onto given basis functions to get functional scores.
  basis=as.matrix(basis)
  N=dim(fdata$coefs)[2]
  grid_point=dim(basis)[1]
  M=dim(basis)[2]

  # int_approx=function(mat){
  #   temp_n=NROW(mat)
  #   return((1/temp_n)*sum(mat))}
  fd_eval=eval.fd(seq(0, 1, length.out=grid_point), fdata)

  y_vec=matrix(NA,N,M)
  for(i in 1:N){
    for(j in 1:M){
      y_vec[i,j] <- dot_integrate(fd_eval[,i]^2*basis[,j])
      # y_vec[i,j]=int_approx(fd_eval[,i]^2*basis[,j])
    }
  }
  y_vec = as.matrix(y_vec)

  # Objective function 2
  obj_fn2_pq_cpp=function(theta,ysq,M,p,q){

    params=get_dABpq_from_theta(theta,M,p,q)
    d=params$d
    As=params$As
    Bs=params$Bs

    #for cpp
    As=do.call(rbind,As)
    Bs=do.call(rbind,Bs)
    d=matrix(d,ncol = 1)

    #cpp function
    return(obj_fn_pq_cpp(d,As,Bs,ysq,p,q))
  }
  # A function to transform theta

  get_dABpq_from_theta=function(theta,M,p,q){
    #first M entries are d
    d=theta[1:M]

    #Next p*M^2 entries are A, A_inds are indices of A matrix params
    num_mat_params=(M^2-M*(M-1)/2)
    A_inds=M+(1:(q*(num_mat_params)))

    #Rest are Bs..
    B_inds=sfsmisc::last(A_inds)+(1:(p*(num_mat_params)))

    theta_A=theta[A_inds]
    theta_B=theta[B_inds]

    As=list()
    for(i in 1:q){
      curr_A_vals=theta_A[1:num_mat_params+(num_mat_params)*(i-1)]
      A=matrix(0,ncol=M,nrow=M)
      diag(A)=curr_A_vals[1:M]
      A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]
      A=t(A)
      A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]

      As=append(As,list(A))
    }

    #see above comment
    Bs=list()
    for(i in 1:p){
      curr_B_vals=theta_B[1:num_mat_params+(num_mat_params)*(i-1)]
      B=matrix(0,ncol=M,nrow=M)
      diag(B)=curr_B_vals[1:M]
      B[lower.tri(B)]=curr_B_vals[(M+1):length(curr_B_vals)]
      B=t(B)
      B[lower.tri(B)]=curr_B_vals[(M+1):length(curr_B_vals)]

      Bs=append(Bs,list(B))
    }

    return(list(d=d,Bs=Bs,As=As))
  }

  #If univariate fix type
  if(is.null(dim(y_vec)))
    y_vec=matrix(y_vec,ncol=1)

  #num of basis functions
  M=dim(y_vec)[2]
  #num params in search space
  num_params=(p+q)*(M^2-M*(M-1)/2)+M; num_params
  #num of params in B and A matrices
  num_mat_params=(M^2-M*(M-1)/2)

  # lower bounds
  d_lower=rep(10^-20,M)
  a_lower=rep(10^-10,q*num_mat_params)
  b_lower=rep(10^-20,p*num_mat_params)

  #constrained optimization
  #formualte objective in terms of theta only
  obj_fn3_pq=function(theta,M,p,q){obj_fn2_pq_cpp(theta=theta,ysq=y_vec,M=M,p=p,q=q)}

  #nonlinear inequality constraints on parameters
  eval_g0<-function(theta,M,p,q){
    h <- numeric(2)
    params=get_dABpq_from_theta(theta,M,p,q)
    upb = 1-(10^-20)
    h[1] <- -(sum(Reduce("+",params$As))-upb)
    h[2] <- -(sum(Reduce("+",params$Bs))-upb)
    #h[1] <- -(sum(Reduce("+",params$As))+sum(Reduce("+",params$Bs))-upb)
    return(h)
  }

  #non linear local optimization
  #see nlopt documentation
  #run 10 times to reduce variability in optimization
  for(i in 1:10){
    init=c(runif(M,10^-10,(1-10^-10)),runif(num_params-M,0,1))
    # https://github.com/astamm/nloptr/pull/157
    ress=nloptr::cobyla(x0 = init,fn = obj_fn3_pq, lower=c(d_lower,a_lower,b_lower),
                        upper=c(rep(1,num_params)),
                        hin = eval_g0, nl.info = FALSE,
                        control = list(maxeval=max_eval), M = M, p = p,q = q)
    # print("1")
    if(i==1){
      curr_val=ress$value
      curr_res=ress
      theta=ress$par; #print(res)
      params=get_dABpq_from_theta(theta,M,p,q)

      d=params$d; d
      if(sum(d==0)==M)
        flag=T
      else
        flag=F
    }

    else{
      theta=ress$par; #print(res)
      params=get_dABpq_from_theta(theta,M,p,q)

      d=params$d; d

      if(((ress$value<curr_val)&(sum(d==0)!=M))||flag){
        # print("3")
        flag=FALSE
        curr_val=ress$value
        curr_res=ress
        # print(curr_res)
      }
    }
  }
  # print("4")
  res=curr_res

  #printed in case there was any problems with the optimization
  theta=res$par#; print(res)
  params=get_dABpq_from_theta(theta,M,p,q)

  d=params$d#; d
  As=params$As#; As
  Bs=params$Bs#; Bs

  return(params)
}

#' Estimating non-negative basis functions
#'
#' @description basis.est function generates/estimates non-negative basis functions, so that they can be used to estimate a Functional GARCH model subsequently.
#'
#' @param yd A (grid_point) x (number of observations) matrix drawn from N discrete evaluation curves.
#' @param M The number/order of basis functions, by setting \eqn{M <= 10}.
#' @param type A string to switch estimation method: "poly" - exponential and Bernstein functions; "tfpca" - truncated functional principal components; "fsnn" - sparse and non-negative functional principal components; "tpf" - truncated predictive factors.
#'
#' @return List of objects:
#' @return basis: a (grid_point) x (M) matrix containing discrete evaluations of M data-driven bases; when "poly" is applied, the basis contains a list with 'exp' for exponential and 'bern' for Bernstein.
#' @return tve: the total variations explained by data-driven bases, not applied for polynomials.
#'
#' @export
#'
#' @import fda
#' @importFrom stats time
#' @importFrom nsprcomp nsprcomp
#' @importFrom QZ qz.geigen
#'
#' @examples
#' # generate discrete evaluations of the FGARCH process.
#' yd = dgp.fgarch(50, 100, "garch")
#' yd = yd$garch_mat
#'
#' # decompose the first truncated non-negative functional principal component.
#' dt = basis.est(yd, M=1, "tfpca")
#' tbasis = dt$basis
#' tve = dt$tve
#'
#' @references
#' Cerovecki, C., Francq, C., Hormann, S., Zakoian, J. M. (2019). Functional GARCH models: The quasi-likelihood approach and its applications. Journal of Econometrics. 209(2), 353-375. <doi:10.1016/j.jeconom.2019.01.006>.\cr
#' Rice, G., Wirjanto, T., Zhao, Y. (2021) Exploring volatility of crude oil intra-day return curves: a functional GARCH-X model. MPRA Paper No. 109231. <https://mpra.ub.uni-muenchen.de/109231>.\cr
basis.est <- function(yd, M, type,NBASIS=20){
  grid_point=nrow(yd)
  N=ncol(yd)
  yd=yd-rowMeans(yd)

  switch(type,
         poly = {
           time=1:grid_point/grid_point
           exp_matrix=matrix(NA,grid_point,M)
           bernstein_matrix=matrix(NA,grid_point,M)
           for (k in 1:M){
             exp_matrix[,k]=exp(1)^(k*time)
             bernstein_matrix[,k]=(factorial(M-1)/(factorial(k-1)*factorial(M-k)))*(time^(k-1))*((1-time)^(M-k))
           }
           basis1 = list(exp = exp_matrix, bern = bernstein_matrix)
           tve = NULL
         },
         tfpca = {
           times=rep(0,grid_point)
           for(i in 1:grid_point){times[i]=i/grid_point}
           squared_y=yd*yd
           basis_obj=create.bspline.basis(rangeval=c(times[1],1),nbasis=NBASIS,norder=4)
           y_sq_fdsmooth=smooth.basis(argvals=times,y=squared_y,fdParobj=basis_obj)
           y_sq_fd=y_sq_fdsmooth$fd
           pca_obj=pca.fd(y_sq_fd,nharm=10,harmfdPar=fdPar(y_sq_fd), centerfns=TRUE)  #centerfns=FALSE.
           eigen_functs=pca_obj$harmonics
           eigen_v=pca_obj$values
           tve=eigen_v/sum(eigen_v)
           ortho_basis_matrix=matrix(0,nrow=grid_point,ncol=10)
           for(j in 1:10){indi_basis=eval.fd(times,eigen_functs[j])
           indi_basis[indi_basis<0]<-0
           ortho_basis_matrix[,j]=indi_basis}
           basis1 = ortho_basis_matrix[,1:M]
         },
         fsnn = {
           times = rep(0,grid_point)
           for(i in 1:grid_point){times[i] = i/grid_point}
           squared_y = yd*yd

           nspca.fd <- function(fdobj, nharm = 2, harmfdPar=fdPar(fdobj), centerfns = TRUE,samplesize){
             meanfd = mean.fd(fdobj)
             if (centerfns) {
               fdobj = center.fd(fdobj)
             }
             #  get coefficient matrix and its dimensions
             coef = fdobj$coefs
             coefd = dim(coef) # dim=[samplesize,nbasis]
             ndim  = length(coefd)
             nrep  = coefd[2]
             coefnames = dimnames(coef)
             if (nrep < 2) stop("PCA not possible without replications.")
             basisobj = fdobj$basis
             nbasis = basisobj$nbasis
             type = basisobj$type
             #  set up HARMBASIS
             harmbasis = harmfdPar$fd$basis
             nhbasis = harmbasis$nbasis
             #  set up LFDOBJ and LAMBDA
             Lfdobj = harmfdPar$Lfd
             lambda = harmfdPar$lambda
             #  compute CTEMP whose cross product is needed
             if (ndim == 3) {
               nvar = coefd[3]
               ctemp = matrix(0, nvar * nbasis, nrep)
               for(j in 1:nvar) {
                 index = 1:nbasis + (j - 1) * nbasis
                 ctemp[index,  ] = coef[,  , j]
               }
             } else {
               nvar = 1
               ctemp = coef
             }
             indi_pc=nsprcomp(t(ctemp),k=nbasis,nneg=TRUE,center = FALSE)
             eigvalc=(indi_pc$sdev)^2
             eigvecc=indi_pc$rotation[, 1:nharm] # check the defination of nsprcomp
             sumvecc = apply(eigvecc, 2, sum)
             eigvecc[,sumvecc < 0] =  - eigvecc[, sumvecc < 0]
             varprop = eigvalc[1:nharm]/sum(eigvalc)
             #  set up harmfd
             if (nvar == 1) {
               harmcoef = eigvecc
             } else {
               harmcoef = array(0, c(nbasis, nharm, nvar))
               for (j in 1:nvar) {
                 index = 1:nbasis + (j - 1) * nbasis
                 temp = eigvecc[index,  ]
                 harmcoef[,  , j] =  temp #slove(Wmat)%*% temp
               }
             }
             harmnames = rep("", nharm)
             for(i in 1:nharm)
               harmnames[i] = paste("PC", i, sep = "")
             if(length(coefd) == 2)
               harmnames = list(coefnames[[1]], harmnames,"values")
             if(length(coefd) == 3)
               harmnames = list(coefnames[[1]], harmnames, coefnames[[3]])
             harmfd = fd(harmcoef, harmbasis, harmnames)
             #  set up harmscr
             if (nvar == 1) {
               harmscr = inprod(fdobj, harmfd)
             } else {
               harmscr = array(0, c(nrep,   nharm, nvar))
               coefarray = fdobj$coefs
               harmcoefarray = harmfd$coefs
               for (j in 1:nvar) {
                 fdobjj= fd(as.matrix(    coefarray[,,j]), basisobj)
                 harmfdj = fd(as.matrix(harmcoefarray[,,j]), basisobj)
                 harmscr[,,j] = inprod(fdobjj, harmfdj)
               }
             }
             pcafd = list(harmfd, eigvalc, harmscr, varprop, meanfd)
             class(pcafd) = "pca.fd"
             names(pcafd) = c("harmonics", "values", "scores", "varprop", "meanfd")
             return(pcafd)
           }

           # smooth data into functional curves
           basis_obj=create.bspline.basis(rangeval=c(times[1],1),nbasis=NBASIS,norder=4)
           y_sq_fdsmooth=smooth.basis(argvals=times,y=squared_y,fdParobj=basis_obj)
           y_sq_fd=y_sq_fdsmooth$fd
           # get sparse and non-negative basis
           pca_obj=nspca.fd(y_sq_fd,nharm=10,harmfdPar=fdPar(y_sq_fd),centerfns = TRUE) # by default to compute the first 10 principal components
           eigen_functs=pca_obj$harmonics
           eigen_v=pca_obj$values
           tve=eigen_v/sum(eigen_v)
           ortho_basis_matrix=matrix(0,nrow=grid_point,ncol=10)
           for(j in 1:10){ortho_basis_matrix[,j]=eval.fd(times,eigen_functs[j])}
           basis1 = ortho_basis_matrix[,1:M]
         },
         tpf = {
           # the regularization parameter (tuning parameter)
           a = 0.75
           # the forecasting horizon
           h = 1
           est_cov<- function(yd,h){
             if (class(yd) != 'matrix'){stop('Input data type is not matrix')}
             # variance operator
             N = ncol(yd)
             C = (yd %*% t(yd)) / N
             # covaraince operator
             CC = (yd[,(h+1):N] %*% t(yd[,1:(N-h)])) / (N-h)
             return(list(C,CC))}

           K = nrow(yd)
           # variance and corvariance operator
           Temp1 = est_cov(yd,h)
           C = as.matrix(data.frame(Temp1[1]))
           CC = as.matrix(data.frame(Temp1[2]))
           # empirical operator calculation
           Ca = C + a * diag(K)
           # Find the predictive factor
           aaa = CC %*% t(CC)

           out = qz.geigen(aaa, B = Ca, only.values = FALSE)
           # eigen values and vectors
           eigen.value = out$ALPHA
           eigen.vector = out$V

           # order and sort the eigen value
           id = order(eigen.value,decreasing=TRUE)
           eigen.value = sort(eigen.value,decreasing=TRUE)
           eigen.vector.new = eigen.vector[,id]

           # Final estimation
           eigen.value.hat = eigen.value[1:M]
           eigen.vector.hat = eigen.vector.new[,1:M]

           vv = sqrt(1/(diag(t(eigen.vector.hat) %*% Ca %*% eigen.vector.hat)))
           vv.matrix = diag(0,length(vv),length(vv))
           diag(vv.matrix) = vv

           eigen.vector.hat = eigen.vector.hat %*% vv.matrix
           R = eigen.vector.hat
           FF = CC %*% R
           # truncate the negative part
           if(any(FF[,1]<0)==TRUE){
             FF[,1]=-FF[,1]}
           if(any(FF[,2]<0)==TRUE){
             FF[,2]=-FF[,2]}
           FF[FF<0]<-0

           tve=as.vector(rep(NA,length(eigen.value)))
           for (i in 1:length(eigen.value)){
             tve[i]=eigen.value[i]/sum(eigen.value)
           }
           basis1 = FF
         },
         stop("Enter something to switch me!"))
  return(list(basis = basis1, tve = tve))
}
