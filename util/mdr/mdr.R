# Author: Yen-Chi Chen, Christopher R. Genovese, Larry Wasserman
# Maintainer: Yen-Chi Chen <yenchic@andrew.cmu.edu>
# Reference: Chen, Yen-Chi, Christopher R. Genovese, Ryan J. Tibshirani, and Larry Wasserman. "Nonparametric Modal Regression." arXiv preprint arXiv:1412.1716 (2014).
# Date: 04/07/2015
library(Rcpp)
library(RANN)

cppFunction('NumericVector pMSCpp(NumericMatrix data, NumericMatrix query, double h, int max_iterations , double eps){
	int n = data.nrow(), d = data.ncol(), m = query.nrow();
	int d_x = d-1;
	double dist_tmp; 
	NumericVector K(n);
	double K_tot;
	
	double pMS;
	NumericVector newPt(d);
	NumericVector oldPt(d);
	int iter_now;
	double err_now;
	
	NumericVector result(d*m);
	

	for(int w=0; w<m; w++){
		// for each point
		for(int j = 0; j<d;j++){
			newPt(j) = query(w,j);
		}
		err_now = 1e14;
		iter_now = 0; 
		
		while((iter_now < max_iterations)&&(err_now > eps)){	// ignore those nearly not change
			pMS =0.0;
			for(int j =0; j<d; j++){
				oldPt(j) = newPt(j);
			}
			
			K_tot = 0;
			for(int i = 0; i<n; i++){
				dist_tmp = 0;
				for(int j =0; j<d; j++){
					dist_tmp += (data(i,j) - newPt(j))*(data(i,j) - newPt(j))/h/h;
				}
				K(i) = exp(-1*dist_tmp/2);
				K_tot += K(i);
				pMS += data(i,d_x)*K(i);	// only update the last element (Y)
				
			}
				
			// updates & errors
			newPt(d_x) = pMS/K_tot;
			err_now = sqrt((newPt(d_x)-oldPt(d_x))* (newPt(d_x)-oldPt(d_x)));
			iter_now++;
		}
		
		for(int j =0; j<d; j++){
			result(w+m*j) = newPt(j);
		}
	}
	
	return result;

}')


#' Partial mean shift algorithm using Rcpp
#' 
#' @param data Input data matrix. The last column is Y.
#' @param query The mesh points that you want to apply mean shift to.
#' @param h Smoothing parameter.
#' @param max.iterations Maximal number of iteration for mean shift.
#' @param eps The tolerance. If mean shift moves less than this value, we will consider it done.
#' @return The mesh points after mean shift.
#' @examples
#' x = matrix(rnorm(1000), ncol=2)
#' x_pMS = pMS(x, x, 0.5)
#'  # mean shift
#' 
#' plot(x, cex=0.5)
#' points(x_MS, col="red", pch=20, cex=2)
#'  # plot the shifted result versus original case
#' 
#' @export
pMS = function(data, query=data, h, max.iterations=200, eps= 1e-15){
	tmp = pMSCpp(data=data, query=as.matrix(query), h=h, max_iterations= max.iterations, eps=eps)
	return(matrix(tmp, ncol=ncol(query)))
}

#' Dual curve dataset.
#' 
#' @param n Sample size.
#' @param err Error size.
#' @return A 2-dimensional dual curve dataset.
#' @export
DualCurve = function(n=500, err=0.25){
  X1 = runif(n)
  X2 = runif(n)
  X=c(X1,X2)
  
  Y1 =rnorm(n,1.5+0.5*sin(X1*3*pi),sd=err)
  Y2 = rnorm(n,0.5*sin(X2*3*pi),sd=err)
  Y =c(Y1,Y2)
  
  r_XY = sd(Y)/sd(X)
  Y.s = Y/r_XY
  return(cbind(X, Y.s))
}

#' Three mixture dataset.
#' 
#' @param n1 Sample size for first mixture.
#' @param n2 Sample size for second mixture.
#' @param n3 Sample size for third mixture.
#' @param err Error size.
#' @return A 2-dimensional three mixture dataset.
#' @export
ThreeMix = function(n1=100, n2=100, n3=100, err=0.2){
  x1 = runif(n1,0,0.4)
  x2 = runif(n2,0.3, 0.7)
  x3 = runif(n3,0.6,1.0)

  y1 = rnorm(n1,3,sd=err)
  y2 = rnorm(n2,2,sd=err)
  y3 = rnorm(n3,1,sd=err)

  X = c(x1,x2,x3)
  Y = c(y1,y2,y3)
  r0 = sd(Y)/sd(X)
  
  return(cbind(X,Y/r0))
}

#' Creating object "MDR".
#' @param data Input data (last column is Y).
#' @param h Smoothing parameter. 
#' @param max.iterations Maximal number of iteration for mean shift.
#' @param eps The tolerance. If mean shift moves less than this value, we will consider it done.
#' @export
MDR = function(data, h, eps=1.0e-8, max.iterations=100, ...) UseMethod("MDR")

#' @param data Input data matrix. The last column is Y.
#' @param query The mesh points that you want to apply mean shift to.
#' @param h Smoothing parameter.
#' @param cut The cut used for clustering the modal manifolds. We will use heirarchial clustering with single linkage and cut the height h*cut.
#' @param max.iterations Maximal number of iteration for mean shift.
#' @param eps The tolerance. If mean shift moves less than this value, we will consider it done.
#' @return An S4 object. A list consisting of 
#' \item{pMS}
#' {The mesh points after partial mean shift.}
#' \item{labels}
#' {The labels for each mesh points by clustering using modal regression.}
#' \item{res}
#' {The residuals for modal regression.}
#' \item{ordered}
#' {(d=2 only) A list of ordered points for each modal curve. The ordering is with respect to the covariate.}
#' @export
MDR.default = function(data, query = data, h=NULL, cut = 1, eps=1.0e-8, max.iterations=100){
  d = ncol(data)
  if(is.null(h)){
    h = (4/(d+4))^(1/(d+6))/n^(1/(d+6))*mean(apply(data,2,sd))
  }
  data.pMS = pMS(data, query=query, h)
  data.hclust = hclust(dist(data.pMS), method="single")
  data.labels = cutree(data.hclust, h=cut*h)
  
  result = list()
  result$pMS = data.pMS
  result$labels = data.labels
  result$res = query[,d] - data.pMS[,d]
  
  if(d==2){ # we can order points according to x
    result$ordered = list()
    for(j in 1:max(data.labels)){
      w_tmp = which(data.labels==j)
      x_tmp = data.pMS[w_tmp,1]
      o_tmp = order(x_tmp)
      D_tmp = data.pMS[w_tmp,][o_tmp,]
      
      result$ordered[[j]] = D_tmp
    }
  }
  return(result)
}

#' Bandwidth selection using size of prediction set (Cross-validation).
#' @param data Input data.
#' @param h_seq Sequence of bandwidth to test.
#' @param n_fold Number of folder for cross validation
#' @param n_Gx Number of grids for covariates.
#' @param n_Gy Number of grids for response.
#' @param alpha The coverage for prediction set.
#' @param h_padding The threshold size for grid points being used. 
#' We only consider the grid points whose distance to nearest data point is less than h_padding.
#' Default is to select from smoothing parameter of normal reference rule.
#' @return A matrix with first column being h_seq and second column being the size for corresponding prediction set.
#' 
#' @export
CV.h = function(data, h_seq, n_fold=5, n_Gx=50, n_Gy=100, alpha=0.9, h_padding = NULL){
  d = ncol(data)
  n = ncol(data)
  if(is.null(h_padding)){
    h_padding = (4/(d+4))^(1/(d+6))/n^(1/(d+6))*mean(apply(data,2,sd))
  }
  
  G0 = seq(from=min(data[,d]), to=max(data[,d]), length.out= n_Gy)
  for(j in (d-1):1){
    G_x = seq(from=min(data[,j]), to=max(data[,j]), length.out= n_Gx)
    G0 = expand.grid(G_x,G0)
  }
  G1 = G0[which(nn2(data,as.matrix(G0),k=1)$nn.dist<h_padding),]

  n_cv = nrow(data)/n_fold
  # number of points for tests set

  k_C = rep(NA, n_Gx)
  # number of clusters

  PI_seq_cv = rep(NA, length(h_seq))

  for(i_tmp in 1:length(h_seq)){
    h_tmp =h_seq[i_tmp]  
    idx = sample(1:nrow(data), size = nrow(data))
  
    PI_cv = rep(NA, n_fold)
    for(i_cv in 1:n_fold){
      idx_test =(1+(i_cv-1)*n_cv):(i_cv*n_cv)
    
      D_test = data[idx[idx_test],]
      D_train = data[idx[-idx_test],]
    
      D_CV.pMS = pMS(data= D_train, query = D_test, h=h_tmp)
      CV_res = D_CV.pMS[,2] - D_test[,2]
      q_RS = quantile(abs(CV_res), alpha)
    
      ## computing number of modal curves at each grid point
      G.pMS = pMS(data=D_train, query=G1, h=h_tmp)
      for(i_gr in 1:n_Gx){
        w_tmp = which(G.pMS[,1]==G_x[i_gr])
        k_C[i_gr] = length(unique(round(G.pMS[w_tmp,2], 4)))
      }
    
      PI_cv[i_cv] = mean(k_C)*q_RS
    }
    PI_seq_cv[i_tmp] = mean(PI_cv)
    # area of prediction set
    print(i_tmp)
  }
  return(cbind(h_seq, PI_seq_cv))
}
