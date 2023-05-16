llt.fun	<- function(x,z,Wmat,cmat,alpha=0) {
	m		<- nrow(cmat)
	n		<- m-1
	N		<- n*(n+1)/2
	w		<- as.matrix(Wmat^.5)
	indices.m	<- which(lower.tri(matrix(, m, m),diag=TRUE) == TRUE, arr.ind=T)
	weight	<- sqrt(2)*rep(1,nrow(indices.m))
	weight[which(indices.m[,1]==indices.m[,2])] <- 1
	d		<- diag(weight)
	w_vector 	<-  w[lower.tri(w,diag=TRUE)]
	d		<- d*w_vector
	dz		<- d%*%z
	b		<- d%*%cmat[lower.tri(cmat,diag=TRUE)]
	scaleval	<- 1/norm(b,type="F")^2

	##Returns a list of index list pairs that encodes G = L*L^T where L = lower triangular into flattened arrays l
	##    and g for the lower triangular parts of L and G, respectively. Element in the output s corresponds to g[i] and
	##    consists of two l indices, i.e.
	indices.n	<- which(lower.tri(matrix(, n, n),diag=TRUE) == TRUE, arr.ind=T)
	k.ind		<- indices.n[,1]
	l.ind		<- indices.n[,2]
	ind		<- matrix(-1,n,n)		
	counter	<- 0
	for(is in 1:nrow(indices.n)) {
		counter	<- counter+1
		ind[indices.n[is,1],indices.n[is,2]] <- counter
	}
	ind.list	<- vector("list",length=nrow(ind))
	for(i in 1:length(ind.list)) {
		ind.list[[i]] <- ind[i,][ind[i,]!= -1]
	} 
	ind.list2	<- vector("list",length=N)
	for(j in 1:N) {
		ind.list2[[j]][[1]] <- ind.list[[k.ind[j]]][1:min(k.ind[j],l.ind[j])]
		ind.list2[[j]][[2]] <- ind.list[[l.ind[j]]][1:min(k.ind[j],l.ind[j])]
	}
	out	<- numeric()
	for(i in 1:N) {
		out[i] <- sum(x[ind.list2[[i]][[1]]]* x[ind.list2[[i]][[2]]])
	}
	r		<- dz %*% out - b
	f1		<- scaleval*norm(r,type="F")^2

	##get second part of function
	xx		<- diag(n)
	xx		<- xx[lower.tri(xx,diag=TRUE)]
	diag_index	<- which(xx==1)
 	xd		<- x[diag_index]
      f2		<- -1*sum(log(xd^2)) + log(sum(x^2))

	f1+alpha*f2	

}

##gradient function 
grad.fun <- function(x,z,Wmat,cmat,alpha=0) {
	m		<- nrow(cmat)
	n		<- m-1
	N		<- n*(n+1)/2
	w		<- as.matrix(Wmat^.5)
	indices.m	<- which(lower.tri(matrix(, m, m),diag=TRUE) == TRUE, arr.ind=T)
	weight	<- sqrt(2)*rep(1,nrow(indices.m))
	weight[which(indices.m[,1]==indices.m[,2])] <- 1
	d		<- diag(weight)
	w_vector 	<- w[lower.tri(w,diag=TRUE)]
	d		<- d*w_vector
	dz		<- d%*%z
	b		<- d%*%cmat[lower.tri(cmat,diag=TRUE)]
	scaleval	<- 1/norm(b,type="F")^2

	##Returns a list of index list pairs that encodes G = L*L^T where L = lower triangular into flattened arrays l
	##    and g for the lower triangular parts of L and G, respectively. Element in the output s corresponds to g[i] and
	##    consists of two l indices, i.e.
	indices.n	<- which(lower.tri(matrix(, n, n),diag=TRUE) == TRUE, arr.ind=T)
	k.ind		<- indices.n[,1]
	l.ind		<- indices.n[,2]
	ind		<- matrix(-1,n,n)		
	counter	<- 0
	for(is in 1:nrow(indices.n)) {
		counter	<- counter+1
		ind[indices.n[is,1],indices.n[is,2]] <- counter
	}
	ind.list	<- vector("list",length=nrow(ind))
	for(i in 1:length(ind.list)) {
		ind.list[[i]] <- ind[i,][ind[i,]!= -1]
	} 
	ind.list2	<- vector("list",length=N)
	for(j in 1:N) {
		ind.list2[[j]][[1]] <- ind.list[[k.ind[j]]][1:min(k.ind[j],l.ind[j])]
		ind.list2[[j]][[2]] <- ind.list[[l.ind[j]]][1:min(k.ind[j],l.ind[j])]
	}
	out	<- numeric()
	for(i in 1:N) {
		out[i] <- sum(x[ind.list2[[i]][[1]]]* x[ind.list2[[i]][[2]]])
	}
	r		<- dz %*% out - b
	
	# Index j for which (r,r) appears in index[j].
	na			<- c(1,cumsum(n:2)+1)
	diagonal_index	<- numeric()
	for(dd in 1:n) {
		diagonal_index	<- c(diagonal_index,na[dd:n])
	}

	# All r's have a cross-term for some h_j(x) except the last (lower-right diagonal entry of L).
	##creating cross term index
	nt		<- c(1,cumsum(n:1)+1)
	nt1		<- nt[1:(length(nt)-1)]
	nt2		<- nt[2:(length(nt))]
	##r-index list
	r_list	<- numeric()
	for(i in 1:n) {
		start		<- nt1[i]
		stop		<- nt2[i]
		rng		<- range(c(start,stop))
		r_list	<- c(r_list,rep(rng[1]:(rng[2]-1),each=stop-start-1))
	}	

	##j-index list
	j_list	<- numeric()
	for(jj in 1:(n-1)) {
		j_list	<- c(j_list,(na[jj]+1):(na[jj+1]-1))
		for(s in 2:(n-jj+1))
		if((jj+s-1)<=n) {	
			if(jj+s >n) { end <- NULL
			} else {
				end <- (na[jj+s-1]+1):(na[jj+s]-1)
			}
			j_list	<- c(j_list,c(na[jj:(jj+s-2)]+s-(1:(s-1)),
							end) )
		}
	}

	##tj-index list
	na1		<- na[1:(length(na)-1)]
	na2		<- na[2:(length(na))]
	tj_list	<- numeric()
	for(i in 1:length(na1)) {
		start		<- na1[i]
		stop		<- na2[i]
		for(ss in (start):(stop-1)) {
			if((ss-1)>=start) { tj0 <- start:(ss-1)
			} else {
				tj0 <- NULL
			}
			if((ss+1)<=(stop-1)) { tjn <- (ss+1):(stop-1)
			} else {
			tjn <- NULL
			}
			tj_list	<- c(tj_list,c(tj0,tjn))
		}	
	}

	# Grad Sparsity pattern remains fixed; reuse the matrix object and only set the data in each grad_h call.
	dh_nz_cross <- as.matrix(cbind(r_list,j_list,tj_list))
	dh_nz_diag	<- as.matrix(cbind(1:N,diagonal_index))
	dh_ij		<- as.matrix(rbind(dh_nz_diag,dh_nz_cross[,-3]))
	dh_row	<- dh_ij[,1]
	dh_col	<- dh_ij[,2]

	##gradient
	# Create a Jacobian matrix with dummy non-zero elements, which are however useful for saving the
	# non-zero element ordering in the CSR's data elements below in self._perm.
	##--had to switch column and row indices to match python
	dh = sparseMatrix(j=dh_col,i=dh_row,x=1:length(dh_row), dims=c(N, N))
	cross_term_index	<- dh_nz_cross[, 3]
	# Save the CSR matrix data vector ordering.
	perm	<- dh@x 
	dh@x	<- c(as.numeric(unlist(2 * x)),as.numeric( x[cross_term_index]))[perm]

	at_r		<- t(dz) %*% r
	g1		<- as.vector(unlist(2*scaleval* dh %*% at_r))

	##get second part of function
	xx		<- diag(n)
	xx		<- xx[lower.tri(xx,diag=TRUE)]
	diag_index	<- which(xx==1)
 	xd		<- x[diag_index]
	g2		<- rep(0,N)
      g2[diag_index]	<- g2[diag_index] - (2 * 1 / xd)
      g2		<- g2 + (2 * x / sum(x^2))

	g1+alpha*g2	
}



rco_fun	<- function(optmethod="NLOPT_LD_LBFGS",zsum,cmat,Wmat,leeway_factor=1.1,num_alpha=5,alpha_init=100,alpha_step=0.01,tol=1e-6,neval=1000) {
	n		<- nrow(cmat)-1
	N		<- n*(n+1)/2
	x_init	<- diag(n)
	x_init	<- x_init[lower.tri(x_init,diag=TRUE)]
	diag_index	<- which(x_init==1)
	xd		<- x_init[diag_index]
	obj.fun	<- llt.fun(x_init,z=zsum,Wmat=Wmat,cmat=cmat,alpha=0)	
	reg.fun	<- -1*sum(log(xd^2)) + log(sum(x_init^2))
	alpha_init_scaled	<- alpha_init*abs(obj.fun)/abs(reg.fun)
	alpha_values	<- alpha_init_scaled*alpha_step^c(0:(num_alpha-1))

	x	<- x_init
	roc	<- vector("list",length=num_alpha)
	for(alph in alpha_values) {
		opt		<- nloptr(x0=x,eval_grad=grad.fun,eval_f=llt.fun,z=zsum,Wmat=Wmat,cmat=cmat,alpha=alph,opts=list("algorithm"=optmethod,xtol_rel=tol,maxeval=neval))
		roc[[which(alpha_values==alph)]][[1]]	<- opt$solution
		resid		<- llt.fun(x=opt$solution,z=zsum,Wmat=Wmat,cmat=cmat,alpha=0)^.5
		l		<- matrix(0,nrow=n,ncol=n)
		l[lower.tri(l,diag=TRUE)] <- opt$solution
		gstar		<- l%*%t(l)
		lam		<- eigen(gstar)$values
		condnum 	<- abs(max(lam))/abs(min(lam))
		roc[[which(alpha_values==alph)]][[2]]	<- c(alph,resid,condnum)
		x		<- opt$solution
	}

	roc_x			<- lapply(roc,FUN=function(x) x[[1]])
	roc_curve		<- as.data.frame(t(sapply(roc,FUN=function(x) x[[2]])))
	names(roc_curve)	<- c("alpha","residual","cond")
	index			<- which(roc_curve$residual < leeway_factor*min(roc_curve$residual))
	if(length(index)==0) {
		index	<- which.min(roc_curve$residual) 
	} else {
		index	<- index[which.min(roc_curve$cond[index])]
	}
	alpha	<- roc_curve[index,1]
	outx	<- roc_x[[index]]

	l	<- matrix(0,nrow=n,ncol=n)
	l[lower.tri(l,diag=TRUE)] <- outx
	gstar	<- l%*%t(l)
	
	out		<- vector("list",length=2)
	out[[1]]	<- gstar
	out[[2]]$info	<- roc_curve[index,]
	out[[2]]$curve	<- roc_curve
	out[[2]]$x	<- outx
	out
}


