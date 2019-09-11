source("/home/jm/Documents/library_dcm/utils/create_X.R")
library("evd")

eulerMasc = 0.5772156649015328

nTime = 6
nAlt = 4
nObs = 5000
nVar = 8

common = list(c("X1","X2","X3","X4"))
specific = list(c("X5"),
		c("X6"),
		c("X7"),
		c("X8"))
#beta = c(-2,2,2,2,2)
beta = c(-1, 1, 1, 1, 1)

D = list()
X = list()
U_det = list()
U = list()

for(i in 1:(nTime + 1)){
	D[[i]] = data.frame(matrix(rnorm(nObs * nVar)
			,nrow = nObs, ncol = nVar))
	D[[i]] = round(1000*D[[i]])/ 1000
	X[[i]] = create_X(common, specific, D[[i]])
	U[[i]] = matrix(X[[i]] %*% as.vector(beta)
			, nrow = nObs, ncol = nAlt, byrow = TRUE)
#	U[[i]] = U_det[[i]] + round(1000*matrix(rgumbel(nVar * nObs)
#			, nrow = nObs, ncol = nAlt))/1000
}

V = matrix(0, nrow = nObs, ncol = nTime + 1)
#V_fut = matrix(0, nrow = nObs, ncol = nTime + 1)
for(t in 1:(nTime + 1)){
	V[,t] = apply(U[[t]], 1, max)
	#V_fut[,t] = apply(U_det[[t]], 1, max)
}



onePer = matrix(0, nrow = nObs, ncol = nTime + 1)

W = matrix(0, nrow = nObs, ncol = nTime)
for(i in 1:nObs)
	for(t in 1:nTime)
		W[i,t] = max(V[i,t], V[i,t+1])

R = colMeans(V)
R = R - eulerMasc
R = R[1:nTime]

diff = W
for(i in 1:nrow(diff))
	diff[i,] = diff[i,] - R

PKeep = exp(-exp(-diff))
PKeep = round(1000*PKeep)/1000
PChange = 1 - PKeep

choices = matrix(-5, nObs, nTime)

for(t in 1:nTime)#{
#	Ut = U[[t]]
	for(i in 1:nObs)#{
		if(runif(1) < PKeep[i,t]){
			choices[i,t] = -1
		} else{
			U_act = U[[t]][i,] + rgumbel(nAlt)
			choices[i,t] = which.max(U_act) - 1
			}
#	}
#}

write.table(choices,"dynData.choices.txt"
		,quote=FALSE,row.names=FALSE,col.names=FALSE)

for(i in 1:(nTime + 1)){
	im1 = i - 1
	name = paste("dynData.",im1,".txt",sep="")
	write.table(D[[i]],name,row.names=FALSE,quote=FALSE)
}


