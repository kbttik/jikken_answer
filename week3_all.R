# 2�T�ڂ̃f�[�^�����Ԑ������肵�Ă���3�T�ڂ̓��e���s��

# �p�b�P�[�W�ǂݍ���
library(markovchain)
library(MDPtoolbox)

# �f�[�^�̓ǂݍ���
source("dumpdata.R")

# �f�[�^�̏ڍ�
n.item = dim(Data)[1]
n.epoch = dim(Data)[2]

if(1){ # �f�[�^�S�̂�\��
plot(c(1,n.epoch),c(min(Data),max(Data)),
     xlab="Time",ylab="Data",type="n")
for( i in c(1:n.item) ) {
  lines(c(1:n.epoch),Data[i,])
}
}

if(0){ # ��Ԃ̕������߂āA��Ԑ������߂�
  status.min = 0.65
  status.max = 1.05
  status.width = 0.15 # ��Ԃ̕�(�w�������߂�ӏ�)
  n.status = ceiling((status.max-status.min)/status.width) # ��Ԑ�
  status.breaks = status.min + c(0:n.status)*status.width # �e��Ԃ̋��E
  
  # �f�[�^�̗��U��
  Data.states = Data
  for( i in c(1:n.status) ) {
    if(sum((Data>status.breaks[i]) & (Data<=status.breaks[i+1]))>0) {
      Data.states[ (Data>status.breaks[i]) & (Data<=status.breaks[i+1]) ] = as.character(n.status-i+1)
    }
  }
  
  plot(c(1,n.epoch),c(min(Data.states),max(Data.states)),
       xlab="Time",ylab="Status",type="n")
  for( i in c(1:n.item) ) {
    lines(c(1:n.epoch),Data.states[i,])
  }
}else{ # ��Ԑ������߂�
  status.breaks = c(0.65,0.7,0.75,0.88,0.92, 1.05)
  n.status = length(status.breaks)-1
  
  # �f�[�^�̗��U��
  Data.states = Data
  for( i in c(1:n.status) ) {
    if(sum((Data>status.breaks[i]) & (Data<=status.breaks[i+1]))>0) {
      Data.states[ (Data>status.breaks[i]) & (Data<=status.breaks[i+1]) ] = as.character(n.status-i+1)
    }
  }
  plot(c(1,n.epoch),c(min(Data.states),max(Data.states)),
       xlab="Time",ylab="Status",type="n")
  for( i in c(1:n.item) ) {
    lines(c(1:n.epoch),Data.states[i,])
  }
}

# �}���R�t�A��
Data.mc = markovchainFit(Data.states)
P.Dgr = Data.mc$estimate@transitionMatrix
rownames(P.Dgr) = c(1:n.status)
colnames(P.Dgr) = c(1:n.status)
P.Rpr = rbind(c(1,rep(0,n.status-1)),
              cbind(diag(1,nrow=n.status-1,ncol=n.status-1),0))
rownames(P.Rpr) = c(1:n.status)
colnames(P.Rpr) = c(1:n.status)
P.Rpl = matrix(rep(c(1,rep(0,n.status-1)),n.status),nrow=n.status,ncol=n.status,byrow=TRUE)
rownames(P.Rpl) = c(1:n.status)
colnames(P.Rpl) = c(1:n.status)


################################
################################
# ��������3�T��

rownames(P.Dgr) = c(0:4)
colnames(P.Dgr) = c(0:4)

expand.grid(c(0:4),c(0:4))

mmdp_expand.P.2(list(P.Dgr,P.Dgr)) # ���s��

mmdp_expand.P.2(list(P.Dgr,mmdp_create.replacement.matrix(c(0:4))))
mmdp_expand.P.2(list(mmdp_create.replacement.matrix(c(0:4)),P.Dgr))
mmdp_expand.P.2(list(mmdp_create.replacement.matrix(c(0:4)),mmdp_create.replacement.matrix(c(0:4))))


mmdp_expand.R.2 = function(R) {
  n.1 = length(R[[1]])
  n.2 = length(R[[2]])
  R.Y = rep(0,mmdp_X.to.Y.n(list(c(0:(n.1-1)),c(0:(n.2-1)))))
  for( i in c(0:(n.1-1)) ) {
    for( j in c(0:(n.2-1)) ) {
      R.Y[mmdp_X.to.Y(c(i,j),list(c(0:(n.1-1)),c(0:(n.2-1))))] = R[[1]][i+1]+R[[2]][j+1]
    }
  }
  return(R.Y)
}

mmdp_expand.R.2(list(c(1,2,3,4,5),c(1,2,3,4,5))) # ���s��


P.Dgr.Dgr = mmdp_expand.P.2(list(P.Dgr,P.Dgr))

P.Rpl.Dgr = mmdp_expand.P.2(list(mmdp_create.replacement.matrix(c(0:(dim(P.Dgr)[1]-1))),
                                 P.Dgr))

P.Dgr.Rpl = mmdp_expand.P.2(list(P.Dgr,
                                 mmdp_create.replacement.matrix(c(0:(dim(P.Dgr)[1]-1)))))

P.Rpl.Rpl = mmdp_expand.P.2(list(
  mmdp_create.replacement.matrix(c(0:(dim(P.Dgr)[1]-1))),
  mmdp_create.replacement.matrix(c(0:(dim(P.Dgr)[1]-1)))))


P = array(0, dim=c(dim(P.Dgr.Dgr)[1],dim(P.Dgr.Dgr)[2],4))
P[,,1] = P.Dgr.Dgr
P[,,2] = P.Rpl.Dgr
P[,,3] = P.Dgr.Rpl
P[,,4] = P.Rpl.Rpl


C.Opr = c(0,0,0,0,2000)
C.Rpl = c(150,150,150,150,150)
C.Opr.Opr = mmdp_expand.R.2(list(C.Opr,C.Opr))
C.Rpl.Opr = mmdp_expand.R.2(list(C.Rpl,C.Opr))
C.Opr.Rpl = mmdp_expand.R.2(list(C.Opr,C.Rpl))
C.Rpl.Rpl = mmdp_expand.R.2(list(C.Rpl,C.Rpl))
Cost = cbind(C.Opr.Opr,C.Rpl.Opr,C.Opr.Rpl,C.Rpl.Rpl)
R = -Cost

library(MDPtoolbox)

mdp_value_iteration(P,R,0.95)

optimal.policy = mdp_value_iteration(P,R,0.95) # ���s��

cbind(expand.grid(c(0:4),c(0:4)),optimal.policy$policy,optimal.policy$V) # ���s��