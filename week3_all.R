# 2週目のデータから状態数を決定してから3週目の内容を行う

# パッケージ読み込み
library(markovchain)
library(MDPtoolbox)


################################################################
################################################################
# ここから2週目

# データの読み込み
source("dumpdata.R")

# データの詳細
n.item = dim(Data)[1]
n.epoch = dim(Data)[2]

if(1){ # データ全体を表示
plot(c(1,n.epoch),c(min(Data),max(Data)),
     xlab="Time",ylab="Data",type="n")
for( i in c(1:n.item) ) {
  lines(c(1:n.epoch),Data[i,])
}
}


#######################
# 決めないといけないところ
if(0){ # 状態の幅を決めて、状態数を決める
  status.min = 0.65
  status.max = 1.05
  status.width = 0.15 # 状態の幅(学生が決める箇所)
  n.status = ceiling((status.max-status.min)/status.width) # 区間数
  status.breaks = status.min + c(0:n.status)*status.width # 各区間の境界
  
  # データの離散化
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
}else{ # 状態数を決める
  status.breaks = c(0.65,0.7,0.75,0.88,0.92, 1.05)
  n.status = length(status.breaks)-1
  
  # データの離散化
  Data.states = Data
  for( i in c(1:n.status) ) {
    if(sum((Data>status.breaks[i]) & (Data<=status.breaks[i+1]))>0) {
      Data.states[ (Data>status.breaks[i]) & (Data<=status.breaks[i+1]) ] = as.character(n.status-i+1)
    }
  }
}

if(1){
  for( i in c(1:n.status) ) {
    abline(h = status.breaks[i], col = "#ff00ff80")
  }
}
if(0){
  plot(c(1,n.epoch),c(min(Data.states),max(Data.states)),
       xlab="Time",ylab="Status",type="n")
  for( i in c(1:n.item) ) {
    lines(c(1:n.epoch),Data.states[i,])
  }
}

# マルコフ連鎖
Data.mc = markovchainFit(Data.states)
P.Dgr = Data.mc$estimate@transitionMatrix
rownames(P.Dgr) = c(1:n.status)
colnames(P.Dgr) = c(1:n.status)
P.Rpr = rbind(c(1,rep(0,(n.status-1))),
              cbind(diag(1,nrow=(n.status-1),ncol=(n.status-1)),0))
rownames(P.Rpr) = c(1:n.status)
colnames(P.Rpr) = c(1:n.status)
P.Rpl = matrix(rep(c(1,rep(0,(n.status-1))),n.status),nrow=n.status,ncol=n.status,byrow=TRUE)
rownames(P.Rpl) = c(1:n.status)
colnames(P.Rpl) = c(1:n.status)


################################################################
################################################################
# ここから3週目


################################
################################
# 作業1
rownames(P.Dgr) = c(0:(n.status-1))
colnames(P.Dgr) = c(0:(n.status-1))

#状態遷移行列の生成
#状態劣化
mmdp_create.replacement.matrix = function(S) {
  R.S = max(S)-min(S)
  n.S = length(S)
  if( R.S != n.S-1 ) {
    stop("state space is not regular and/or does not begin with 0.")
  }
  P = matrix(c(1,rep(0,n.S-1)),nrow=n.S,ncol=n.S,byrow=TRUE)
  rownames(P)=S
  colnames(P)=S
  return(P)
}

#直積の状態空間Sと状態ベクトルXを与えると、一番左の数字Yを返す関数
mmdp_X.to.Y = function(X,S) {
  mdp_check.state = function(S) {
    if( min(S)!=0 ) {
      stop("state space should begin with 0.")
    } else if( max(S)!=(length(S)-1) ) {
      stop("state space should end with N.")
    } else {
      return(TRUE)
    }
  }
  m = length(X)
  n = NULL
  for( i in c(1:m) ) {
    mdp_check.state(S[[i]])
    if( sum(X[i] %in% S[[i]]) == 0 ) {
      stop(paste("invalid first state. first state should be in {",paste(S.1, sep=" ", collapse=" "),"}."));
    }
    n = append(n,length(S[[i]]))
  }
  n = cumprod(n)
  #  print(n)
  #  print(X)
  y = n[-m]*X[-1]
  y = sum(y)+X[1]
  return(y+1)
}

#その直積空間の要素数をすべて数える関数
mmdp_X.to.Y.n = function(S) {
  mdp_check.state = function(S) {
    if( min(S)!=0 ) {
      stop("state space should begin with 0.")
    } else if( max(S)!=(length(S)-1) ) {
      stop("state space should end with N.")
    } else {
      return(TRUE)
    }
  }
  m = length(S)
  x.max = NULL
  for( i in c(1:m) ) {
    mdp_check.state(S[[i]])
    x.max = append(x.max,max(S[[i]]))
  }
  return(mmdp_X.to.Y(x.max,S))
}

expand.grid(c(0:(n.status-1)),c(0:(n.status-1)))

# 状態遷移行列の組み合わせ
mmdp_expand.P.2 = function(P) {
  n.1 = dim(P[[1]])[1]
  n.2 = dim(P[[2]])[1]
  P.Y = matrix(0,
               nrow=mmdp_X.to.Y.n(list(c(0:(n.1-1)),c(0:(n.2-1)))),
               ncol=mmdp_X.to.Y.n(list(c(0:(n.1-1)),c(0:(n.2-1)))))
  for( i in c(0:(n.1-1)) ) {
    for( j in c(0:(n.2-1)) ) {
      for( k in c(0:(n.1-1)) ) {
        for( l in c(0:(n.2-1)) ) {
          P.Y[mmdp_X.to.Y(c(i,j),list(c(0:(n.1-1)),c(0:(n.2-1)))),
              mmdp_X.to.Y(c(k,l),list(c(0:(n.1-1)),c(0:(n.2-1))))] = P[[1]][i+1,k+1]*P[[2]][j+1,l+1]
        }
      }
    }
  }
  return(P.Y)
}

mmdp_expand.P.2(list(P.Dgr,P.Dgr)) # 実行例
mmdp_expand.P.2(list(P.Dgr,mmdp_create.replacement.matrix(c(0:(n.status-1)))))
mmdp_expand.P.2(list(mmdp_create.replacement.matrix(c(0:(n.status-1))),P.Dgr))
mmdp_expand.P.2(list(mmdp_create.replacement.matrix(c(0:(n.status-1))),mmdp_create.replacement.matrix(c(0:(n.status-1)))))


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

mmdp_expand.R.2(list(c(1:n.status),c(1:n.status))) # 実行例


P.Dgr.Dgr = mmdp_expand.P.2(list(P.Dgr,P.Dgr))

P.Rpl.Dgr = mmdp_expand.P.2(list(mmdp_create.replacement.matrix(c(0:(dim(P.Dgr)[1]-1))),
                                 P.Dgr))

P.Dgr.Rpl = mmdp_expand.P.2(list(P.Dgr,
                                 mmdp_create.replacement.matrix(c(0:(dim(P.Dgr)[1]-1)))))

P.Rpl.Rpl = mmdp_expand.P.2(list(
  mmdp_create.replacement.matrix(c(0:(dim(P.Dgr)[1]-1))),
  mmdp_create.replacement.matrix(c(0:(dim(P.Dgr)[1]-1)))))


P = array(0, dim=c(dim(P.Dgr.Dgr)[1],dim(P.Dgr.Dgr)[2], 4))
P[,,1] = P.Dgr.Dgr
P[,,2] = P.Rpl.Dgr
P[,,3] = P.Dgr.Rpl
P[,,4] = P.Rpl.Rpl


#######################
# 決めないといけないところ
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

optimal.policy = mdp_value_iteration(P,R,0.95) # 実行例

cbind(expand.grid(c(0:(n.status-1)),c(0:(n.status-1))),optimal.policy$policy,optimal.policy$V) # 実行例


################################
################################
# 作業2: 年齢取り換え

# 年齢取替の行動の状態遷移行列
mmdp_create.ageing.matrix = function(S) {
  R = max(S)-min(S)
  n = length(S)
  if( R != n-1 ) {
    stop("age space is not regular and/or does not begin with 0.")
  }
  max.S = max(S)
  P = rbind(cbind(0,diag(rep(1,max.S))),0)
  P[max.S+1,max.S+1] = 1
  rownames(P)=S
  colnames(P)=S
  return(P)
}

# 状態指定取替 P.Age
mmdp_create.age.replacement.matrix = function(S,T.ast) {
  R.S = max(S)-min(S)
  n.S = length(S)
  if( R.S != n.S-1 ) {
    stop("state space is not regular and/or does not begin with 0.")
  }
  P = diag(rep(1,n.S))
  P[T.ast+1,1] = 1
  P[T.ast+1,T.ast+1] = 0
  rownames(P) = S
  colnames(P) = S
  return(P)
}

mmdp_create.ageing.matrix(c(0:20))
mmdp_create.age.replacement.matrix(c(0:(n.status-1)),2)

# 故障したら取り換え必須
P.Dgr[dim(P.Dgr)[1],] = c(1,rep(0,dim(P.Dgr)[1]-1))

# 状態遷移行列の準備
S.Dgr = c(0:(dim(P.Dgr)[1]-1))
S.Age = c(0:10)
P.Age = mmdp_create.ageing.matrix(S.Age)
P.Hrd = mmdp_create.age.replacement.matrix(S.Age,4)
P.Dgr.2 = mmdp_expand.P.2(list(P.Dgr,P.Age))
P.Hrd.2 = mmdp_expand.P.2(list(mmdp_create.replacement.matrix(S.Dgr),P.Hrd))
P = array(0,dim=c(dim(P.Dgr.2)[1],dim(P.Dgr.2)[2],2))
P[,,1] = P.Dgr.2
P[,,2] = P.Hrd.2

# 費用の準備
C.Opr.Dgr = c(0,0,0,0,2000)
C.Opr.Age = rep(0,length(S.Age))
C.Rpl.Dgr = rep(150,length(S.Dgr))
C.Rpl.Age = rep(0,length(S.Age))
C.Opr.Opr = mmdp_expand.R.2(list(C.Opr.Dgr,C.Opr.Age))
C.Rpl.Rpl = mmdp_expand.R.2(list(C.Rpl.Dgr,C.Rpl.Age))
Cost = cbind(C.Opr.Opr,C.Rpl.Rpl)
R = -Cost

# 年齢が6になったら、交換を行うという保全基準の方策
cbind(expand.grid(c(0:4),c(0:10)),c(rep(1,5*6),rep(2,5*5)))

# 反復によって総期待割引き費用
mdp_eval_policy_iterative(P,R,0.95,c(rep(1,5*6),rep(2,5*5)))

# 価値反復法によって、故障したら事後取替、そして年齢による予防取替、の２つの最適な方策を求める
mdp_value_iteration(P,R,0.95)
