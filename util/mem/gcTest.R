# test the gc
# testing scenario : 
#.rs.restartR()

topInfo <- system("top -n 1 -b -u qs.chou|grep rsession")
testFunobj <- function() {testFunObjdt <- read.csv("~/data/olb-test.csv")}
#scenario 1 wo assign object name
topInfo <- system("top -n 1 -b -u qs.chou|grep rsession")
n <- 0
while(n<10) {
  testFunobj()
  n <- n+1
}
print("scenario 1: do unreferenced obj 10 times")
topInfo <- system("top -n 1 -b -u qs.chou|grep rsession")
print("after gc")
gc()
topInfo <- system("top -n 1 -b -u qs.chou|grep rsession")

print("scenario 2: do referenced obj 10 times")
topInfo <- system("top -n 1 -b -u qs.chou|grep rsession")
n <- 0
while(n<10) {
  objn <- paste0("test",n)
  assign(objn,testFunobj())
  n <- n+1
}
topInfo <- system("top -n 1 -b -u qs.chou|grep rsession")
gc()
topInfo <- system("top -n 1 -b -u qs.chou|grep rsession")

print("scenario 3: do referenced distinct obj 10 times")
topInfo <- system("top -n 1 -b -u qs.chou|grep rsession")
n <- 0
while(n<10) {
  assign(objn,testFunobj())
  n <- n+1
}

topInfo <- system("top -n 1 -b -u qs.chou|grep rsession")
gc()
topInfo <- system("top -n 1 -b -u qs.chou|grep rsession")
