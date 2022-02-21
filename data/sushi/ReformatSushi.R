rawsushi<-read.csv("sushi3a.5000.10.order",sep=' ',header=FALSE) #Isn't most sushi raw?

perms=as.matrix(rawsushi)[,3:12]

write.csv(perms,'SushiClean.txt')