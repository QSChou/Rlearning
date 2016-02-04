## use par("mar") to check the current margin setting, then try to adjust the margin by par(mar=c(3,1,3,1))
par("mar")
par(mar=c(3,1,3,1)
## use par("mfrow") to check the current layout, then use par(mfrow = c(3,1)) to do partition
par("mfrow")
par(mfrow = c(1,3))

## Plot with point lable
attach(mtcars)
plot(wt, mpg, main="Milage vs. Car Weight", xlab="Weight", ylab="Mileage", pch=18, col="blue")
text(wt, mpg, row.names(mtcars), cex=0.6, pos=4, col="red")
## Case when figure margins too large

## Plot the value distribution
attach(iris)
plot(table(Sepal.Length)
