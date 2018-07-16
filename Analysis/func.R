## Script contains the necessary functions to plot various pavo objects 

## Script for JNDbar plots. Note that this functin requires the proper formatting of the data object to be a matrix. By defaut it will assume that all related bars are stacked within columns. The same is true for the error argument. Make use of ddply and plotrix on pavo objects to take averages and standard errors of the various dS and dL columns.

JNDBarplot <- function(data, name = "region", error = 1, pos, fontsize = 2, ...){
	barplot(data, beside = TRUE, legend.text = rownames(data), args.legend = list(bty = "n"), ...) -> bp.out
	abline(h = 0)
	arrows(x0 = bp.out, y0 = data, x1 = bp.out, y1 = data+error, length = 0.05, angle = 90)
	arrows(x0 = bp.out, y0 = data, x1 = bp.out, y1 = data-error, length = 0.05, angle = 90)
	text(name, x = max(bp.out)/2 + 0.5, y = pos, cex = fontsize)
	bp.out
}

## Function below will plot various spectral curves of aggregated spectral reflectance for grouping variables.

plotcol <- function(x, error, x0, x1, y0, y1, region, ...){
	plot(x[,"f"]~x[,1], col = "brown", type = "l", ...)
	lines(x[,"f"] + error[,"f"] ~x[,1], col = "brown", lty = 2)
	lines(x[,"f"] - error[,"f"] ~x[,1], col = "brown", lty = 2)

	lines(x[,"m"]~x[,1], col = "blue")
	lines(x[,"m"] + error[,"m"] ~x[,1], col = "blue", lty = 2)
	lines(x[,"m"] - error[,"m"] ~x[,1], col = "blue", lty = 2)

	#legend
	arrows(x0, y0, x1, y1, length = 0, col = "brown", lwd = 2)
	arrows(x0, y0-3, x1, y1-3, length = 0, col = "blue", lwd = 2)
	text("\\VE", vfont = c("sans serif", "bold"), xpd = TRUE, x = x1+12, y = y1, cex = 2)
	text("\\MA", vfont = c("sans serif", "bold"), xpd = TRUE, x = x1+12, y = y1-3, cex = 2)

	text(x = 500, y = 45, region, cex = 2)
}

