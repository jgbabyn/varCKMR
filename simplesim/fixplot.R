sdf = readRDS("./sameVsDiff.rds")

pdf(file="sameVsDiff3.pdf")
plot(sdf$nus,sdf$same_case,type="l",lty=1,xlab=expression(paste("log(",nu,")")),
     ylab = "Expected number of pairs sharing a mother",
     main=expression(paste("Expected number of pairs sharing a mother under varying ",nu)))
lines(sdf$nus,sdf$same_case_theo,type="l",lty=2,lwd=3)
lines(sdf$nus,sdf$diff_case_theo,type="l",lty=3,lwd=3)
lines(sdf$nus,sdf$diff_case,type="l",lty=4)
legend(0,3,legend =c("Within-cohort theoretical",
                     "Within-cohort simulated",
                     "Between-cohort theoretical",
                     "Between-cohort simulated"),
       lty=c(2,1,3,4),lwd=c(3,1,3,1))
dev.off()


