library(wavethresh)


blocks  = DJ.EX()$blocks

#kind of RNA seq
noisy_blocks=blocks+rnorm(1024)
plot(noisy_blocks , pch=19)
points(blocks, pch=19, size=.5, col="green")
bumps  = DJ.EX()$bumps
  plot(bumps)
noisy_bumps= bumps+rnorm(1024)
#kind of ATACseq
plot(noisy_bumps , pch=19)
points(bumps, pch=19, size=.5, col="green")
