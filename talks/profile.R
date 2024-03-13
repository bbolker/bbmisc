## https://www.stat.auckland.ac.nz//~ihaka/downloads/Taupo-handouts.pdf
rw2s1 <- function(n) {
    xpos = ypos = numeric(n)
    xdir = c(TRUE, FALSE)
    pm1 = c(1, -1)
    for (i in 2:n) if (sample(xdir, 1)) {
        xpos[i] = xpos[i - 1] + sample(pm1, 1)
        ypos[i] = ypos[i - 1]
    } else {
        xpos[i] = xpos[i - 1]
        ypos[i] = ypos[i - 1] + sample(pm1, 1)
    }
    list(x = xpos, y = ypos)
}

Rprof("out.out")
for (i in 1:1000) pos = rw2s1(1000)
Rprof(NULL)
summaryRprof("out.out")


rw2d2 <- function(n) {
    steps <- sample(c(-1, 1), n - 1, replace = TRUE)
    xdir <- sample(c(TRUE, FALSE), n - 1, replace = TRUE)
    xpos <- c(0, cumsum(ifelse(xdir, steps, 0)))
    ypos <- c(0, cumsum(ifelse(xdir, 0, steps)))
    return(list(x = xpos, y = ypos))
}
Rprof()
for (i in 1:100) {
    pos <- rw2d2(1e5)
}
Rprof(NULL)
summaryRprof()$by.self[1:5,]

library(profvis)
profvis(replicate(100,rw2d2(1e5)))
