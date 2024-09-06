## https://en.wikipedia.org/wiki/Double_pendulum
## https://stackoverflow.com/questions/68339204/double-pendulum-rk4
library(deSolve)
library(sonify)
library(tuneR)

## th1: upper angle, th2: lower angle
## om1, om2: angular velocities
## lower pendulum center of mass (distance below suspension point):
##    y = l*(cos(th1) +cos(th2)/2)

grad <- function(t, y, parms) {
    g <- with(as.list(c(y, parms)), {
        dth <- th1-th2
        num1 <- m2*l1*om1^2*sin(2*dth)+2*m2*     l2*om2^2*sin(dth)+
            2*g*m2*     cos(th2)*sin(dth)+2*g*m1*sin(th1)
        num2 <- m2*l2*om2^2*sin(2*dth)+2*(m1+m2)*l1*om1^2*sin(dth)+
            2*g*(m1+m2)*cos(th1)*sin(dth)
        den1 <- -2*l1*(m1+m2*(sin(dth))^2)
        den2 <-  2*l2*(m1+m2*(sin(dth))^2)
        c(th1=om1,
          th2=om2,
          om1=num1/den1,
          om2=num2/den2)
    })
    list(g)
}

pars <- c(g=9.8, l1=1, l2=1, m1 = 1, m2 = 1)
y0 <- c(th1=pi/4, th2=pi/4, om1=0, om2=0)
## need more extreme starting values for chaos)
y1 <- c(th1=3*pi/4, th2=-3*pi/4, om1=0, om2=0)
## grad(0, y1, pars)

pert <- c(0.001, 0, 0, 0)
res1 <- as.data.frame(ode(y1, seq(0, 20, by = 0.1), grad, pars))
res2 <- as.data.frame(ode(y1+pert, seq(0, 20, by = 0.1), grad, pars))

calc_y <- function(d) with(d, cos(th1)+cos(th2)/2)

matplot(cbind(res1[,-1], res2[,-1]), type = "l",
        col = rep(1:4, 2), lty= rep(1:2, each = 4))

matplot(cbind(calc_y(res1), calc_y(res2)), type = "l", lty = 1:2)
tvec <- res1$time

## stereo = FALSE?
s1 <- sonify(tvec, calc_y(res1), waveform = "square", play = FALSE)
s2 <- sonify(tvec, calc_y(res2), waveform = "sawtooth", play = FALSE)

## how do we combine these?? multichannel/etc. ?
tuneR::writeWave(s1, "tmp1.wav")
tuneR::writeWave(s2, "tmp2.wav")
system("xdg-open tmp1.wav")
system("xdg-open tmp2.wav")

s_comb <- function(s1, s2) {
    ## adding seems to work ...
    s12 <- s1
    s12@.Data <- s1@.Data + s2@.Data
    return(normalize(s12, unit = "16"))
}

s_lrcomb <- function(s1, s2) {
    ## adding seems to work ...
    s12 <- s1
    s12@.Data[] <- cbind(s1@.Data[,"FR"], s2@.Data[,"FL"])
    return(normalize(s12, unit = "16"))
}

tuneR::writeWave(s_comb(s1, s2), "timbres.wav")
system("xdg-open timbres.wav")

s1L <- sonify(tvec, calc_y(res1), waveform = "square", play = FALSE,
             stereo = FALSE)
s2R <- sonify(tvec, calc_y(res2), waveform = "sawtooth", play = FALSE,
              stereo = FALSE)
tuneR::writeWave(s_comb(s1L, s2R), "timbres_LR.wav")
system("xdg-open timbres_LR.wav")
## * distinguish pendula via octave/ left vs right channel / volume /timbre?

s1L_low <- sonify(tvec, calc_y(res1), waveform = "square", play = FALSE,
             stereo = FALSE, flim = c(220, 440))
tuneR::writeWave(s_comb(s1L_low, s2R), "timbres_LR_oct.wav")
system("xdg-open timbres_LR_oct.wav")

