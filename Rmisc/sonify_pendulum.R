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
grad(0, y1, pars)

pert <- c(0.001, 0, 0, 0)
res1 <- as.data.frame(ode(y1, seq(0, 20, by = 0.1), grad, pars))
res2 <- as.data.frame(ode(y1+pert, seq(0, 20, by = 0.1), grad, pars))

f <- function(d) with(d, cos(th1)+cos(th2)/2)

matplot(cbind(res1[,-1], res2[,-1]), type = "l",
        col = rep(1:4, 2), lty= rep(1:2, each = 4))

matplot(cbind(f(res1), f(res2)), type = "l", lty = 1:2)

## stereo = FALSE?
s1 <- with(res1, sonify(time, th1, play = FALSE))
s2 <- with(res2, sonify(time, th1, waveform = "sawtooth", play = FALSE))

## how do we combine these?? multichannel/etc. ?
tuneR::writeWave(s1, "tmp1.wav")
tuneR::writeWave(s2, "tmp2.wav")
system("xdg-open tmp1.wav")
system("xdg-open tmp2.wav")

## adding seems to work ...
s12 <- s1
s12@.Data <- s1@.Data + s2@.Data
s12 <- normalize(s12, unit = "16")
tuneR::writeWave(s12, "tmp12.wav")
system("xdg-open tmp12.wav")

## FIXME
## * pitches seem different -- different harmonics emphasized?
## * explain/understand meanings of variables (see W'pedia page)
## * normalize angles somehow/switch back to x, y coordinates?
## * distinguish pendula via octave/ left vs right channel / volume /timbre?

