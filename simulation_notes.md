# simulation (best?) practices

## general

- start simple: before running time-consuming simulations, run some small examples to make sure that your code is working the way you think it should be
- do some simple benchmarks to estimate how long your code is going to take to run. Try establishing rough scaling rules: if you run your simulations at size 20, 40, and 80, how does the elapsed time and/or memory use scale up? Use that to predict how long you can expect to wait (and to establish wall-clock times if you need to run on an HPC system like a SHARCNET cluster)
- a single thoughtfully chosen simulation can tell you a lot, but you can also be misled by stochasticity ("what the heck am I doing wrong? Oh, I just happened to hit a weird sample")
- when debugging, it's good (as always) to do some hugely oversimplified simulations where you know *exactly* what should happen (e.g., a simple case where you can calculate an analytical solution to the dynamics). Or pick some extreme cases where you should *definitely* be able to see the effect you're looking for, if your expectations are correct and your code isn't buggy
- increasing the size/number of observations of your system will reduce the distracting effect of process stochasticity (your estimated parameters should be closer to the true parameters), but also takes longer. Sometimes going in the other direction (simulating *smaller* samples/worlds) can be informative.

## pseudo-random number generators (PRNGs)

- *always* set the random number seed explicitly
- pick a standard rule to set your seed (so that you and others will know that you didn't cheat by hunting for the seed that gave you good results); e.g. I usually use 101, but you could use the YYYYMMDD date that you started a project - anything so long as it's determined in advance. (If you want to be clever you can [use a hash of an arbitrary object or string](https://stackoverflow.com/questions/52334409/set-rs-random-seed-with-a-hash), but this is probably not worth the trouble ...)
- a useful trick when you are iterating over a large number of simulations and experience a problem/crash after a long runtime (e.g. on run 6664 out of 10,000) is to  `set.seed(base_seed+i)` on every iteration `i`; then you can start with interactive debugging (or verbose debugging output) at iteration `i` and know where to start the PRNG stream. (Alternately you could [dump the current state of the PRNG](http://www.cookbook-r.com/Numbers/Saving_the_state_of_the_random_number_generator/) to a checkpointing file at the beginning of each iteration, but `set.seed()` is easier ...
- these days the random number generator that comes with whatever programming language/system you're using will probably be good enough for most simulation purposes; however, it's good to know what PRNG you're actually using (in R, `RNGkind()`, see `?RNGkind` for details) (if you are doing something that requires cryptographic security you need to be **much** more careful about your choice of PRNGs)
- if you are running code/simulations in parallel you may need to be careful about how you use PRNGs  (now that you know this is a potential issue, [Google search on "R parallel RNG"](https://www.google.com/search?channel=fs&q=R+parallel+RNG) will probably get you enough information ...)
- don't reset your seed *too* frequently  (e.g. to the same value within each iteration of a simulation loop) - if you do it thoughtlessly you can end up running the *same* simulation thousands of times in a row (I've done this)

## design

Figuring out what combinations of parameters to try is like designing an experiment. Suppose you are exploring $d$ parameters (i.e., a $d$-dimensional parameter space). Some possible designs:

- *sensitivity analysis*: if you have a well-established set of "typical parameters" you can run linear transects through this point, e.g. for $d=3$ and a central point $\{x_1^*, x_2^*, x_3^*\}$ you would explore $x_1$ by simulating for parameter combinations $\{\{x_{1,\text{min}},x_2^*, x_3^*\}, \ldots, \{x_{1,\text{max}},x_2^*, x_3^*\}\}$, then the analogous sets for $x_2$ and $x_3$. This is computationally cheap (for transects of size $n$ it only takes $nd$ simulations), but doesn't tell you anything about possible interactions between parameters.
- regular grid (hypercube) of total size $n^d$ (only works if $d$ is reasonably small)
- random samples over a hypercube
- random samples over a multivariate parameter distribution (typically this is a product of independent marginal distributions, e.g. $N(\mu_1, \sigma^2_1) \times \ldots \times N(\mu_d, \sigma^2_d)$). (You can use different distributions for every margin if it seems appropriate.) The relationship between this (product of non-uniform margins) and sampling over the hypercube takes us partway down the road to *copula models* ...
- *Latin hypercube sampling* to sample (approximately) evenly over a multidimensional space (uniform or margin-product)
- *Sobol sampling* to sample (non-randomly and approximately evenly over a multidimensional space) (this approach has the nice property that successive samples 'fill in' the parameter space as you go along)


Your decisions about how many parameters to explore, how many values for each parameter, etc., will also be driven by the computational cost of a single simulation (and possibly how computational cost varies with parameter values, e.g. if you want to explore the dependence of outcomes on the size of the system itself, larger simulations will take longer ...)

### how many simulations for each parameter combination?

This depends on your goals. If you are trying to estimate power, coverage, or some other proportion-based value, you'll typically need to run at least a few hundred simulations per value (since you're usually trying to estimate these values to the precision of near a percentage point); use rules of thumb about binomial samples (standard error of a proportion = $\sqrt{p(1-p)/n}$ ) to figure out how many samples you need for a given precision. For example, if you want to estimate a proportion between 0.94 and 0.96, this translates to a standard error of approximately 0.005 (assuming that the 95% CI is $\pm 2 \sigma$), so

$$
\begin{split}
\sqrt{\frac{1}{n} \cdot p(1-p)} & = 0.005 \\
\frac{1}{n} & = \frac{0.005^2}{0.0475} \\
n & = 1900
\end{split}
$$


If you're not interested in power etc., should you run one or multiple simulations per parameter combination? It depends.  If  you want to separate *process variability* (i.e., the variation due to the stochastic processes in your simulation) from *parameter variability* (what it sounds like), then you will need multiple simulations. Otherwise, it may be sufficient to run one simulation per parameter value (which will combine/confound these two sources of variability)

## Estimation based on true model, true parameters vs more complex model

Simulating a model and then testing your estimates on the output is a best-case scenario, but probably one you should try before moving on to harder cases.

If doing serious evaluation of estimation quality you should consider simulating from a model that makes *different* assumptions from your estimation model (more complex, or even completely structurally different); you can also explore sensitivity to *particular* forms of model misspecification. Think about the fact that in real life you are always estimating in an [M-open](https://danmackinlay.name/notebook/m_open.html) scenario (in other words, the truth is always different from and more complex than your model)

Similarly, starting estimation procedures from *known, true* parameters is always the best-case scenario; worth trying first.