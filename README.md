# Simulating the Mirage Distribution
Part of the work on quantal response dashboards.

In this repository, we're going to be simulating how prices change the mirage distribution for a given value distribution.

Running the simulation is performed by navigating to the project directory in a terminal and running the command

```
julia init.jl
```

To do this, you'll first need to install Julia. Go to the [Julia downloads page](https://julialang.org/downloads/#current_stable_release) and select the download type corresponding to your operating system. During the installation process, select the option to add Julia to the path. *If you do not select the option to add Julia to the path during the installation process, you will have to do so manually afterwards.*

Once the process is complete, you should be able to run the simulation from your preferred terminal emulator.

For basics on the Julia language, see the [getting started page](https://docs.julialang.org/en/v1/manual/getting-started/).

> [!TIP]
> We've carefully written the `init.jl` script so that you can run it from any directory. If you wanted to, you could use `julia mirage-distribution/init.jl` from the parent directory of this GitHub repository and it would still work.


## The Model
Inside `MirageSimulation/src/MirageSimulation.jl`, you'll find two different functions for the CDF of the value distribution: the Kumaraswamy distribution and the beta distribution. In the code, these are named `kumaraswamyCDF` and `betaCDF`, respectively. These both have a continuous support on the interval $[0, 1]$ and are vaguely bell-curve-shaped, making them ideal candidates for a flexible family of agent value distributions.

However, in order to compute the distribution of agents' bids, we need to first discretize these continuous value distributions. How do we do this?

First, we assume that *bid space is discretized.* That is, there exists a set of bids $b_{0}, ..., b_{n}$ such that an agent faced with a dashboard can only select one of those $n + 1$ discrete options to bid.

Furthermore, we define each of these bids $b_{i}$ to have an associated *value region* $R_{i} \subseteq [0, 1]$, such that all value regions are pairwise disjoint and

$$\biguplus_{i = 0}^{n} R_{i} = [0, 1] \text{.}$$



> [!NOTE]  
> This section is unfinished. The plan is to document the major functions in `MirageSimulation/src/MirageSimulation.jl` so that an unfamiliar user can tinker around with them.
