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

However, these value distributions do not, for any arbitrary dashboard, produce a closed form for the CDF of a quantal responding agent's bid distribution. The bid distribution's CDF must instead by discretely approximated.

In order to compute a discrete approximation to the bid distribution's CDF, we need to first find a suitable way to discretize our continuous value distributions, and a suitable way to discretize the dashboard mechanism itself. We also need to make sure that the discretized version of our (continuous, IC) dashboard has a discretized pricing rule that both (a) retains the same incentive compatibility properties as the continuous version and (b) retains the same functional form regardless of how finely/coarsely discretized we make the type space. How do we do this?

First, we assume that *bid space is discretized.* That is, there exists a set of bids $b_{0} < ... < b_{n}$ such that an agent faced with a dashboard can only select one of those $n + 1$ discrete options to bid.

Furthermore, we define each of these bids $b_{i}$ to have an associated *value region* $R_{i} \subseteq [0, 1]$, such that all value regions are pairwise disjoint, adjacent value regions are "connected" such that $\sup R_{i}  = \inf R_{i + 1}$, and

```math
\biguplus_{i = 0}^{n} R_{i} = [0, 1] \text{.}
```

We redefine a dashboard and pricing rule $(x, p)$ as *truthful* if for any $v \in R_{i}$,

```math
v x(b_{i}) - p(b_{i}) \geq v x(b_{j}) - p(b_{j})
```

for all $j \neq i$. In other words, a dashboard is truthful if all agents weakly prefer reporting the bid corresponding to their own value region, rather than a different value region.

We can write the local downward/upward IC constraints in this model as

```math
\forall i \in \{ 1, ..., n \}, \forall v_{i}^{\downarrow} \in R_{i}, v_{i}^{\downarrow} x(b_{i}) - p(b_{i}) \geq v_{i}^{\downarrow} x(b_{i - 1}) - p(b_{i - 1})
```

```math
\forall i \in \{ 0, ..., n - 1 \}, \forall v_{i}^{\uparrow} \in R_{i}, v_{i}^{\uparrow} x(b_{i}) - p(b_{i}) \geq v_{i}^{\uparrow} x(b_{i + 1}) - p(b_{i + 1})
```

Combining these two constraints into one proposition and rearranging so as to isolate the pricing rule, we can write

```math
\displaylines{
\forall i \in \{ 1, ..., n \}, \forall v_{i}^{\downarrow} \in R_{i}, v_{i - 1}^{\uparrow} \in R_{i - 1},
\\ v_{i}^{\downarrow} [x(b_{i}) - x(b_{i - 1})] \geq p(b_{i}) - p(b_{i - 1}) \geq v_{i - 1}^{\uparrow} [x(b_{i}) - x(b_{i - 1})] \text{.}
}
```

Observe that the upper and lower constraints are the tightest when $v_{i}^{\downarrow}$ is as small as possible, and when $v_{i - 1}^{\uparrow}$ is as large as possible. Then, summing the constraints together from $i = 1$ to $i = k$ and applying summation by parts, we obtain

```math
\displaylines{
\forall k \in \{ 1, ..., n \}, \forall v_{1}^{\downarrow} \in R_{1}, ..., v_{k}^{\downarrow} \in R_{k}, \forall v_{0}^{\uparrow} \in R_{0}, ..., v_{k - 1}^{\uparrow} \in R_{k - 1},
\\ v_{k}^{\downarrow} x(b_{k}) - v_{1}^{\downarrow} x(b_{0}) - \sum_{i = 1}^{k - 1} x(b_{i}) [v_{i + 1}^{\downarrow} - v_{i}^{\downarrow}] \geq p(b_{k}) - p(b_{0}) \geq v_{k - 1}^{\uparrow} x(b_{k}) - v_{0}^{\uparrow} x(b_{0}) - \sum_{i = 0}^{k - 2} x(b_{i + 1}) [v_{i + 1}^{\uparrow} - v_{i}^{\uparrow}]
}
```

A pricing rule that satisfies all of these constraints must also do so for the smallest possible $v_{i}^{\downarrow}\text{s}$ and the largest possible $v_{i}^{\uparrow}\text{s}$. Conversely, a pricing rule that satisfies these constraints for the smallest possible $v_{i}^{\downarrow}\text{s}$ and the largest possible $v_{i}^{\uparrow}\text{s}$ will also satisfy them for any $v_{i}^{\downarrow}\text{s}$ and $v_{i}^{\uparrow}\text{s}$. In the limit, as each $v_{i}^{\downarrow}$ and $v_{i}^{\uparrow}$ moves toward these minimal/maximal values, we obtain

```math
v_{i}^{\uparrow}, v_{i + 1}^{\downarrow} \rightarrow \sup R_{i}
```

since the value regions $R_{i}$ were assumed to be "connected" such that $\sup R_{i}  = \inf R_{i + 1}$. It turns out that this *exactly pins down the pricing rule up to an additive constant*. Using $\partial_{i}^{\uparrow} \triangleq \sup R_{i}$, we can write this expression for the pricing rule as

```math
\displaylines{
\forall k \in \{ 1, ..., n \},
\\ p(b_{k}) - p(b_{0}) = \partial_{k}^{\uparrow} x(b_{k}) - \partial_{0}^{\uparrow} x(b_{0}) - \sum_{i = 0}^{k - 1} x(b_{i + 1}) [\partial_{i + 1}^{\uparrow} - \partial_{i}^{\uparrow}]
}
```

This motivates a very natural discretization of the value space. If all the probability mass in value region $R_{i}$ were reallocated to the location $\sup R_{i}$ and $\sup R_{i} \in R_{i}$, observe that the above expression for the pricing rule would not change. Thus, we have found a dashboard and pricing rule that are incentive compatible in a discretized bid space yet continuous value space, which maintain the same functional form if the value space were then discretized to only place nonzero probabilities on each $\sup R_{i}$.

The simplest way to do this is to split the unit interval of possible values $[0, 1]$ into $n$ half-open regions, such that for all $i \in [n]$ we have the region
```math
(\frac{i - 1}{n}, \frac{i}{n}]
```

and then define one additional region for the zero type, equal to the singleton set solely containing $0 \in [0, 1]$. We then define the bids as

```math
b_{i} = \sup (\frac{i - 1}{n}, \frac{i}{n}] = i / n
```

for all $i \in [n]$ and

```math
b_{0} = \sup \{ 0 \} = 0
```

and say that in our discretized type space, an agent has probability of having value $b_{i}$ equal to $F(i / n) - F((i - 1) / n)$, where $F$ is the CDF of our original continuous value distribution.

Observe that including the existence of a zero type in this discretization ensures that $b_{0} x(b_{0}) - p(b_{0}) \leq 0$, and that this would not necessarily be the case if $b_{0} \neq 0$. Additionally, note that the IR constraint means that we must have $b_{0} x(b_{0}) - p(b_{0}) \geq 0$ in this discretized model. Thus, we know $b_{0} x(b_{0}) - p(b_{0}) = 0$, and this further simplifies our pricing rule to be

```math
\displaylines{
\forall k \in \{ 1, ..., n \},
\\ p(b_{k}) = b_{k} x(b_{k}) - \sum_{i = 0}^{k - 1} x(b_{i + 1}) [b_{i + 1} - b_{i}] \text{.}
}
```

Note that the above formula for the pricing rule also happens to remain valid for when $k = 0$.

This particular discretization of the type space also allows us to define two vital quantities: a discretized version of the conditional mirage distribution's CDF, and a discretized version of the mirage distribution's CDF. They are as follows:

```math
\mathbb{P}(\text{bid index} \leq \hat{i} \mid \text{value index} = i) = \frac{\sum\limits_{j = 0}^{\min(\hat{i}, i)} \exp(\lambda [b_{i} x(b_{j}) - p(b_{j})])}{\sum\limits_{j = 0}^{i} \exp(\lambda [b_{i} x(b_{j}) - p(b_{j})])}
```
and
```math
\mathbb{P}(\text{bid index} \leq \hat{i}) = \sum_{i = 0}^{n} \mathbb{P}(\text{bid index} \leq \hat{i} \mid \text{value index} = i) [F(\frac{i}{n}) - F(\frac{i - 1}{n})]
```

Furthermore, these can be even more easily computed by observing that the agent's utility from having value index $i$ and bidding at index $j$, which we have previously written as $b_{i} x(b_{j}) - p(b_{j})$, can be rewritten to be more numerically stable when using floating point numbers as

```math
b_{i} x(b_{j}) - p(b_{j}) = [(i - j) x(j / n) + \sum_{\ell = 1}^{j} x(\ell / n)] / n
```

Note that this formula, as with the one for the discretized conditional mirage distribution, also works when any of the indices are $0$. They can thus be fearlessly applied without having to worry about whether they're being written for the non-zero types' indices (i.e. $1, ..., n$) or the zero type's index (i.e. $0$).



> [!NOTE]  
> This section is unfinished. The plan is to document the major functions in `MirageSimulation/src/MirageSimulation.jl` so that an unfamiliar user can tinker around with them.
