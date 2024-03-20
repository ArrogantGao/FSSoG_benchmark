# Accuracy Results

Summarize the accuracy results, for parameters choosing.

## Short-range part

The following figure shows the decaying rate of short range interaction with different preset parameters.
For all five cases, the absolute error decay to minimum at about $r_c = 10$, which is a proper short range cutoff.

![](figs/accuracy_short_naive.png)

To verify that, we take a random system, assume that there are $1000$ particles in the system, we can find that the error of total short-range interaction energy vary as follow.
Here we compared the summation error of different cutoffs, comparing to the case where $r_c = L / 2$, where $L$ is the box length.
The results indicate that the error is minimized when $r_c = 10$.

![](figs/accuracy_short_naive_absolute.png)

To further reduce the cost of calculating the short range part, we can use the Chebyshev interpolation to approximate the sum of Gaussians.
The error is given as below, which indicates that terms of Chebyshev polynomial to reach the high accuracy on this finite region $x \in [0.5, 10]$ is much smaller that the naive summation.

![](figs/accuracy_short_cheb.png)

Then fix cutoff as $r_c = 10$, we can compare the error of Chebyshev interpolation with different number of terms.
It is shown that with a low order of Chebyshev polynomial, the error is fine, and the correspond parameters are given in the form.

![](figs/accuracy_short_cheb_total.png)

|  preset   | Cheb Order  |
|  ----  | ----  |
| 1  | 6 |
| 2  | 8 |
| 3  | 10 |
| 4  | 20 |
| 5  | 26 |

