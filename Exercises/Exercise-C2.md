Exercise-C2
================

# 2E1

Probability of rain on Monday: (2) + (4)

# 2E2

3.  The probability that it is Monday, given that it is raining

# 2E3

The probability that it is Monday, given that it is raining: (1)

# 2E4

While water is uncountable, given the model and our prior information,
we believe that water covers 70% of the globe’s surface

# 2M1

``` r
# (1). W, W, W

p_grid <- seq(from = 0, to = 1, length.out = 20)
prior <- rep(1,20)
likelihood <- dbinom(x = 3, size = 3, prob = p_grid)
unstd.posterior <- likelihood * prior
posterior <- unstd.posterior / sum(unstd.posterior)

plot(x = p_grid, y = posterior, type = "b",
     xlab = "probability of water", ylab = "posterior probability")
mtext("20 points")
```

![](Exercise-C2_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# (2) W, W, W, L
p_grid <- seq(from = 0, to = 1, length.out = 20)
prior <- rep(1,20)
likelihood <- dbinom(x = 3, size = 4, prob = p_grid)
unstd.posterior <- likelihood * prior
posterior <- unstd.posterior / sum(unstd.posterior)

plot(x = p_grid, y = posterior, type = "b",
     xlab = "probability of water", ylab = "posterior probability")
mtext("20 points")
```

![](Exercise-C2_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# (3) L, W, W, L, W, W, W
p_grid <- seq(from = 0, to = 1, length.out = 20)
prior <- rep(1,20)
likelihood <- dbinom(x = 5, size = 7, prob = p_grid)
unstd.posterior <- likelihood * prior
posterior <- unstd.posterior / sum(unstd.posterior)

plot(x = p_grid, y = posterior, type = "b",
     xlab = "probability of water", ylab = "posterior probability")
mtext("20 points")
```

![](Exercise-C2_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

#2M2

``` r
# (1). W, W, W

p_grid <- seq(from = 0, to = 1, length.out = 20)
prior <- ifelse(p_grid < 0.5, 0, 1)
likelihood <- dbinom(x = 3, size = 3, prob = p_grid)
unstd.posterior <- likelihood * prior
posterior <- unstd.posterior / sum(unstd.posterior)

plot(x = p_grid, y = posterior, type = "b",
     xlab = "probability of water", ylab = "posterior probability")
mtext("20 points")
```

![](Exercise-C2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# (2) W, W, W, L
p_grid <- seq(from = 0, to = 1, length.out = 20)
prior <- ifelse(p_grid < 0.5, 0, 1)
likelihood <- dbinom(x = 3, size = 4, prob = p_grid)
unstd.posterior <- likelihood * prior
posterior <- unstd.posterior / sum(unstd.posterior)

plot(x = p_grid, y = posterior, type = "b",
     xlab = "probability of water", ylab = "posterior probability")
mtext("20 points")
```

![](Exercise-C2_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# (3) L, W, W, L, W, W, W
p_grid <- seq(from = 0, to = 1, length.out = 20)
prior <- ifelse(p_grid < 0.5, 0, 1)
likelihood <- dbinom(x = 5, size = 7, prob = p_grid)
unstd.posterior <- likelihood * prior
posterior <- unstd.posterior / sum(unstd.posterior)

plot(x = p_grid, y = posterior, type = "b",
     xlab = "probability of water", ylab = "posterior probability")
mtext("20 points")
```

![](Exercise-C2_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

# 2M3

$$
Pr(E\|L) = \\frac{Pr(L\|E)Pr(E)}{Pr(L)}\\\\
= \\frac{Pr(L\|E)Pr(E)}{Pr(L\|E)Pr(E) + Pr(L\|M)Pr(M)}\\\\
= \\frac{0.3\*0.5}{0.3\*0.5 + 1.0\*0.5}\\\\\\
= 0.23
$$
#2M4

$$
Pr(BB\|B) = \\frac{Pr(B\|BB)Pr(BB)}{Pr(B)}\\\\
= \\frac{Pr(B\|BB)Pr(BB)}{Pr(B\|BB)Pr(BB)+Pr(B\|BW)Pr(BW)+Pr(B\|WW)Pr(WW)}\\\\
= \\frac{1\*\\frac{1}{3}}{(1 + \\frac{1}{2} + 0)\*\\frac{1}{3}}\\\\
= \\frac{2}{3}
$$

# 2M5

$$
Pr(BB\|B) = \\frac{Pr(B\|BB)Pr(BB)}{Pr(B)}\\\\
= \\frac{Pr(B\|BB)Pr(BB)}{Pr(B\|BB)Pr(BB)+Pr(B\|BW)Pr(BW)+Pr(B\|WW)Pr(WW)}\\\\
= \\frac{1\*\\frac{2}{3}}{1\*\\frac{2}{3} + \\frac{1}{2}\*\\frac{1}{3} + 0\*\\frac{1}{3}}\\\\
= \\frac{4}{5}
$$

# 2M6

$$
Pr(BB\|B) = \\frac{Pr(B\|BB)Pr(BB)}{Pr(B)}\\\\
= \\frac{Pr(B\|BB)Pr(BB)}{Pr(B\|BB)Pr(BB)+Pr(B\|BW)Pr(BW)+Pr(B\|WW)Pr(WW)}\\\\
= \\frac{1\*\\frac{1}{6}}{1\*\\frac{1}{6} + \\frac{1}{2}\*\\frac{2}{6} + 0\*\\frac{3}{6}}\\\\
= \\frac{1}{2}
$$
# 2M7

$$
Pr(W\|BB) = Pr(W\|BB,BW)Pr(BW) + Pr(W\|BB,WW)Pr(WW)\\\\
= \\frac{1}{2}\\frac{1}{2} + 1\\frac{1}{2}\\\\=
\\frac{3}{4}
$$
# 2H1

$$
Pr(TT)= Pr(TT\|A)Pr(A) + Pr(TT\|B)Pr(B)\\\\
= 0.1^2\*0.5 + 0.2^2\*0.5=0.025
$$
# 2H2

$$
Pr(A\|TT) = \\frac{Pr(TT\|A)Pr(A)}{Pr(TT)}\\\\
= \\frac{0.1^2\*0.5}{0.025} = 0.2
$$
# 2H3

$$
Pr(A\|ST)=\\frac{Pr(ST\|A)Pr(A)}{Pr(ST)}\\\\
=\\frac{0.1\*0.9\*0.5}{0.1\*0.9\*0.5+0.2\*0.8\*0.5}\\\\
= 0.36
$$
# 2H4
$$
Pr(A\|P) = \\frac{Pr(P\|A)Pr(A)}{Pr(P)}\\\\
= \\frac{0.8\*0.5}{0.8\*0.5+0.65\*0.5}\\\\
= 0.55
$$
