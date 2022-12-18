Stash
================
Noah Scheider
2022-12-16

## R Markdown

After the end of \[2. Cross-validation on PCA methods\] we will be
equipped with multiple approaches which enables us to recover the best
amount of dimensions needed to represent data “accurately” of a data
set. Intuitively, every method must have its different advantages and
draw backs. In order to get a better feeling on how our Cross-Validation
methods behave in different circumstances, we will run a simulation
study on a variety of synthesized data sets and try to cover a good
variety of covariance matrices. This enables us to discover the limits
of the different methods used and allow us to conclude and recommend
different Algorithms with respect to its underlying data. As by
construction of our Alogrithms and every single one of them besides the
\[2.1 Naïve Approach\], the \[2.5 Matrix completion method\] they are
only applicable to Multivariate Gaussian distributed data. Therefore we
restrict ourselves to this special distribution and create corresponding
data in order to get meaningful results. All in all we consider $12$
different data sets and compare the outcome of our $5$ Algorithms on
every single one of them. Furthermore, we only consider centered
Gaussian Multivariate Random variables with mean $0$. We create three
base data sets, on which we infer $4$ kinds of different noises. In
order to compare those we start with a standard Gaussian data set. The
second and third base data sets will rely on a random and structured
covariance matrix. We denote our base data sets as $D^{i}_0, D^{i}_1$
and $D^{i}_2$ and $i \in \{1,\dots,4\}$. It consists of data
$X=\{X^1,\dots,X^n\}$ and noise $\sigma$.

\[ D^{i}\_0 = X + ^{i}*p, X^j(0, *{pp}), j = 1,,n\\

D^{i}\_1 = X + ^{i}\_p, X^j(0, *1), j = 1,,n *1=MM^{T}, (M*{ij})*{1i,jp}
(0,1)\\

D^{i}\_2 = X + ^{i}\_p, X^j(0, *2), j = 1,,n *2=MM^{T},
(M*{ij})*{1i,jp}= \]

For the covariance $\Sigma_2$ of $D^{i}_2$ we will get a structured the
matrix that will emphasize the meaningful position of the last
variables. To make this more clear, we will give an example for $p=3$:

$$
M = 
\frac{1}{9}\begin{pmatrix}
  1 & 2 & 3 \\
  4 & 5 & 6 \\
  7 & 8 & 9 \\
\end{pmatrix}
\implies
\Sigma_2 = \frac{1}{81}\begin{pmatrix}
0.81 & 0.96 & 1.11 \\
0.96 & 1.15 & 1.33 \\
1.11 & 1.33 & 1.56 \\
\end{pmatrix}
$$

We infer either uniform noise, differing noise or noise that will
increase for each variable.

$$
\epsilon^{1}_p\stackrel{\text{iid}}{\sim}\mathcal{N}(0,\sigma^2\mathbb{1}_{p\times p}), \sigma^2=1.5\\
\epsilon^{2}_p\stackrel{\text{iid}}{\sim}\mathcal{N}(0,\sigma^2\mathbb{1}_{p\times p}), \sigma^2=0.2\\
\epsilon^{3}_p\stackrel{\text{iid}}{\sim}\mathcal{N}(0, U), U_{ii}\sim \mathcal{U}([0,2]) \text{ and } U_{ij}=0, \text{for } i\neq j\\
\epsilon^{4}_p\stackrel{\text{iid}}{\sim}\mathcal{N}(0, U), U_{ii}= \frac{2i}{p}\text{ and } U_{ij}=0, \text{for } i\neq j
$$
