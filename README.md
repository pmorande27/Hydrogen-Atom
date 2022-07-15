# Hydrogen-Atom
In one dimension, the time-independent Schrödinger equation is given by

$
\large
\begin{align}
\mathbf{H}\ \mathbf{\Psi} = E\ \mathbf{\Psi}
\end{align}
$,

where $\mathbf{H}$ is the Hamiltonian. Here, $E$ and $\mathbf{\Psi}$ are the eigenvectors and eigenvalues of $\mathbf{H}$, respectively. The Hamiltonian is expressed as

$
\Large
\begin{align}
H = -\frac{\hbar^2}{2m} \nabla^2 + V(r),
\end{align}
$

where $V(r)$ is the electric potential energy, given by

$
\Large
\begin{align}
V(r) = -\frac{e^{2}}{4 \pi \epsilon_{0} r}.
\end{align}
$

In matrix form, the Schrödinger equation is solved for N equally spaced values of r, such that r goes from ($r_{max}$/N) to $r_{max}$, where $r_{max} \sim 1.5$ nm is a sensible choice. To turn the Schrödinger equation into a matrix, $\textbf{V(r)}$ should be a diagonal matrix with the values of the potential at each r along the diagonal.
