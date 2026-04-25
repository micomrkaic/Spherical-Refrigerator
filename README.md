# Coupled Spherical Heat Transfer: Water, Glass, Stagnant Air, and Convection

This project solves transient heat transfer for a spherical glass container filled with warm water and placed in a cold fridge.

The key point is that the **water is not lumped**. The model allows the water temperature to vary with radius:

\[
T_w = T_w(\rho,t).
\]

Two exterior models are compared:

1. **Stagnant-air conduction model**  
   Water, glass, and fridge air are all treated as radial heat-conduction PDE domains.

2. **Robin convection model**  
   Water and glass are treated as radial PDE domains, while the external fridge air is represented by a convective heat-transfer boundary condition at the outer glass surface.

The Robin model is usually the more realistic approximation for an actual fridge.

---

## 1. Geometry

Let \(\rho\) denote radial position.

The water radius is

\[
r.
\]

The glass thickness is

\[
d.
\]

The outer glass radius is

\[
b = r+d.
\]

The fridge-wall radius used in the stagnant-air model is

\[
R.
\]

The regions are:

\[
0 < \rho < r
\]

for water,

\[
r < \rho < b
\]

for glass, and, in the stagnant-air model only,

\[
b < \rho < R
\]

for fridge air.

We assume

\[
R>b.
\]

---

## 2. Temperature Variables

Let

\[
T(\rho,t)
\]

be temperature.

The initial water temperature is

\[
T_0.
\]

The fridge bulk temperature is

\[
T_f.
\]

Usually,

\[
T_0 > T_f.
\]

It is convenient to define excess temperature:

\[
\theta(\rho,t)=T(\rho,t)-T_f.
\]

Then the final equilibrium is

\[
\theta(\rho,\infty)=0.
\]

The initial condition is

\[
\theta(\rho,0)=
\begin{cases}
T_0-T_f, & 0<\rho<r,\\
0, & r<\rho.
\end{cases}
\]

For example, room-temperature water at \(20^\circ\mathrm C\) is

\[
T_0 = 293.15\ \mathrm K,
\]

and a typical fridge at \(4^\circ\mathrm C\) is

\[
T_f = 277.15\ \mathrm K.
\]

The freezing point of pure water is

\[
273.15\ \mathrm K.
\]

---

## 3. Governing Heat Equation

In spherical symmetry, the heat equation is

\[
C(\rho)\frac{\partial \theta}{\partial t}
=
\frac{1}{\rho^2}
\frac{\partial}{\partial \rho}
\left[
k(\rho)\rho^2
\frac{\partial \theta}{\partial \rho}
\right],
\]

where

\[
C(\rho)=\rho_m(\rho)c_m(\rho)
\]

is volumetric heat capacity.

Here:

- \(k\) is thermal conductivity,
- \(\rho_m\) is material density,
- \(c_m\) is specific heat capacity,
- \(C=\rho_m c_m\) is volumetric heat capacity.

The spherical Laplacian is

\[
\nabla^2 T
=
\frac{1}{\rho^2}
\frac{\partial}{\partial \rho}
\left(
\rho^2\frac{\partial T}{\partial \rho}
\right).
\]

---

## 4. Material Properties

The coefficients are piecewise constant.

For the stagnant-air model,

\[
k(\rho)=
\begin{cases}
k_w, & 0<\rho<r,\\
k_g, & r<\rho<b,\\
k_a, & b<\rho<R,
\end{cases}
\]

and

\[
C(\rho)=
\begin{cases}
C_w=\rho_w c_w, & 0<\rho<r,\\
C_g=\rho_g c_g, & r<\rho<b,\\
C_a=\rho_a c_a, & b<\rho<R.
\end{cases}
\]

For the Robin convection model, the computational domain stops at the outer glass surface:

\[
0<\rho<b.
\]

Then

\[
k(\rho)=
\begin{cases}
k_w, & 0<\rho<r,\\
k_g, & r<\rho<b,
\end{cases}
\]

and

\[
C(\rho)=
\begin{cases}
C_w, & 0<\rho<r,\\
C_g, & r<\rho<b.
\end{cases}
\]

---

## 5. Boundary Conditions at the Center

At the center of the sphere, spherical symmetry implies zero radial heat flux:

\[
\frac{\partial \theta}{\partial \rho}(0,t)=0.
\]

Equivalently,

\[
\theta_\rho(0,t)=0.
\]

This condition applies to both models.

---

## 6. Interface Conditions

At the water-glass interface,

\[
\rho=r,
\]

temperature is continuous:

\[
\theta_w(r,t)=\theta_g(r,t),
\]

and heat flux is continuous:

\[
k_w\frac{\partial \theta_w}{\partial \rho}(r,t)
=
k_g\frac{\partial \theta_g}{\partial \rho}(r,t).
\]

Thus the temperature itself is continuous, but its slope may jump because \(k_w\) and \(k_g\) differ.

The physical condition is continuity of heat flux:

\[
k\theta_\rho
\]

not necessarily continuity of

\[
\theta_\rho.
\]

---

# Part I: Stagnant-Air Conduction Model

## 7. Stagnant-Air Physics

The stagnant-air model treats the air between the glass and fridge wall as a motionless conducting medium.

The air region is

\[
b<\rho<R.
\]

In this model, heat leaves the glass and diffuses through air by molecular conduction only.

The air satisfies

\[
C_a\frac{\partial \theta_a}{\partial t}
=
\frac{1}{\rho^2}
\frac{\partial}{\partial \rho}
\left[
k_a\rho^2
\frac{\partial \theta_a}{\partial \rho}
\right].
\]

The outer wall is fixed at the fridge temperature:

\[
T(R,t)=T_f.
\]

In excess-temperature form,

\[
\theta(R,t)=0.
\]

At the glass-air interface,

\[
\rho=b,
\]

we impose temperature continuity:

\[
\theta_g(b,t)=\theta_a(b,t),
\]

and heat-flux continuity:

\[
k_g\theta_g'(b,t)=k_a\theta_a'(b,t).
\]

---

## 8. Why Stagnant Air Cools Too Slowly

Air has very low thermal conductivity:

\[
k_a \approx 0.026\ \mathrm{W/(m\,K)}.
\]

If the model forces heat to conduct through centimeters of stagnant air, it creates a large thermal resistance.

The spherical conduction resistance of the air shell is

\[
\mathcal R_a
=
\frac{1}{4\pi k_a}
\left(
\frac{1}{b}
-
\frac{1}{R}
\right).
\]

For typical numbers, this can be tens of \(\mathrm{K/W}\), producing unrealistically long cooling times.

Therefore, the stagnant-air model is useful as a limiting case, but it is not usually a good model of a real fridge.

---

# Part II: Robin Convection Model

## 9. Convective Physics

In a real fridge, air is not a static spherical shell. It circulates due to:

- natural convection,
- fan-driven circulation,
- thermostat cycling,
- interaction with shelves and walls,
- thermal plumes near warm objects.

The true air-side equation would be an advection-diffusion problem:

\[
C_a
\left(
\frac{\partial T_a}{\partial t}
+
\mathbf u\cdot\nabla T_a
\right)
=
\nabla\cdot(k_a\nabla T_a),
\]

where \(\mathbf u\) is the air velocity field.

Solving that requires a fluid-flow model. The Robin approximation avoids that by replacing the external air domain with a boundary condition at the glass surface.

---

## 10. Robin Boundary Condition

At the outer glass surface,

\[
\rho=b,
\]

the heat conducted through the glass must equal the heat removed by convection into the fridge air.

Fourier conduction inside the glass gives surface heat flux

\[
q_{\text{cond}}
=
-k_g\frac{\partial T}{\partial \rho}(b,t).
\]

Newton cooling gives convective heat flux

\[
q_{\text{conv}}
=
h\left[T(b,t)-T_f\right],
\]

where \(h\) is the convective heat-transfer coefficient.

Equating the two fluxes gives

\[
-k_g\frac{\partial T}{\partial \rho}(b,t)
=
h\left[T(b,t)-T_f\right].
\]

In excess-temperature form,

\[
\boxed{
-k_g\frac{\partial \theta}{\partial \rho}(b,t)
=
h\theta(b,t).
}
\]

This is a **Robin boundary condition**, also called a mixed boundary condition.

---

## 11. Interpretation of \(h\)

The coefficient \(h\) has units

\[
\mathrm{W/(m^2\,K)}.
\]

It is not the same object as thermal conductivity.

Thermal conductivity has units

\[
\mathrm{W/(m\,K)}.
\]

The coefficient \(h\) summarizes the external air-side heat transfer:

\[
\text{near-surface conduction}
+
\text{air motion}
+
\text{boundary-layer thickness}
+
\text{geometry}
+
\text{fridge circulation}.
\]

A crude boundary-layer interpretation is

\[
h \sim \frac{k_a}{\delta_T},
\]

where \(\delta_T\) is an effective thermal boundary-layer thickness.

For air,

\[
k_a \approx 0.026\ \mathrm{W/(m\,K)}.
\]

So if

\[
h=5\ \mathrm{W/(m^2K)},
\]

then

\[
\delta_T\sim \frac{0.026}{5}\approx 5.2\ \mathrm{mm}.
\]

If

\[
h=10\ \mathrm{W/(m^2K)},
\]

then

\[
\delta_T\sim 2.6\ \mathrm{mm}.
\]

These are plausible effective boundary-layer scales for air around a small object.

---

## 12. Typical Values of \(h\)

For weak natural convection in air,

\[
h \sim 2\text{--}8\ \mathrm{W/(m^2K)}.
\]

For fan-assisted fridge circulation,

\[
h \sim 8\text{--}25\ \mathrm{W/(m^2K)}.
\]

A sensible baseline is

\[
h=5\ \mathrm{W/(m^2K)}
\]

for weak circulation, and

\[
h=8\text{--}10\ \mathrm{W/(m^2K)}
\]

for a moderately mixed fridge.

The script defaults to

\[
h=8\ \mathrm{W/(m^2K)}.
\]

It also includes a sensitivity comparison for

\[
h\in\{2,5,8,12,20\}.
\]

---

## 13. Limiting Cases of the Robin Boundary

The Robin condition is

\[
-k_g\theta_\rho(b,t)=h\theta(b,t).
\]

If

\[
h=0,
\]

then

\[
\theta_\rho(b,t)=0.
\]

The outer boundary is insulated. No heat leaves.

If

\[
h\to\infty,
\]

then the only way to keep the boundary flux finite is

\[
\theta(b,t)\to 0.
\]

So the Robin boundary approaches a fixed-temperature boundary:

\[
T(b,t)=T_f.
\]

Thus Robin convection interpolates between:

\[
\text{insulated boundary}
\]

and

\[
\text{fixed-temperature boundary}.
\]

---

## 14. Total Heat Flow at the Robin Boundary

The surface area at \(\rho=b\) is

\[
A_b=4\pi b^2.
\]

The total convective heat loss is

\[
Q(t)=hA_b\left[T(b,t)-T_f\right].
\]

In excess-temperature form,

\[
Q(t)=hA_b\theta(b,t).
\]

The same total heat flow through the glass surface is

\[
Q(t)
=
-4\pi b^2 k_g \theta_\rho(b,t).
\]

Equating these gives

\[
-4\pi b^2 k_g \theta_\rho(b,t)
=
4\pi b^2h\theta(b,t).
\]

Canceling the area gives the local Robin condition:

\[
-k_g\theta_\rho(b,t)=h\theta(b,t).
\]

---

# Part III: Eigenfunction Formulation

## 15. Continuous Eigenproblem

We seek separated solutions

\[
\theta(\rho,t)=\phi(\rho)e^{-\lambda t}.
\]

Substituting into the PDE gives

\[
-\frac{1}{\rho^2}
\frac{d}{d\rho}
\left[
k(\rho)\rho^2
\frac{d\phi}{d\rho}
\right]
=
\lambda C(\rho)\phi.
\]

This is a Sturm-Liouville eigenproblem with discontinuous coefficients.

The natural weighted inner product is

\[
\langle f,g\rangle_C
=
\int C(\rho)f(\rho)g(\rho)\rho^2\,d\rho.
\]

The eigenfunctions are orthogonal under this inner product.

---

## 16. Local Form of Eigenfunctions

Inside each homogeneous layer,

\[
k_i=\text{constant},
\qquad
C_i=\text{constant}.
\]

The eigenproblem reduces to

\[
\phi_i''+\frac{2}{\rho}\phi_i'
+
\mu_i^2\phi_i=0,
\]

where

\[
\mu_i^2=\frac{\lambda C_i}{k_i}
=
\frac{\lambda}{\alpha_i},
\]

and

\[
\alpha_i=\frac{k_i}{C_i}
\]

is thermal diffusivity.

The local radial eigenfunctions are

\[
\phi_i(\rho)
=
\frac{A_i\sin(\mu_i\rho)+B_i\cos(\mu_i\rho)}{\rho}.
\]

At the origin, regularity eliminates the singular term. Therefore in water,

\[
\phi_w(\rho)
=
A_w j_0(\mu_w\rho),
\]

where

\[
j_0(x)=\frac{\sin x}{x}.
\]

In glass,

\[
\phi_g(\rho)
=
\frac{A_g\sin(\mu_g\rho)+B_g\cos(\mu_g\rho)}{\rho}.
\]

In the stagnant-air model, air has the analogous form:

\[
\phi_a(\rho)
=
\frac{A_a\sin(\mu_a\rho)+B_a\cos(\mu_a\rho)}{\rho}.
\]

---

## 17. Eigenconditions for the Stagnant-Air Model

The stagnant-air eigenproblem is posed on

\[
0<\rho<R.
\]

Boundary conditions:

\[
\phi'(0)=0,
\]

\[
\phi(R)=0.
\]

Interface conditions at \(\rho=r\):

\[
\phi_w(r)=\phi_g(r),
\]

\[
k_w\phi_w'(r)=k_g\phi_g'(r).
\]

Interface conditions at \(\rho=b\):

\[
\phi_g(b)=\phi_a(b),
\]

\[
k_g\phi_g'(b)=k_a\phi_a'(b).
\]

The eigenvalues \(\lambda_n\) are the values for which all these conditions can be satisfied simultaneously.

---

## 18. Eigenconditions for the Robin Model

The Robin eigenproblem is posed on

\[
0<\rho<b.
\]

Boundary condition at the center:

\[
\phi'(0)=0.
\]

Interface conditions at \(\rho=r\):

\[
\phi_w(r)=\phi_g(r),
\]

\[
k_w\phi_w'(r)=k_g\phi_g'(r).
\]

Outer convective boundary condition at \(\rho=b\):

\[
-k_g\phi_g'(b)=h\phi_g(b).
\]

The Robin condition changes the eigenvalues and eigenfunctions. It does not merely change the conductivity. It changes the mathematical boundary condition representing the external air-side physics.

---

## 19. Eigenfunction Expansion

For either model, the solution has the modal form

\[
\theta(\rho,t)
=
\sum_{n=1}^{\infty}
c_n\phi_n(\rho)e^{-\lambda_n t}.
\]

The initial condition is

\[
\theta(\rho,0)
=
(T_0-T_f)\mathbf 1_{[0,r]}(\rho).
\]

The projection coefficients are

\[
c_n
=
\frac{
\int C(\rho)\theta(\rho,0)\phi_n(\rho)\rho^2\,d\rho
}{
\int C(\rho)\phi_n^2(\rho)\rho^2\,d\rho
}.
\]

Because the initial excess temperature is nonzero only in the water,

\[
c_n=
\frac{
C_w(T_0-T_f)
\int_0^r \phi_{wn}(\rho)\rho^2\,d\rho
}{
\int C(\rho)\phi_n^2(\rho)\rho^2\,d\rho
}.
\]

---

# Part IV: Finite-Volume Discretization

## 20. Discrete Heat Equation

The Octave code uses a finite-volume discretization.

For a spherical shell cell with edges

\[
\rho_{j-1/2},
\qquad
\rho_{j+1/2},
\]

the cell volume is

\[
V_j
=
\frac{4\pi}{3}
\left(
\rho_{j+1/2}^3-\rho_{j-1/2}^3
\right).
\]

The heat capacity of cell \(j\) is

\[
M_j=C_jV_j.
\]

The radial face area is

\[
A_{j+1/2}=4\pi\rho_{j+1/2}^2.
\]

The conductance between two neighboring cells is

\[
G_{j+1/2}
=
\frac{1}{
\dfrac{\Delta \rho_j/2}{k_j A_{j+1/2}}
+
\dfrac{\Delta \rho_{j+1}/2}{k_{j+1} A_{j+1/2}}
}.
\]

This expression is a thermal-resistance formula. It handles discontinuous conductivity correctly.

The semi-discrete system is

\[
M\dot{\theta}=-K\theta.
\]

---

## 21. Discrete Generalized Eigenproblem

The code solves

\[
K v_n=\lambda_n M v_n.
\]

The eigenvectors are normalized so that

\[
v_n^\top M v_n=1.
\]

Then the modal coefficients are

\[
c_n=v_n^\top M\theta(0).
\]

The discrete modal solution is

\[
\theta(t)
=
\sum_{n=1}^{N_m}
c_nv_ne^{-\lambda_nt}.
\]

Equivalently,

\[
\theta(t)
=
V\left(c\odot e^{-\lambda t}\right),
\]

where \(\odot\) denotes elementwise multiplication.

---

## 22. Finite-Volume Robin Boundary

At the outer glass surface,

\[
\rho=b,
\]

the Robin boundary is

\[
-k_g\theta_\rho(b,t)=h\theta(b,t).
\]

In the finite-volume code, the last glass cell is connected to the bulk fridge temperature \(\theta=0\) through two thermal resistances in series:

1. conduction from the last cell center to the glass surface,
2. convection from the glass surface to the bulk fridge air.

The surface area is

\[
A_b=4\pi b^2.
\]

The half-cell conduction resistance is

\[
\mathcal R_{\text{cond}}
=
\frac{\Delta \rho_N/2}{k_gA_b}.
\]

The convection resistance is

\[
\mathcal R_{\text{conv}}
=
\frac{1}{hA_b}.
\]

So the Robin boundary conductance is

\[
G_{\text{Robin}}
=
\frac{1}{
\mathcal R_{\text{cond}}+\mathcal R_{\text{conv}}
}.
\]

That is,

\[
G_{\text{Robin}}
=
\frac{1}{
\dfrac{\Delta \rho_N/2}{k_gA_b}
+
\dfrac{1}{hA_b}
}.
\]

The code adds this conductance to the last diagonal entry of the stiffness matrix:

```octave
K(N,N) = K(N,N) + Grobin;
