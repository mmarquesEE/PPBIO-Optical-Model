\documentclass{article}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{geometry}
\geometry{a4paper, margin=1in}

\begin{document}

\section*{Mathematical Explanation of MATLAB Code Block}

This section details the mathematical operations performed by the provided MATLAB code, which simulates the propagation of a reflected light beam and visualizes its intensity on a plane orthogonal to its propagation direction.

\subsection*{1. Reflected Angular Spectrum ($A_r$)}
The angular spectrum of the reflected field, $A_r(f_x)$, is calculated by multiplying the incident angular spectrum, $A_i(f_x)$, with the frequency-dependent reflection coefficient for TM polarization, $r_p(f_x)$.
$$ A_r(f_x) = r_p(f_x) \cdot A_i(f_x) $$
Here, $f_x$ represents the spatial frequency components corresponding to `fs1` in the code. $A_i(f_x)$ corresponds to `Ai1` and $r_p(f_x)$ corresponds to `rp1`.

\subsection*{2. Propagation Parameters}
Parameters for propagation in the initial medium (refractive index $n_1$, wavelength $\lambda_1 = \lambda_0/n_1$) are defined.
\begin{itemize}
    \item Wavenumber in medium 1: $k = \frac{2\pi}{\lambda_1}$
    \item Propagation distance: $L$
    \item Angle of incidence/reflection: $\theta$
    \item Propagation vector components in global coordinates: $(\Delta x, 0, \Delta z)$
    $$ \Delta z = -L \cos(\theta) $$
    $$ \Delta x = L \sin(\theta) $$
\end{itemize}

\subsection*{3. Angular Spectrum Propagation}
The propagation is performed in the frequency domain using the Angular Spectrum Method (ASM).
\begin{itemize}
    \item The ASM propagation kernel for distance $\Delta z$:
    $$ H_{prop}(f_x) = \exp\left( j k_z \Delta z \right) = \exp\left( j k \Delta z \sqrt{1 - (\lambda_1 f_x)^2} \right) $$
    where $k_z = k \sqrt{1 - (\lambda_1 f_x)^2}$ is the z-component of the wavevector and $j = \sqrt{-1}$.
    \item The phase shift operator for lateral displacement $\Delta x$:
    $$ H_{shift}(f_x) = \exp\left( j k_x \Delta x \right) = \exp\left( j (2\pi f_x) \Delta x \right) $$
    where $k_x = 2\pi f_x$.
    \item Evanescent Wave Filter: Only propagating waves, where $k_z$ is real, are kept. This requires $1 - (\lambda_1 f_x)^2 \ge 0$, or $|f_x| \le 1/\lambda_1$. A mask function, $\text{Mask}(f_x)$, is applied:
    $$ \text{Mask}(f_x) = \begin{cases} 1 & \text{if } |f_x| \le 1/\lambda_1 \\ 0 & \text{if } |f_x| > 1/\lambda_1 \end{cases} $$
    This corresponds to `prop_mask` in the code.
    \item The propagated angular spectrum, $A_{r,prop}(f_x)$, is the product of these factors:
    $$ A_{r,prop}(f_x) = A_r(f_x) \cdot H_{prop}(f_x) \cdot H_{shift}(f_x) \cdot \text{Mask}(f_x) $$
    This corresponds to `Ar1_propagated`.
\end{itemize}

\subsection*{4. Field Reconstruction and Intensity}
The complex electric field amplitude in real space, $U_{far}(x')$, is obtained by the inverse Fourier transform of the propagated angular spectrum. The spatial coordinates $x'$ correspond to `x_cam1`.
$$ U_{far}(x') = \mathcal{F}^{-1}\{A_{r,prop}(f_x)\}(x') $$
The optical intensity $I(x')$ is the squared magnitude of the complex field:
$$ I(x') = |U_{far}(x')|^2 $$
This corresponds to `Intensity`.

\subsection*{5. 2D Intensity from 1D Data}
The code assumes the beam is invariant along the Y-direction to create a 2D intensity map from the 1D result $I(x')$.
$$ I_{2D}(x', y') = I(x') $$
The `meshgrid` and `repmat` functions generate the 2D grid and replicate the intensity data.

\subsection*{6. Orthogonal Plane Definition}
A plane orthogonal to the reflected beam direction is defined for visualization.
\begin{itemize}
    \item Reflected beam direction vector (unit vector): $\vec{n} = (\sin\theta_r, 0, -\cos\theta_r)$, where $\theta_r = \theta$.
    \item Basis vectors spanning the orthogonal plane: $\vec{v}_1 = (0, 1, 0)$ and $\vec{v}_2 = \frac{\vec{n} \times \vec{v}_1}{\|\vec{n} \times \vec{v}_1\|}$.
    \item Origin (reference point) of the plane in global coordinates: $\vec{r}_0 = (0, 0, -200 \times 10^{-6})$.
    \item Parametric equation for points $\vec{R}_{plane}(u, v)$ on the plane, using local coordinates $(u, v)$:
    $$ \vec{R}_{plane}(u, v) = \vec{r}_0 + u \vec{v}_2 + v \vec{v}_1 $$
    The components $(X_{plane}, Y_{plane}, Z_{plane})$ of $\vec{R}_{plane}$ give the global coordinates of the grid points on the orthogonal plane.
\end{itemize}

\subsection*{7. Intensity Interpolation}
The calculated intensity $I_{2D}$, defined on the grid $(X_{orig}, Y_{orig}) = (x'_{cam} + \Delta x, x'_{cam})$, is interpolated onto the coordinates $(X_{plane}, Y_{plane})$ of the orthogonal plane.
$$ I_{rot}(u,v) = \text{Interpolate}_{2D}[I_{2D}(X_{orig}, Y_{orig}) \text{ at } (X_{plane}(u,v), Y_{plane}(u,v))] $$
The function `interp2` with 'linear' interpolation is used. The result corresponds to `RotatedIntensity`. Finally, the intensity is normalized:
$$ I_{norm}(u,v) = \frac{I_{rot}(u,v)}{\max_{(u',v')}(I_{rot}(u',v'))} $$

\subsection*{8. Plotting}
The normalized, interpolated intensity $I_{norm}(u,v)$ is visualized as a colored surface plot (`surf`) located at the global coordinates $(X_{plane}, Y_{plane}, Z_{plane})$.

\end{document}