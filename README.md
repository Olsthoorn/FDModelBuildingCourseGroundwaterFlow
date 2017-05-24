# Finite Difference Modeling course for groundwater flow
Graduate coarse given at TUDelft until end 2014

The course shows how a finite difference groundwater model can be built up from base in Python and how it is used.

A 3D steady state finite difference model (FDM) is built up from ground and later extended to transient.
The models are implemented as Python functions that can be called with a simple signature. A number of examples show how they can be used and also verifies their correctness by comparison with analytical solutions.

The FDM models work for 3D as well as for axially symmetric cases. Axially symmetric cases are simulated in a grid row (x-z). One can simulationed a large number of axially-symmetric cases simultaneously, one in each row, like cross sections.

The FDM grid/mesh is conveniently stored and handled by the Grid class. Only grids of block form are possible, yet their size is practially unlimited and can be adapted by making parts inactive like it is done in MODFLOW. The more cell that are inactive, the faster the model will run. The switch axial=True in the call of the grid suffices to treat model rows as separated axially symmetric models in the xz plane and the x-coordinate as radial distance to the center.

In the original version, uploaded December 2016, 3D arrays were [Ny, Nx, Nz] oriented as it is the case in Matlab. In the renewed upload of May 24, 2017, the arrays are [Nx, Ny, Nx] oriented, as they are in the FDM MODFLOW of the United States Geological Survay (USGS). The reason is that this interpretation of the 3D arrays representing FDM grid data is smarteds in Python. It implies that if the head Phi is an [Nz, Ny, Nx] 3D array, Phi[i] (=Phi[i, : :]) is the a 2D array of the ith-1 layer and Phi[i][j] (=Phi[i, j, :]) is a 1D array representing the jth-1 row of that layer and, of course, Phi[i][j][k] (=Phi[i, j, k] is zero-dimension kth cell of that row, equal to Phi[i, j, k]. Likewise in transient modeling we interprete the resulting 4D arrays as [Nt, Nz, Ny, Nx] so that the Phi[m] of the transient 4D array is the entire 3D array of stress period Phi[m-1] and Phi[m, i] the head in layer i-1 of that stress period and so on. This is elegant and intuitive. Morover, for cross sections along the x-axis we can select an entire cross section of heads of a steady-state model like so: Phi[:, j, :] which is a 2D array of size [Nz, Nx] that can immediately be contoured withou any need to transpose it. A cross section along column k-1 is selected as Phi[:, :, k]. And selecting the heads of stress period m-1 of a transient model is Phi[m, :, :, k]. Again, this is a 2D array of size [Nz, Ny] that can immediately be contoured.

Streamlines are derived and implemented as the stream function that is contoured, which is useful in divergence-free 3D flow. In practice this is in cross sections, where the stream function and the stream lines are extremely useful.

Additionally a 3D particle tracker is made and used. It allows to readily track hundreds of thousands of particles simultaneously. These particles can be shown in 3D or 2D.

Simple, complex, small and large models can readily be constructed. As sparse matrices are used, the only limit is the size of the computer's memory.

The flow boundary conditions are limited to fixed-heads and prescribed flows and cells may be labeled `inactive`. However the theory of how to implement special packages like MODFLOW's DRAIN, RIVER and GHB is explained. Implemtation is straightforward, but was not done to keep the signature of the FDM functions lean. The more cells with a head fixed, the faster the model.

For the same reasons, cells falling (partially) dry under declining water tables have not been implemented. However, variably filled cells can be tackled in a simple loop around the solver.

I consider the course as one to provide hands-on insight in FDM model construction and FDM groundwater modeling. For professional regional modeling using special packages, is probably better off using the USGS's programs MODFLOW and MODPATH, perhaps using the USGS's Python interface called FlowPy. The current course may be considered as an in depth background introduction for those who want to use FlowPy.

The idea to use Python is my conviction that students should nowadays learn two "things"  Python and QGIS. Not only can one compute
and do practically everything that is possible to code and visualize, but Python also allows one to cooperate worldwide across any boundaries, physical or political. Python and QGIS can be considered for engineers and scientists what is English for worldwide communcation. Python and QGIS are free, so nobody can take it from you, no matter where you are and who you are. It's up to you and only you to learn to use these tools for any computational and visualization task. Doing so will make you a master in it, promise you a liflong benefit. The posibilities provided to each individuals
by these tools are simply unprecedented.

Theo Olsthoorn, 30-Dec-2016
