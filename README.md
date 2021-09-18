# thermaNet
A small modelling and simulation tool for thermal lumped models, using graphs
This small project comes from a course realised at PoliMi in Modelling and Simulation of Aerospace systems. The original scripts were realised in MATLAB, with the aim of simulating the thermal evolution of a CubeSat structure, in the context of the LUMIO mission. My role was to develop the thermal simulation, from which the following code is inspired.

The goal was to develop a simple tool for thermal simulations, using graphs to easily increase the number of lumps for a CubeSat structure.

The CubeSat graph script will be reimplemented in Python with time.

A small improvement of the method used is the computation of the matrix corresponding to the linear application corresponding to thermal conduction, making the eigenvalues of the problem readily accessible for thermal conduction problems. The same principle is applied to radiation, with work to be done on the computation of the eigenvalues. The principles used for this computation are summed up in the demo.pdf file.
The computation of such a matrix can become quite burdensome for large number of lumps and therefore may be used with care.

## Use 

These scripts are making use of  Numpy (general computations), NetworkX (graph library), Matplotlib (plotting) and Scipy (ODE integration).

To test the simulation simply launch the graph_sim_simple.py script, simulating the thermal evolution of an adiabatic beam with temperature boundary conditions set at each ends, with an initial parabolic thermal profile.

## To be done 

- CubeSat graphs
- Heat source/sink
- Convection
- GUI for graph management
- Orbit definition and Sun radiation
- Albedo
- Use of sparse matrices ?

## Finally

I hope it can be useful some day for someone, and if you wish to contribute, you are welcome !