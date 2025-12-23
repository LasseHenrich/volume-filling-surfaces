# Volume-Filling Surfaces (and Curves)
<div class="one-half column category" style="text-align: center;">
    <video class="u-max-full-width" style="padding-left: 10rem; padding-right: 10rem" autoplay loop muted>
      <source src="_readme-sources/volume_filling_curve.mp4" type="video/mp4">
    Your browser does not support the video tag.
    </video>
  <p>Fig. 1: Filling a volume with a curve</p>
</div>
  <div class="one-half column category" style="text-align: center;">
    <video class="u-max-full-width" style="padding-left: 10rem; padding-right: 10rem" autoplay loop muted>
    <source src="_readme-sources/volume_filling_surface.mp4" type="video/mp4">
    Your browser does not support the video tag.
    </video>
  <p>Fig. 2: Filling a volume with a surface</p>
</div>

We present a novel algorithm to generate closed surfaces that fill a given volume with a given radius. Our algorithm starts with some arbitrary surface and then minimizes an energy that is based on an implicit medial surface representation via standard gradient descent. As an intermediary step, we also developed an approximation for volume-filling curves. Our results, especially for the case of surfaces, prove effective in constructing volume-filling manifolds that look asthetically pleasing. While we do not conduct a thorough comparison to the current state of the art, we have reason to believe that our algorithm can generate comparable results in far less time.

A publication on _Volume-Filling Surfaces_ will follow.<br>
Note that our _Volume-Filling Curves_ are fairly well described in [this recently published paper](https://www.dgp.toronto.edu/projects/medial-sphere-preconditioning/), though their approach is seems polished. (Our project existed before this paper came out.)

## Contents
1. [Project overview and technical background](#project-overview-and-technical-background): Detailed report on what this project 
2. [Setup](#setup)
3. [Controls](#controls)


## Project overview and technical background
### Introduction
Given some two-dimensional manifold, a <i>surface-filling curve</i> is a curve that lies on the manifold and fills it uniformly.
In other words, every point on the surface has a given maximum distance to the curve, measured geodesically.
Surface-filling curves are primarily used for artistic creation and content generation [1-3],
but also have real-world use cases such as manufacturing [4, 5].<br><br>
Many attempts at an efficient algorithm for generating surface-filling curves have been published (ref. [6, sec. 2]). <i>Repulsive Curves</i> [2]
from Yu et al. was published in 2021 and presents an approach for minimizing a curve's repulsive energy under constraints. While this was originally aimed at solving
different problems in the realm of computational design (global shape optimization while avoiding self-intersections),
the resulting geometric flow could fairly easily be adapted to what they
call <i>curve packing</i>. Curve packing works by growing curves under self-repulsion and constrained to some subspace,
e.g. a surface or a volume, which we label as surface- or volume-filling in this paper respectively.
However, as surface- and volume-filling was merely a bi-product of
Yu et al.'s work, their algorithm was not optimized for that specific task. For the case of surface-filling curves,
Noma et al. therefore developed a much faster solution in 2024 [6],
using an implicit representation of the curve's medial axis to satisfy the surface-filling constraint. However, prior to our paper,
<i>Repulsive Curves</i> remained the state-of-the-art work for volume-filling curves, as Noma et al. have not applied their geometric flow to volumetric domains.<br><br>
Transitioning from curves to surfaces, the currently most important paper for volume-filling surfaces was again published by Yu et al. in 2021
and is called <i> Repulsive Surfaces</i> [7]. Here, Yu et al. once more research global shape optimization without self-intersections,
but this time for surfaces, achieving the then most efficient method for generating volume-filling surfaces as a bi-product.
In contrast to surface-filling curves, to the best of our knowledge, no other paper has been published that specifically tackles
volume-filling surfaces and surpasses the results of [7].<br><br>
In this paper, we do for volume-filling surfaces what Noma et al. did for surface-filling curves: Developing an efficient algorithm
specifically for this use case. We do this by simply applying their idea of using an implicit medial axis, but move from curves to surfaces and from
two-dimensional to three-dimensional manifolds as containing domains, obtaining a fast flow for generating volume-filling surfaces. As intermediate steps,
we first build <i>space-filling</i> curves and then introduce a volumetric boundary-constraint, which results in a faily good approximation
of <i>volume-filling curves</i>. Afterwards, transitioning from curves to surfaces is not difficult.

### Technical Background
What follows are the basic concepts needed to develop our approach for volume-filling curves and surfaces.
We start by defining the Frenet-Serret frame and then transition to medial axes and surfaces.

#### Frenet-Serret frame
<div style="
    display: flex;
    flex-direction: column;
    align-items: center;
">
<img style="width: 50%" src="_readme-sources/documentation/frenet_frame.png">
    <p>Fig. 3: Visualization of the Frenet-Serret frame. Tangents are red, normals are blue and bitangents are yellow.</p>
</div>

#### Medial axes of curves and surfaces
#### Calculating distances from the medial surface
### Method
#### Medial surface energy
#### Volumetric energy
#### Length energy
#### From curves to surfaces
### Implementation
#### Approximating the medial energy for curves
#### Computing distances to the medial surface
#### Computing distances to the volume's boundary
#### Evaluating $\alpha$ and $\beta$
### Results

## Setup
The project is a CMake project. It has, however, been developed with Visual Studio 2022, and there are some dependencies to that development environment and the build tools.
1. If you have vcpkg globally installed, proceed with the next step. If not: vcpkg is a library manager that we need to use OpenVDB, the library that generates the SDF for us. If you have not installed vcpkg at all, install it globally on your system. (You can also install it locally for the project, but I'm not sure how to use it then.)
2. Use vcpkg to install OpenVDB. It should work by just running `vcpkg install openvdb`, but you might need to play around with this. Afterwards, when running `vcpkg list` in the terminal, you should see _openvdb_ listed. Now, `find_package(OpenVDB, ...)` in CMakeLists.txt should work.
1. Before opening the project with Visual Studio, open the CMakePresets.json file and replace _PATH_TO_VCPKG_ with the path to your vcpkg installation.
1. Run Visual Studio as an __administrator__.
1. Open the _development_ folder (just as a folder, not a VS project).
1. Make sure that volume-filling-curves.exe is selected in the target dropdown. Also, you may want to switch to release mode.
![alt text](_readme-sources/setup/image-2.png)
1. Wait for the CMake generation to finish and then build the project.

You can run the project directly from within Visual Studio, but you have to specify a scene file. One way of doing that is to specify a launch.vs.json file, which you can do by right-clicking the CMakeLists.txt and opening the _Debug and Launch Settings_.

![alt text](_readme-sources/setup/image-1.png)

Afterwards, paste the following contents:
```json
{
  "version": "0.2.1",
  "defaults": {},
  "configurations": [
    {
      "type": "default",
      "project": "CMakeLists.txt",
      "projectTarget": "volume-filling-curves.exe",
      "name": "volume-filling-curves.exe",
      "args": [
        "PATH_TO_PROJECT/models/98480/scene_small.txt"
      ]
    }
  ]
}
```
Alternatively (e.g. if you don't see that option in the context-window), you can create a launch.vs.json file in the .vs folder. 

Note that while the entry file is named _volume-filling-curves.cpp_, it does also contain the code for volume-filling _surfaces_.

## Controls
The application is based on polyscope and the UI should be familiar. On the right, you can run the algorithm step by step or for a certain number of iterations. You can also click the space bar to run a single iteration. On the left, multiple visualizations can be controlled.

The scene file that you have to pass can contain a lot of different options. Refer to the many examples in this project (e.g. 72286, 98480, 68380) for guidance, but note that some scenes might not work anymore because of changes in accepted inputs and the algorithm itself. Also refer to _scene_file.cpp_ for more insights.