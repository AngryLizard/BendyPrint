# BendyPrint

This project uses [libIGL](https://github.com/libigl/libigl). To install just clone this project including submodules, or clone [libIGL](https://github.com/libigl/libigl) yourself into the project folder. To run, create a build/ folder and run cmake.

Bendyprint is split in 4 modes that can be switched at the top of the menu.

## Selection

Selection lets you select which vertices are fixed (black) and which have gravity applied to them (green).
Leftclick and drag on the mesh for adding, rightclick for removing. Selection mode (fixed vs gravity) can be changed with the combo-box at the bottom.
The selection doesn't reset on rebuilds, only when new meshes are loaded.

## Simulation

FEM simulation. The shear and bulk as well as the shell thickness can be changed for the mesh. Rebuild mesh to apply the changes.
Refresh rerenders the mesh.
Test FEM applies a finite difference test on all elements, results are printed into stdout.
Pick vertex displaces vertex 0 downwards.
Iterate applies 100 iterations of gradient descent on the FEM simulation.

## Wrap

Wrap shows the unwrapped version of the mesh. Fill can be turned on or off to show the transformed shell elements.

## Slice

White dots signify printhead jumps.
The G-Code generator works specifically for the Ultimaker 2 and uses a Cuda format. Might work on other printers but no guaranties!

Slice contains all the G-Code generator settings. Multiple copies can be printed, make sure they don't overlap the printing bed before printing. If "Slice via outlines" is disabled, triangles are sliced instead. Line stride changes the distance between fill-lines. Line union threshold changes the threshold between two fill-lines for when the printhead jumps.
