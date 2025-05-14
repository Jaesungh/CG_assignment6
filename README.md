# CG_assignment6

### Introduction
This project implements a software-based triangle rasterizer that renders a unit sphere using a custom rendering pipeline in software. The pipeline includes :

  - Model Transformation : scales the sphere to radius 2 and translates it to z = -7
  - Camera transformation : positions the eye at (0,0,0) looking down -z
  - Perspective projection : maps 3D eye-space into normalized device coordinates
  - Viewport transformation : maps NDC to a 512x512 image
  - Triangle rasterization : uses barycentric coordinates, depth testing, and per-pixel shading algorithms

It finally writes out three PPM images (output_flat.ppm, output_gouraud.ppm, output_phong.ppm) showing Flat, Gouraud, and Phong shading respectively. 

### Result 

Below are the rendered spheres for each shading mode : Flat, Gouraud, and Phong. 

Flat Shading : 
![output_flat](https://github.com/user-attachments/assets/308acf4f-4894-4d7e-9afe-f5e082d926d1)

Gouraud Shading : 
![output_gouraud](https://github.com/user-attachments/assets/82c2373f-ec44-41d0-b566-1eeb3aadfd5c)

Phong Shading :
![output_phong](https://github.com/user-attachments/assets/eedd9a7f-ba4a-4478-90b5-a1b5c9b8c449)

### Compilation & Run Instructions

Requirements : Visual Studio 2022 (Windows)

Build Steps :
1. Open the provided .sln file in Visual Studio 2022.
2. Set platform to x64.
3. Build the solution (Ctrl+Shift+B)
4. Run without debugging (Ctrl+F5)

### Description
  - create_scene() : Constructs a triangular mesh of a unit sphere using spherical coordinates.
  - modelTransform() : Applies uniform sclae and translation to position the sphere in front of the camera.
  - cameraTransform() : Leaves eye-space coordinates unchanged (eye at origin)
  - perspectiveProject() : Applies a canoncial perspective projection.
  - drawTriangle() : Rasterizes each triangle via barycentric interpolation, depth testing, and chooses shading per pixel based on the selected mode (Flat, Gouraud, Phong)

### Files 
  - rasterized_wshading.cpp : Main rendering logic
  - output_flat.ppm : Flat shaded result
  - output_gouraud.ppm : Gouraud-shaded result
  - output_phong.ppm : Phong shaded result

