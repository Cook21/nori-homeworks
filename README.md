Advanced Computer Graphics — Homeworks
## Dependencies
### 1. Install Packages
Ubuntu:
```
sudo apt install cmake-gui zlib1g-dev xorg-dev libglu1-mesa-dev
```
Frdora:
```
sudo dnf install cmake-gui zlib-devel mesa-libGLU-devel libXi-devel libXcursor-devel libXinerama-devel libXrandr-devel xorg-x11-server-devel
```
### 2. Download Submodules
Use the following command to download all 3rd-patry libs from github.
```
git submodule update --init --recursive
```
### 3. Optinal: Update the Starter Code
```
$ git remote add upstream https://github.com/cs440-epfl/nori-base-2021.git
$ git pull upstream master
```
## How To Build
```
$ cd path-to-nori
$ mkdir build
$ cd build
$ cmake-gui ..
```
In the pop up cmake window,click configure, then click generate, run `$ make -j ` afterward.

## How to Use
Nori can be run in the command line by specifying the path to an XML scene file:
```
$ ./nori path/to/scene.xml
```



======================================

Student name:

Sciper number:


## Build status

**Insert your build badge URL here**

## Homework results

| Homework   |  Links
| ---------: | ---------------------------------------------
| 1          | [report.html](results/homework-1/report.html)
| 2          | [report.html](results/homework-2/report.html)
| 3          | [report.html](results/homework-3/report.html)
| 4          | [report.html](results/homework-4/report.html)
| 5          | [report.html](results/homework-5/report.html)


## Featured result

Feel free to show off your best render here!
