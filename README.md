# C++ Rasteriser and Raytracer using SDL

```RedNoise``` contains the original implementation for coursework submission.\
```2.0``` contains a revised and re-structured version to assist in being a teaching assistant on the unit.

## Building

Both versions can be compiled and run with ```make``` or ```cmake```. Each Makefile contains a 'speedy' rule for optimal performance:
```
cd 2.0/
make speedy
```
```
cd RedNoise/
make speedy
```


## 2.0 Controls

``WASD`` to move the camera forward-backward / left-right.\
``EQ`` to move the camera up-down.\
``Arrow Keys`` to rotate the camera.\
``LCTRL`` to orbit the camera around the model.\
``SPACE`` to turn the camera towards the model.\
``LALT`` to reset the camera's rotation and position.\
``123`` to switch between render modes.\
``ESC`` to exit.
