# Pràctica 2: Crank-Nicolson method for diffusion equation

This repository contains the program that implements the Crank-Nicolson method to solve the equation
```
du/dt - ∆u = f(t,x,y)
u = g(t,x,y)
```
with _f_ and _g_ given in each exercise.

## Repository organisation
This repository consists on three library files of extension **.c* with ```main.c``` file which implements all the functions. The file ```methods.h```contains the important stuff, that is, the Crank-Nicolson method and other methods to solve the problem. Then we have the ```auxFunctions.h```file which contains auxiliar *C* functions that help make the code more readable. Finally we have the ```functions.h```file which contains just the functions of the exercises.

As there are two exercises which demand the same but with different functions, I created a variable named ```option```which is an integer that tells my program which exercise must be executed. This variable is asked for by terminal at the beggining of the execution. I will explain a little bit about each library file.

### File ```methods.h```
This file we have the main function, which is Crank-Nicolson. I followed what we were told in class and also some random information on the web such as [Wikipedia][https://en.wikipedia.org/wiki/Crank–Nicolson_method] or the given bibliography of the course.
