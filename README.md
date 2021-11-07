# Pràctica 2: Crank-Nicolson method for diffusion equation

This repository contains the program that implements the Crank-Nicolson method to solve the equation, given a domain 
```
du/dt - ∆u = f(t,x,y) for (x,y) in 
u = g(t,x,y)          for (x,y) in ∂
```
with _f_ and _g_ given in each exercise.

## Repository organisation
This repository consists on three library files of extension **.c* with ```main.c``` file which implements all the functions. The file ```methods.h``` contains the important stuff, that is, the Crank-Nicolson method and other methods to solve the problem. Then we have the ```auxFunctions.h``` file which contains auxiliar *C* functions that help make the code more readable. Finally we have the ```functions.h``` file which contains just the functions of the exercises.

As there are two exercises which demand the same but with different functions, I created a variable named ```option``` which is an integer that tells my program which exercise must be executed. This variable is asked for by terminal at the beggining of the execution. I will explain a little bit about each library file.

### File ```methods.h```
This file we have the main function, which is Crank-Nicolson. I followed what we were told in class and also some random information on the web such as [Wikipedia](https://en.wikipedia.org/wiki/Crank–Nicolson_method) or the given bibliography of the course.

The aim is to discretize the given space (which, in this case, is the 2-dimensional 1-square _[0,1]x[0,1]_) in a relatively large number, which in the program I usually call ```dim``` or ```dimension```. Then we create the system of equations which is obtained by the iteration of the Crank-Nicolson. This system of equations must be solved using another numerical method such as Jacobi or Gauss-Seidel. As in the previous exercise we saw that Conjugate Gradient was largely the most efficient method for this purpose, I rescued the function from this first practice. The code executes like this:

1. First we store the matrices and save the boundary values, given by _g_ function, with the function called ```boundaries```.
2. Then we use Crank-Nicolson to create the system of equations with the function called ```crank_nicolson```.
3. We solve the equation system with ```conjugate_gradinet```. Steps 1 and 2 are executed inside this function indeed.

### File ```auxFunctions.h```
This is a boring file which contains helper functions. Most of them have been rescued from practice one, so I will not explain them because they are already known.

### File ```main.c```
This is the main file, where everything is executed. I created an intermediate function called ```main_process``` which implements the algorithm, so in ```main()``` function we just have to do the boring work (declare variables and initialize them, store matrices, etc.) and then call the ```main_process```. In the ```main_process``` function I loop over the time interval given by each exercise, which has been discretized in _N_ steps. This _N_ variable is asked for by terminal at the beggining, just like _dim_ and the tolerance. Then we iterate over all these times and execute the ```conjugate_gradinent``` each iteration, so we can see the progress it makes and how it changes through time. The results can be printed out, but I prefered to print them out in an external file, so that it is easier to me to give the results as a solution. 