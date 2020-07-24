# Truss
This is a programm for solving every 2D truss problem. This solving method is based on deformation of the elements of the truss.

This is a project Saleh Faghfoorian, mechanical engineering undergraduate student at Sharif University of Technology.

## Input Nodes File's Format
you should enter your nodes information in this format :
```
1st Column : number of the node
2nd Column : x component of position the node (in the cartesian coordinate system)
3rd Column : y component of position the node (in the cartesian coordinate system)
4th Column : x component of Force applied on the node
5th Column : y component of Force applied on the node
6th Column : x component of constraint of the node (0 for no restriction , 1 for restricion in displacement)
7th Column : y component of constraint of the node (0 for no restriction , 1 for restricion in displacement)
```
## Input Elements File's Format
you should enter your elements information in this format :
```
1st Column : number of the element
2nd Column : number of node on the buttom of the element
3rd Column : number of node on the top of the element
4th Column : Area of the element
5th Column : Young's modulus of the element
6th Column : Yield Strength of the element
```
* Now you just need to run the code
* Job's done! :sunglasses:
* Now you have a report which tells you the forces, displacements, stresses, and the factor of safety.
* Also you have 2 diagrams.
* 1st diagram is the truss before and after deformation. ( Before : "gray" , After : "black" )
* 2nd diagram is the force of each element.

## Developers
* **Saleh Faghfoorian** [Saleh Faghfoorian](https://github.com/saleh-faghfoorian)
