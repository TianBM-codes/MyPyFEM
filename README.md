## MyPyFEM
用Python写的有限元程序,供大家学习交流

## Support Input File Type:
* Nastran *.bdf
* Abaqus *.inp
* Ansys *.cdb

## Support Element
* Truss
* Plane
* Beam
* MITC4 Shell
* DKT Shell
* DKQ Shell

## Support Analysis Type
* Static
* Nonlinear

## 代码约定
* 全局内容字母全部大写

## To Do List
- [ ] Dynamic Analysis
- [ ] Stress
- [ ] Compare With Commercial Software

## Reference Book
1. Javier Bonet 《Nonlinear Solid Mechanics for Finite Element Analysis: Statics》
2. Bath 《Finite Element Procedures》
3. DRJ Owen 《Computation Methods For Plasticity Theory and Application》
4. Eduardo 《Notes on Continuum Mechanics》

## Question
- [ ] 在弧长法计算的过程中，均布载荷为什么在外边计算
- [ ] 在均布载荷积分的时候不乘以面积吗？PressureElementLoadAndStiffness函数



## 待参考实现

* [Welcome to Numerical tours of Computational Mechanics using FEniCS — Numerical tours of continuum mechanics using FEniCS master documentation (comet-fenics.readthedocs.io)](https://comet-fenics.readthedocs.io/en/latest/index.html)

  