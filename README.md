----------------------------------------------------------------
최적설계 과목 프로젝트에서 작성한 코드(5명 공동)
-----------------------------------------------------------------

# Optimization

> Lets solve optimization problem using this C codes. 

In this project, we firstly defined an engineering problem with objective function, variables, and constraints.
Secondly, We wrote them on a C program. We used Quasi Newton Method(BFGS) to find direction for optimization, and used "Equal interval section method" or "Golden Section method" or "Quadratic approximation method" to search initial interval and decide final interval reducing steps.

My contributions to this code are following. -> functions named "Equal interval section method", "Golden section method", "Quadratic approximation method" were fully coded by me.

> Matlab codes contains fmincon function. It is an usage of original function "fmincon" in Matlab library. 

> Fotran code using steepest descent method, conjugate gradient method, and BFGS method is also following. It specifies 3 direction searching method, and use only 1 method,"Golden section method" to decide optimization step size. These Fotran code was scripted by my coworker, Hong.

> To sum up, We could use this Fotran code to compare the performances from three different direction search method: Steepest Descent Method, Conjugate Gradient Method, and Quasi Newton Method(BFGS). And also we could use the C code to compare the performances from three different step size decision method: Equal Interval Method, Golden Section Method, Quadratic Approximation Method.
 Result is that in simple case with less iteration calculation, Conjugate Gradient Method could be more effective than BFGS. But mostly Quasi Newton Method is the best way within those three methods. Similarly, Golden Section Method could be more effective in some cases, but mostly Quadratic Approximation Method is the best way among those.



[We referred algorithms and main ideas by a textbook, [Jasbir_S._Arora]_Introduction_to_optimum_design].

