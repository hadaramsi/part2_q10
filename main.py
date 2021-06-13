# https://github.com/hadaramsi/part2_q10.git
import sympy
import sympy as sp
from sympy import ln
from sympy.utilities.lambdify import lambdify
# from numpy import log as ln
import math


def calcDerivative(func):
    """
    :param func: original function
    :return: None
    """
    x = sp.symbols('x')
    f_prime = func.diff(x)
    return f_prime


def simpson(f, startPoint, endPoint, parts):
    """
    :param f: original function
    :param startPoint: start of range
    :param endPoint: end of range
    :param parts: amount of segments
    :return: approximate area of the integral
    """
    if parts % 2 == 1:
        print("Amount of parts must be even")
        return None
    x = sp.symbols('x')
    func = lambdify(x, f)
    gap = abs(endPoint - startPoint) / parts
    appr = func(startPoint)
    for i in range(1, parts):
        if i % 2 == 0:
            appr += 2 * func((i * gap) + startPoint)
        else:
            appr += 4 * func((i * gap) + startPoint)
    appr += func(endPoint)
    appr *= 1 / 3 * gap
    return appr


def rombergMethod(f, a, b, end, epsilon):
    """
    :param f: Original function
    :param a: start of the range
    :param b: end of the range
    :param end: limit of iteration
    :param epsilon: allowed error
    :return: The area in the range
    """
    results = [[0 for i in range(end + 1)] for j in range(end + 1)]
    for k in range(0, end):
        res = trapezoidMethod(f, a, b, 2 ** k)
        results[k][1] = res
    for j in range(2, end + 1):
        for k in range(2, end + 1):
            results[k][j] = results[k][j - 1] + (
                    (1 / ((4 ** (j - 1)) - 1)) * (results[k][j - 1] - results[k - 1][j - 1]))
            if abs(results[k][j] - results[k - 1][j]) < epsilon:
                return results[k][j]


def trapezoidMethod(f, a, b, n):
    """
    :param f: Original function
    :param a: start of the range
    :param b: end of the range
    :param n: the number of the segments
    :return: The area in the range
    """
    x = sp.symbols('x')
    f = lambdify(x, f)
    h = (b - a) / n
    sum = 0
    while a < b:
        sum += 0.5 * ((a + h) - a) * (f(a) + f(a + h))
        a += h
    return sum


def rangeDivision(polinom, start_point, end_point, epsilon, function):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The excepted error
    :param function: The function that needs to be activated
    :return: None
    """
    results = []
    sPoint = start_point
    ePoint = end_point
    flag = False
    temp = start_point + 0.1
    while temp <= end_point:  # while we dont reach to the end point of the range
        print("Range: [" + str(start_point) + ", " + str(temp) + "]")
        res, iter = function(polinom, start_point, temp,
                             epsilon)  # activates the requested function with the original function
        if iter is not None:  # if the return iteration value is not None
            if (res > -epsilon) and (res < epsilon):  # check if the result is very close to 0
                res = 0
            print("The root is " + str(res) + "\nNumber of iteration is: " + str(iter))
            results.append(res)
            flag = True
        else:
            print("* There is no Intersection root in range ")
        der = calcDerivative(polinom)  # calculate the derivative
        x = sp.symbols('x')
        res, iter = function(der, start_point, temp,
                             epsilon)  # activates the requested function with the derivative function.
        if iter is not None:  # if the return iteration value is not None
            if (res > -epsilon) and (res < epsilon):  # check if the result is very close to 0
                res = 0
            if (lambdify(x, polinom)(res) > - epsilon) and (
                    lambdify(x, polinom)(res) < epsilon):  # only if the res is root
                print("The root is " + str(res) + "\nNumber of iteration is: " + str(iter))
                results.append(res)
                flag = True
            else:
                print("This point is a root of the derivative but is not a root of the function..\n")
        else:
            print("* There is no touching root in range\n")
        start_point = temp  # increase the start point of the range
        temp += 0.1  # increase the end point of the range
    if flag is False:
        print("There is no root in the range " + "[" + str(sPoint) + ", " + str(ePoint) + "]")
    return results


def NewtonRaphson(polinom, start_point, end_point, epsilon):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The excepted error
    :return: None
    """
    return rangeDivision(polinom, start_point, end_point, epsilon, calcByNewtonRaphson)


def calcByNewtonRaphson(pol, startPoint, endPoint, epsilon):
    """
    :param pol: Original function
    :param startPoint: int value, the start point of the range
    :param endPoint: int value, the end point of the range
    :param epsilon: The excepted error
    :return: The result and the number of the iteration for getting that result
    """
    Xr = (startPoint + endPoint) / 2  # middle of the range
    iteration = 1
    der = calcDerivative(pol)  # calculate the derivative
    x = sp.symbols('x')
    pol = lambdify(x, pol)
    der = lambdify(x, der)
    if pol(startPoint) * pol(endPoint) > 0:  # check if the is change in the sign in the function
        return None, None
    print("==== Iterations ====")
    while iteration < 100:
        res1 = pol(Xr)
        res2 = der(Xr)
        print(
            "Iteration number: " + str(iteration) + ", Xr = " + str(Xr) + ", f(x) = " + str(res1) + ", f`(x) = " + str(
                res2))
        iteration += 1
        Xnext = Xr - (res1 / res2)
        if abs(Xnext - Xr) < epsilon:
            print("Iteration number: " + str(iteration) + ", Xr = " + str(Xr) + ", f(x) = " + str(
                res1) + ", f`(x) = " + str(res2))
            print("==== End of iterations ====\n")
            return Xnext, iteration
        Xr = Xnext
    print("The system does not converge... :(")
    return None, None


def secant_method(polinom, start_point, end_point, epsilon):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The excepted error
    :return: None
    """
    return rangeDivision(polinom, start_point, end_point, epsilon, calcBySecant)


def calcBySecant(polinom, start_point, end_point, epsilon):
    """
    :param polinom: Original function
    :param start_point: int value, the start point of the range
    :param end_point: int value, the end point of the range
    :param epsilon: The excepted error
    :return: The result and the number of the iteration for getting that result
    """
    Xr = start_point
    Xnext = end_point
    iteration = 1
    x = sp.symbols('x')
    polinom = lambdify(x, polinom)
    if polinom(start_point) * polinom(end_point) > 0:  # check if the is change in the sign in the function
        return None, None
    print("==== Iterations ====")
    while iteration < 100:
        res1 = polinom(Xr)
        res2 = polinom(Xnext)
        print(
            "Iteration number: " + str(iteration) + ", Xr = " + str(Xr) + ", f(x) = " + str(res1) + ", f`(x) = " + str(
                res2))
        iteration += 1
        temp = Xnext
        Xnext = ((Xr * res2) - (Xnext * res1)) / (res2 - res1)
        Xr = temp
        if abs(Xnext - Xr) < epsilon:
            print("Iteration number: " + str(iteration) + ", Xr = " + str(Xr) + ", f(x) = " + str(
                res1) + ", f`(x) = " + str(res2))
            print("==== End of iterations ====\n")
            return Xnext, iteration
    print("The system does not converge... :(")
    return None, None


def checkDiffer(l, d, epsilon):
    print("check the difference between the methods:")
    flag = True
    for _ in range(len(l)):
        print("Root " + str() + ":\nSecant: " + str(l[_]) + ", Newton Raphson: " + str(d[_]))
        if abs(l[_] - d[_]) > epsilon:
            flag = False
            print("The difference is bigger than the epsilon for some of the roots")
            return
    print("The difference is smaller than the epsilon for all the roots")


def driver():
    """
    the main program
    :return: print the results
    """
    x = sp.symbols('x')
    f = sympy.ln(x ** 2)
    # f = (x * math.exp(1) ** (-x) + sympy.ln(x ** 2)) * (2 * x ** 3 + 2 * x ** 2 - 3 * x - 5)
    startRange = 0
    endRange = 1.5
    epsilon = math.pow(10, -4)
    print("Newton Raphson method")
    NewtonRaphson(f, startRange, endRange, epsilon)


driver()
