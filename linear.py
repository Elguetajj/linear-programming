import numpy as np
import sympy as sp
from sympy.solvers.inequalities import reduce_rational_inequalities

x,y = sp.symbols('x,y')

class LinearSolver:

    def __init__(self,objective,constraints,maxmin):
        '''inicializes an instance of the class Linear Solver.

            Args:
                objective (string): Objective function of the problem as a string in numpy format.
                constraints (list of strings) : Problem constraints as a list of strings in numpy format.
                maxmin (boolean): 0 if the problem is a maximization, 1 if it is a minimization. 

        '''
        self.__constraints: list = constraints
        self.__constraints_as_eqs: list = []
        self.__objective: str = objective
        self.__intercepts: list = []
        self.__posible_solutions = []
        self.__solution_points: dict = None
        self.__solution: int = None
        self.__maxmin: str = maxmin
        self.__binding_constraints: list = None
        self.__binding_slopes:list=[]
        self.__sensitivity:dict = {'x': None , 'y': None}

    def constraints_to_eqs(self):
        ''' Converts restrictions to equations equaled to 0 and saves them in the private field "__constraints_as_eqs" as an array of strings.
        '''
        for constraint in self.__constraints:
            if '<=' in constraint:
                self.__constraints_as_eqs.append(constraint.replace('<=','-'))
            elif '>=' in constraint:
                self.__constraints_as_eqs.append(constraint.replace('>=','-'))
            elif '==' in constraint:
                self.__constraints_as_eqs.append(constraint.replace('=','-'))
            elif '<' in constraint:
                self.__constraints_as_eqs.append(constraint.replace('<','-'))
            elif '>' in constraint:
                self.__constraints_as_eqs.append(constraint.replace('>','-'))  
        print(self.__constraints_as_eqs)
    
    def find_intercepts(self):
        ''' Uses "__constraints_as_eqs" to find the intercepts of the problem and saves them to "__intercepts"
        '''
        for i in range(0,len(self.__constraints_as_eqs)-1):
            for j in range(i+1,len(self.__constraints_as_eqs)):
                self.__intercepts.append(sp.solve([eval(self.__constraints_as_eqs[i]),eval(self.__constraints_as_eqs[j])]))
        print(self.__intercepts)
    
    def eval_possible_solutions(self):
        ''' Loops through "__intercepts" and evaluates the points. If a point satisfies all restrictions, it is saved to "__posible_solutions"
        '''
        for intercept in self.__intercepts:
            check = []
            # print(intercept[x],intercept[y])
            for i in self.__constraints:
                res = i.replace('x',str(intercept[x])).replace('y',str(intercept[y]))
                check.append(eval(res))
            print(check)
            if all(check):
                self.__posible_solutions.append(intercept)
    
    def find_solution(self):
        '''Loops throgh solution points and evaluates them in the objective function to find the answer to the problem.
        '''
        
        if self.__posible_solutions == []:
            self.__solution = "No possible solution for this problem"
        else:
            solutions = []
            for posible_solution in self.__posible_solutions:
                res = self.__objective.replace('x',str(posible_solution[x])).replace('y',str(posible_solution[y]))
                solutions.append(eval(res))
            if self.__maxmin == 0:
                self.__solution = max(solutions)
                self.__solution_points = self.__posible_solutions[solutions.index(max(solutions))]
                
            else:
                self.__solution = min(solutions)
                self.__solution_points = self.__posible_solutions[solutions.index(min(solutions))]
        print(f"Solution Points: {self.__solution_points}")
        print(f'Solution: {self.__solution}')
    
    def solve(self):
        self.constraints_to_eqs()
        self.find_intercepts()
        self.eval_possible_solutions()
        self.find_solution()

    def find_binding_constraints(self):
        for i in range(0,len(self.__constraints_as_eqs)-1):
            for j in range(i+1,len(self.__constraints_as_eqs)):
                if (self.__solution_points == sp.solve([eval(self.__constraints_as_eqs[i]),eval(self.__constraints_as_eqs[j])])):
                    self.__binding_constraints = [self.__constraints_as_eqs[i],self.__constraints_as_eqs[j]]
                    print(f"Binding Constraints: {self.__constraints_as_eqs[i]} , {self.__constraints_as_eqs[j]}")

    def find_binding_slopes(self):
        for constraint in self.__binding_constraints:
            solved = sp.solve(eval(constraint),y)
            print(solved)
            if solved != []:
                slope = sp.diff(solved[0])
                self.__binding_slopes.append(slope)
        self.__binding_slopes.sort()
        print(self.__binding_slopes)

    def find_sensitivity(self):
        x_coeficient = str(eval(self.__objective.replace('y','0').replace('x','1')))
        y_coeficient = str(eval(self.__objective.replace('x','0').replace('y','1')))
        print(f'x coeficient: {x_coeficient}')
        print(f'y coeficient: {y_coeficient}')
        slope = '-x/y'
        x_ssss = slope.replace('y',y_coeficient)
        x_increase = sp.solve([eval(f'{self.__binding_slopes[0]}<={x_ssss}')],x)
        x_decrease = sp.solve([eval(f'{x_ssss}<={self.__binding_slopes[1]}')],x)
        print(f'{self.__binding_slopes[0]}<={x_ssss}: {x_increase}')
        print(f'{x_ssss}<={self.__binding_slopes[1]}: {x_decrease}')
        y_ssss = slope.replace('x',x_coeficient)
        y_increase = sp.solve([eval(f'{self.__binding_slopes[0]}<={y_ssss}')],y)
        y_decrease = sp.solve([eval(f'{y_ssss}<={self.__binding_slopes[1]}')],y)
        print(f'{self.__binding_slopes[0]}<={y_ssss}: {y_increase}')
        print(f'{y_ssss}<={self.__binding_slopes[1]}: {y_decrease}')

