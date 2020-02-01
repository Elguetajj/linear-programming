import numpy as np
import sympy as sp

x,y = sp.symbols('x,y')

class LinearSolver:

    def __init__(self, constraints, objective,maxmin):
        '''inicializa la clase.

            Args:
                constraints (list of strings) : restricciones del problema en una lista de strings igual al formato numpy
                objective (string): la función objetivo en string con formato de numpy 

            Returns:
                Instancia de la clase Linear Solver 
        '''
        self.__constraints: list = constraints
        self.__constraints_as_eqs: list = []
        self.__objective: str = objective
        self.__intercepts: list = []
        self.__posible_solutions = []
        self.__solution_points: dict = None
        self.__solution: int = None
        self.__maxmin: str = maxmin

    def constraints_to_eqs(self):
        ''' Convierte las restricciones en igualdades a 0 y las guarda en el atributo privado "__constraints_as_eqs" como un array de strings.
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
        ''' Utiliza los elementos del atributo privado "__constraints_as_eqs" para encontrar los interceptos entre todas las restricciones y los guarda en el atributo privado "__intercepts" del objeto
        '''
        for i in range(0,len(self.__constraints_as_eqs)-1):
            for j in range(i+1,len(self.__constraints_as_eqs)):
                self.__intercepts.append(sp.solve([eval(self.__constraints_as_eqs[i]),eval(self.__constraints_as_eqs[j])]))
        print(self.__intercepts)
    
    def eval_possible_solutions(self):
        ''' Recorre los interceptos en "__intercepts" y los evalua contra las restricciones. Si los interceptos cumplen se guardan en el atributo privado "__posible_solutions"
        '''
        for intercept in self.__intercepts:
            check = []
            print(intercept[x],intercept[y])
            for i in self.__constraints:
                res = i.replace('x',str(intercept[x])).replace('y',str(intercept[y]))
                check.append(eval(res))
            print(check)
            if all(check):
                self.__posible_solutions.append(intercept)
    
    def solution_points(self):
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
        print('solution points:')
        print(self.__solution_points)
        print('solution:')
        print(self.__solution)
    
    def solve(self):
        self.constraints_to_eqs()
        self.find_intercepts()
        self.eval_possible_solutions()
        self.solution_points()



# print(reduce_rational_inequalities([[eval(eq),eval(eq2)]],x))





