# Juan Jose Elgueta
import sys
import numpy as np
import sympy as sp
import math



class Simplex:
    def __init__(self, objective:str, constraints:list, noneg: set, maxmin:bool, binary:bool, entero: bool):
        ''' Constructor for the class.
            Args:
                objective (string): Objective function of the problem as a string in numpy format.
                constraints (list of strings) : Problem constraints as a list of strings in numpy format.
                noneg (set of strings): the variables with non negativity constraints
                maxmin (boolean): 1 if the problem is a maximization, 0 if it is a minimization.
                binary (boolean): 1 if the problem is suposed to be binary 0 if not.
                entero (boolean): 1 if the problem is suposed to be integwe 0 if not.
        '''
        self.__min_max: bool = maxmin
        self.__simplex_matrix: np.array=None
        self.__constraints: list = constraints
        self.__standard: list = []
        self.__objective: str = objective
        self.__constraints_arr: np.array= []
        self.__objective_arr: np.array= []
        self.__symbols = []
        self.__symbols_dict = {}
        self.__noneg= noneg
        self.__pivot_column=None
        self.__pivot_row=None
        self.tableaus = 0 
        self.stage = 1
        self.artificials = False
        self.binary = binary
        self.entero = entero
        self.cutting_planes = 0
        self.gomory_tableaus = 0

    @property
    def simplex_matrix(self):
        return self.__simplex_matrix

    def get_symbols_dict(self):
        '''Gets the symbols in the equations and restrictions and saves them in a dictionary to be used by sympy.
        '''
        for constraint in self.__standard:
            eq = sp.Eq(sp.sympify(constraint),0)
            for i in eq.free_symbols:
                self.__symbols_dict[i.name] = i

    def standard_form(self):
            ''' Rewrites restrictions in standard form and stores them in  "__constraints_standard" as an array of strings.
                Rewrites objective function in standard form and stores it in __objective
            '''
            count = 0
            art_count = 0
            lhs = []
            symbols_dict= {}

            for constraint in self.__constraints:
                count = count + 1
                if '<=' in constraint:
                    constraint = f's{count} +' + constraint + "*zconst"
                    lhs.append(constraint.replace('<=','-'))
                elif '>=' in constraint:
                    constraint = f'-s{count} + a{count} +' + constraint + "*zconst"
                    lhs.append(constraint.replace('>=','-'))
                elif '=' in constraint:
                    constraint = f'a{count} + ' + constraint + "*zconst"
                    lhs.append(constraint.replace('=','-'))

            if (self.__min_max == 1):
                self.__objective = f"-1*({self.__objective})"

            lhs.append(self.__objective)


            for constraint in lhs:
                eq = sp.Eq(sp.sympify(constraint),0)
                for i in eq.free_symbols:
                    symbols_dict[i.name] = i


            for variable in symbols_dict.keys():
                if 's' not in variable and 'a' not in variable and variable is not 'zconst':
                    if variable not in self.__noneg:
                        art_count = art_count+1
                        for constraint in lhs:
                            lhs[lhs.index(constraint)] = constraint.replace(variable, f'(y{art_count} - w{art_count})')

            for constraint in lhs:
                if sp.diff(constraint, symbols_dict['zconst']) > 0:
                    self.__standard.append(f"-({constraint})")
                else:
                    self.__standard.append(f"{constraint}")

            z_prime = ''
            for constraint in self.__standard[:-1]:
                if 'a' in constraint:
                    z_prime = z_prime+'+'+constraint
            
            if z_prime is not '':
                self.artificials = True
                self.__standard.append(z_prime)
            

            # print(self.__standard)


    def get_symbols(self):
        ''' Makes a sorted array of the names of the symbols to be used as a guide.
        '''
        syms = sorted(self.__symbols_dict.keys())
        self.__symbols = syms
        print(syms)

    def standard_to_array(self):
        ''' Rewrites the entire problem but only storing the coefficients of the variables.
            The constraints are stored in __constraints_arr and the objective function in __objective_arr.
        '''
        for constraint in self.__standard:
            eq = sp.Eq(sp.sympify(constraint),0)
            const_syms = eq.free_symbols
            const_syms = [str(i) for i in const_syms]
            const_syms = sorted(const_syms)
            constraint_arr = []
            for sym in self.__symbols:
                if sym not in const_syms:
                    constraint_arr.append(0)
                else:
                    coefficient = sp.diff(constraint,self.__symbols_dict[sym])
                    if sym is 'zconst':
                        constraint_arr.append(-(coefficient))
                    else:
                        constraint_arr.append((coefficient))
            self.__constraints_arr.append(constraint_arr)
        

    def fix_zprime(self):
        '''Adds needed zeroes for artificial vatiables to the zprime line 
        '''
        for i in range(len(self.__symbols)):
            if 'a' in self.__symbols[i]:
                self.__constraints_arr[-1][i] = 0
        # print(self.__constraints_arr)

    def create_simplex_matrix(self):
        '''Creates the simplex matrix based on the data found in __restrictions_arr and __objective.
        '''
        self.__simplex_matrix = np.array(self.__constraints_arr)
        self.__simplex_matrix[-1] = self.__simplex_matrix[-1]*-1

    def check_binary(self):
        ''' Checks if the coefficients were correctly entered as binary.
        '''
        if self.binary:
            for n in self.__simplex_matrix:
                print(n[:-1])
                for i in n[:-1]:
                    if(i != 1 and i != 0 and i != -1):
                        print(i)
                        return False
            return True
        else:
            return True



    def pivot_column(self,line):
        ''' Gets the pivot column
            Args:
                line (int): the index of the line to be used to find the pivot column
        '''
        minimo = min(self.__simplex_matrix[line][:-1])
        if minimo < 0:
            for num in range(len(self.__simplex_matrix[line])-1):
                if self.__simplex_matrix[line][num] == minimo:
                    column = num
                    break
            return column

    def pivot_row(self,line):
        ''' Gets the pivot row
            Args:
                line (int): the index of the line to be used to find the pivot row
        '''
        ratios = []
        for arr in self.__simplex_matrix[0:line]:
            if arr[self.pivot_column(line)] > 0 and arr[-1]>0:
                ratios.append((arr[-1]/arr[self.pivot_column(line)]))
            else:
                ratios.append(sys.maxsize)

        minimo = min(ratios)
        for num in range(len(ratios)):
            if ratios[num] == minimo:
                row = num
        return row


    
        

    def dual_row(self):
        minimo = min(self.__simplex_matrix[:-1,-1])
        if minimo < 0:
            for num in range(len(self.__simplex_matrix[:-1,-1])):
                if self.__simplex_matrix[num,-1] == minimo:
                    column = num
                    break
            return column
        

    def dual_column(self, column):
        ratios = []
        for i in range(len(self.simplex_matrix[column])-1):
            if self.simplex_matrix[column,i] < 0:
                ratios.append(abs(self.simplex_matrix[-1,i]/self.simplex_matrix[column,i]))
            else:
                ratios.append(sys.maxsize)

        minimo = min(ratios)
        for num in range(len(ratios)):
            if ratios[num] == minimo:
                row = num
        return row

    def choose_restriction(self):
        if self.artificials:
            line = -2
        else:
            line = -1 
        
        non_integers = []
        for i in self.simplex_matrix[:line,-1]:
            if isinstance(i,int):
                non_integers.append(0)
            else:
                non_integers.append(i)

        for i in range(len(non_integers)):
            non_integers[i] = non_integers[i] - (non_integers[i].numerator()//non_integers[i].denominator())

        maximo = max(non_integers)
        for i in range(len(self.simplex_matrix[:line,-1])):
            if (self.simplex_matrix[i,-1]- (self.simplex_matrix[i,-1].numerator()//self.simplex_matrix[i,-1].denominator())) == maximo:
                return self.simplex_matrix[i]


    def add_cutting_plane(self, restriction):

        if self.artificials:
            line = -2
        else:
            line = -1

        cutting_plane = []

        for i in restriction:
            cutting_plane.append(-(i - (i.numerator()//i.denominator())))
        
        count = 0
        for i in self.__symbols:
            if "s" in i:
                count = count + 1

        self.__symbols.insert(-1,f"s{count}")
        cutting_plane.insert(-1, 1)
        self.__simplex_matrix = np.insert(self.__simplex_matrix, -1, 0, axis = 1)
        self.__simplex_matrix = np.insert(self.__simplex_matrix, line, cutting_plane, axis = 0)

    def cutting_plane_needed(self):
        condition = []

        for i in self.simplex_matrix[:, -1]:
            condition.append(i%1 == 0)
        if  all(condition):
            return False
        else:
            return True
             

    def dual_pivot_needed(self):
        if all(i>=0 for i in self.simplex_matrix[:-1,-1]):
            return False
        else:
            return True


    def pivot_needed(self, line):
        ''' Checks if a new pivot is needed
            Args:
                line (int): the line to be checked for a new pivot.
        '''
        if all(i >=0 for i in self.simplex_matrix[line][:-1]):
            return False
        else:
            return True


    def pivot(self):
        ''' Performs the whole pivot operation and advances the simplex method.
        '''
        np.printoptions(precision=4, suppress=True, formatter={'float': '{:0.4f}'.format}, linewidth=100) 
        if self.tableaus == 0:
            print("--------------------------------------------- Start of Simplex -----------------------------------------")
        line = -self.stage
        print(f"Simplex Tableau {self.tableaus}")
        print(self.__symbols)
        print(self.simplex_matrix)
        self.tableaus = self.tableaus+1


        if self.pivot_needed(line):
            prow = self.pivot_row(line)
            pcolumn = self.pivot_column(line)
            print(f'Fila pivote: {prow}')
            print(f'Columna pivote : {pcolumn}')

            pivot_element = self.simplex_matrix[prow][pcolumn]
            if pivot_element != 1:
                self.__simplex_matrix[prow] = self.simplex_matrix[prow]/pivot_element


            for row in range(len(self.simplex_matrix)):
                if row != prow:
                    subs = self.simplex_matrix[prow]*self.__simplex_matrix[row][pcolumn]
                    self.__simplex_matrix[row] = self.__simplex_matrix[row] - subs
            self.pivot()

            
        else:
            if (self.stage==1 and self.artificials):
                self.stage = self.stage+1
                self.pivot()
            else:
                print("finished")
                return True





    def gomory(self):   
        self.gomory_tableaus = self.gomory_tableaus+1
        print(f"Gomory Tableau {self.gomory_tableaus}")
        print(self.__symbols)
        print(self.simplex_matrix)
 
        print(self.dual_pivot_needed())

        if self.dual_pivot_needed():
            prow = self.dual_row()
            pcolumn = self.dual_column(self.dual_row())
            print(self.simplex_matrix[prow,pcolumn])
            pivot_element = self.simplex_matrix[prow][pcolumn]
            if pivot_element != 1:
                self.__simplex_matrix[prow] = self.simplex_matrix[prow]/pivot_element
            
            for row in range(len(self.simplex_matrix)):
                if row != prow:
                    subs = self.simplex_matrix[prow]*self.__simplex_matrix[row][pcolumn]
                    self.__simplex_matrix[row] = self.__simplex_matrix[row] - subs

            self.gomory()

        else:
            print("finished gomory")


    def cut_that_plane(self):
        if self.cutting_plane_needed():
            self.cutting_planes = self.cutting_planes+1
            print(f"--------------------------------------- new cutting plane: {self.cutting_planes} --------------------------------------------")
            restriction = self.choose_restriction()
            self.add_cutting_plane(restriction)
            self.gomory()
            self.cut_that_plane()
        else:
            print("")


    def get_variables(self):
        if self.artificials:
            line = -2
        else:
            line = -1
        

        values = {}
        for i in range(len(self.__symbols)):
            if "x" in self.__symbols[i]:
                for j in range(len(self.simplex_matrix[:line,i])):
                    if self.simplex_matrix[j,i] == 1:
                        values[self.__symbols[i]]= self.simplex_matrix[j,-1]
        print("Variable Values:")
        print(values)
        print("Optimal Value:")
        print(self.simplex_matrix[line,-1])



    def simplex(self):
        ''' Performs the whole algorithm of the simplex method.
        '''
        self.standard_form()
        self.get_symbols_dict()
        self.get_symbols()
        self.standard_to_array()
        self.fix_zprime()
        self.create_simplex_matrix()
        if self.check_binary():
            self.pivot()
        else:
            print("not binary :(")
        if self.entero:
            self.cut_that_plane()
        self.get_variables()
            

#------------ ejmplo 1 -------------
# obj = '6*x1 + 4*x2'
# cons1 = '2*x1 + 3*x2 <= 30'
# cons2 = '1*x1 +  1*x2 >= 3'
# cons3 = '3*x1 + 2*x2 <= 24'

#------------ ejmplo 2 -------------
# obj = '4*x1 + 3*x2 +  2*x3'
# cons1 = '1*x1 + 2*x2 + 3*x3 <= 6'
# cons2 = '2*x1 + 1*x2 + 1*x3 <= 3'
# cons3 = '1*x1 + 1*x2 + 1*x3 <= 2'

#------------ ejmplo 3 -------------
# obj = '3*x1 + 9*x2'
# cons1 = '2*x1 + 1*x2 >= 8'
# cons2 = '1*x1 + 2*x2 >= 8'


#----------- ejemplo 4 -------------
# obj = '1*x1 + 1*x2 + 1*x3 + 1*x4'
# cons1 = '1*x1 + 1*x3 + 1*x4 <= 20'
# cons2 = '1*x2 + 1*x3 <= 40'
# cons3 = '1*x1 + 1*x2 >= 8'
# cons4 = '1*x2 + 1*x3 >= 15'

#----------- ejemplo 5 ------------
# obj = '1*x1 + 1*x2 + 1*x3'
# cons1 = '1*x1 + 1*x2 <= 20'
# cons2 = '1*x1 + 1*x3 <= 40'
# cons3 = '1*x1 + 1*x2 + 1*x3 >= 10'

# constr = [cons1,cons2,cons3]
# aa = Simplex(obj, constr,{'x1','x2','x3'} ,1,1,0)
# aa.simplex()
# print(aa.simplex_matrix[:-1,-1])



#------------- ejemplo enteras --------------
obj = "10*x1 + 15*x2"
cons1 = '282*x1 + 400*x2 <= 2000'
cons2 = '4*x1 + 40*x2 <= 140'
cons3 = '1*x1 <= 5'
constr = [cons1,cons2,cons3]
aa = Simplex(obj, constr,{'x1','x2','x3'} ,0,0,1)
aa.simplex()

