import numpy as np
"""read a tab seperated data matrix and a file with row names and provide a mapping between  the two of them"""
# TODO make dtype a parameter to allow to read continuous valued arrays as well


class named_matrix():

    """class to represent a numpy array with row name and row index access"""

    def __init__(self, matrix_f, row_names_f):
        self.array = np.genfromtxt(matrix_f, filling_values="-1")
        if len(self.array.shape) == 1:
            #TODO not sure if change from dtype int to float might cause problems
            self.array = np.transpose(np.array(self.array, ndmin=2, dtype=float))
        self.rn_dict = self.parse_row_names(row_names_f)

    def parse_row_names(self, row_names_f):
        """read and parse a file with row names and map that to a integer index for self.array"""
        f = open(row_names_f, 'r')
        ll = f.readlines()
        f.close()
        rn_dict = dict([(ll[i].strip(), i) for i in range(len(ll))])
        return rn_dict

    def __getitem__(self, a):
        """provides row name and row index access, which is dynamically determined"""
        if isinstance(a, type(0)):
            return self.array[a]
        if isinstance(a, type('')):
            return self.array[self.rn_dict[a]]
        else:
            raise IndexError


if __name__ == "__main__":
    nm = named_matrix("test_matrix.txt", "test_names.txt")
    print nm['a']
    print nm[1]
