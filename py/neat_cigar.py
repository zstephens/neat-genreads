import re


class CigarString:
    """
    Creates and modifies CigarString for output
    """
    def __init__(self, string_in=None, list_in=None):

        if string_in is None and list_in is None:
            print('\nError: CigarString object not initialized.\n')
            exit(1)

        self.cigar_data = []

        if string_in is not None:
            self.join_cigar(j_string_in=string_in)

        if list_in is not None:
            self.join_cigar(j_list_in=list_in)

    @staticmethod
    def string_to_list(string_to_covert: str) -> list:
        """
        Converts cigar string to a full list

        :param string_to_covert: Cigar string to convert to list
        :return: cigar string in list format
        """

        cigar_dat = []
        letters = re.split(r"\d+", string_to_covert)[1:]
        numbers = [int(n) for n in re.findall(r"\d+", string_to_covert)]
        d_reserve = 0
        for i in range(len(letters)):
            if letters[i] == 'D':
                d_reserve = numbers[i]
            if letters[i] == 'M' or letters[i] == 'I':
                if d_reserve:
                    cigar_dat += ['D' * d_reserve + letters[i]] + [letters[i]] * (int(numbers[i]) - 1)
                else:
                    cigar_dat += [letters[i]] * int(numbers[i])
                d_reserve = 0
        return cigar_dat

    @staticmethod
    def list_to_string(list_to_convert: list) -> str:
        """
        Convert a cigar string in list format to a standard cigar string

        :param list_to_convert: Cigar string in list format
        :return: cigar string in string format
        """

        symbols = ''
        current_sym = list_to_convert[0]
        current_count = 1
        if 'D' in current_sym:
            current_sym = current_sym[-1]
        for k in range(1, len(list_to_convert)):
            next_sym = list_to_convert[k]
            if len(next_sym) == 1 and next_sym == current_sym:
                current_count += 1
            else:
                symbols += str(current_count) + current_sym
                if 'D' in next_sym:
                    symbols += str(next_sym.count('D')) + 'D'
                    current_sym = next_sym[-1]
                else:
                    current_sym = next_sym
                current_count = 1
        symbols += str(current_count) + current_sym
        return symbols

    def get_list(self):

        return self.cigar_data

    def get_string(self):

        return self.list_to_string(self.cigar_data)

    def join_cigar(self, j_string_in: str = None, j_list_in: list = None):
        """
        Adds an inputted cigar string in string or list format to the existing cigar string

        :param j_string_in: cigar string in string format
        :param j_list_in: cigar string in list format
        :return: None
        """

        if j_string_in is None and j_list_in is None:
            print('\nError: Invalid join operation in CigarString\n')
            exit(1)

        if j_string_in is not None:
            self.cigar_data += self.string_to_list(j_string_in)

        if j_list_in is not None:
            self.cigar_data += j_list_in

    def insert_cigar_element(self, pos: int, i_string_in: str = None, i_list_in: list = None):
        """
        Inserts a cigar string in either string or list format to the existing cigar string at position pos.

        :param pos: integer position where to insert the input cigar string
        :param i_string_in: input cigar string in string format
        :param i_list_in: input cigar string in list format
        :return:
        """

        if i_string_in is None and i_list_in is None:
            print('\nError: Invalid insertion operation in CigarString\n')
            exit(1)

        if pos < 0 or pos >= len(self.cigar_data):
            print('\nError: Invalid insertion position in CigarString\n')
            exit(1)

        if i_string_in is not None:
            self.cigar_data = self.cigar_data[:pos] + self.string_to_list(i_string_in) + self.cigar_data[pos:]

        if i_list_in is not None:
            self.cigar_data = self.cigar_data[:pos] + i_list_in + self.cigar_data[pos:]


if __name__ == '__main__':
    print('testing CigarString class...')

    str1 = '50M10D7I23M'
    str2 = '10I25M'
    iPos = 20
    my_cigar = CigarString(string_in=str1)
    my_cigar.insert_cigar_element(iPos, i_string_in=str2)
    assert(my_cigar.get_string() == "20M10I55M10D7I23M")
    print("passed")
