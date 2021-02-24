from itertools import groupby


class CigarString:
    """"
    Now we're testing out a list method of CigarString that we're hoping is faster
    """
    _read_consuming_ops = ("M", "I", "S", "=", "X")
    _ref_consuming_ops = ("M", "D", "N", "=", "X")

    @staticmethod
    def items(string_in: str) -> iter:
        """
        iterator for cigar string items
        :return: Creates an iterator object
        """
        if string_in == "*":
            yield 0, None
            raise StopIteration
        cig_iter = groupby(string_in, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int("".join(n)), "".join(next(cig_iter)[1])

    @staticmethod
    def string_to_list(string_in: str) -> list:
        """
        This will convert a cigar string into a list of elements
        :param string_in: a valid cigar string.
        :return: a list version of that string.
        """
        cigar_dat = []
        d_reserve = 0
        for item in CigarString.items(string_in):
            if item[1] == 'D':
                d_reserve = item[0]
            if item[1] in ['M', 'I']:
                if d_reserve:
                    cigar_dat += ['D' * d_reserve + item[1]] + [item[1]] * (item[0] - 1)
                else:
                    cigar_dat += [item[1]] * item[0]
                d_reserve = 0
        return cigar_dat

    @staticmethod
    def list_to_string(input_list: list) -> str:
        """
        Convert a cigar string in list format to a standard cigar string
        :param input_list: Cigar string in list format
        :return: cigar string in string format
        """

        symbols = ''
        current_sym = input_list[0]
        current_count = 1
        if 'D' in current_sym:
            current_sym = current_sym[-1]
        for k in range(1, len(input_list)):
            next_sym = input_list[k]
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


if __name__ == '__main__':
    cigar = "10M1I3D1M"
    lst = CigarString.string_to_list(cigar)
    print(lst)
    st = CigarString.list_to_string(lst)
    print(st)
    