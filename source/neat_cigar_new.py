import sys
from itertools import groupby
from operator import itemgetter
from source.neat_cigar_old import CigarString as CigarStringOld
from source.neat_cigar import CigarString as CigarStringIntermediate
import timeit


class CigarString:
    """"
    Now we're testing out a list method of CigarString that we're hoping is faster
    """
    _read_consuming_ops = ("M", "I", "S", "=", "X")
    _ref_consuming_ops = ("M", "D", "N", "=", "X")

    def __init__(self, cigar_string: str = None, list_in: list = None):
        """
        Converting the cigar string to a list type method thingy as a test
        :param cigar_string: An actual string version of the cigar
        """
        self._needs_update = False
        self._cigar_string = "*"
        if cigar_string:
            self._cigar_string = cigar_string
        self._cigar = []
        if list_in:
            self._cigar = list_in
            self._needs_update = True
        else:
            d_reserve = 0
            for item in self.items():
                if item[1] == 'D':
                    d_reserve = item[0]
                if item[1] in ['M', 'I']:
                    if d_reserve:
                        self._cigar += ['D' * d_reserve + item[1]] + [item[1]] * (item[0] - 1)
                    else:
                        self._cigar += [item[1]] * item[0]
                    d_reserve = 0
            self.update_cig_string()

    def __str__(self):
        if self._needs_update:
            return "Needs Update, old string: " + self._cigar_string
        else:
            return self._cigar_string

    def __repr__(self):
        if self._needs_update:
            return "Needs update, old: CigarString('%s')" % self._cigar_string
        else:
            return "CigarString('%s')" % self._cigar_string

    def __len__(self):
        return self._cigar.__len__()

    def __getitem__(self, item):
        if type(item) == slice:
            ret = self._cigar.__getitem__(item)
            return CigarString(list_in=ret)
        return self._cigar.__getitem__(item)

    def __setitem__(self, key, value):
        self._needs_update = True
        self._cigar[key] = value

    def __getslice__(self, i, j):
        ret = self._cigar.__getslice__(i, j)
        print(type(ret))
        ret2 = CigarString(list_in=ret)
        return ret2

    def append(self, item):
        self._needs_update = True
        self._cigar.append(item)

    def __add__(self, other):
        self._needs_update = True
        result = self._cigar.__add__(other)
        try:
            return CigarString(list_in=result)
        except TypeError:
            return result

    def __mul__(self, other):
        return CigarString(list_in=self._cigar.__mul__(self, other))

    def items(self):
        """
        iterator for cigar string items
        :return: Creates an iterator object
        """
        if self._cigar_string == "*" or self._needs_update:
            yield 0, None
            raise StopIteration
        cig_iter = groupby(self._cigar_string, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int("".join(n)), "".join(next(cig_iter)[1])

    @property
    def cigar(self):
        """
        getter for cigar item
        :return: None
        """
        return self._cigar

    # setter for cigar item
    @cigar.setter
    def cigar(self, new):
        """
        Setter that updates cigar list to new list
        :param new: new cigar list
        :return: None
        """
        # TODO this validation isn't correct and may not be needed
        # if type(input) is not list or [l for l in input if l not in self._read_consuming_ops]:
        #     raise ValueError("Invalid input string to cigar")
        #     sys.exit(1)
        self._cigar = new
        self._needs_update = True

    def update_cig_string(self):
        """
        If the list has changed, this will update the corresponding string.
        Due to speed problems with converting list to string, we won't do this unless
        instructed, so note that that cigar_string is only valid if needs_update is also false
        :return: None

        >>> test_cig = CigarString("100M")
        >>> test_cig[10] = "I"
        >>> test_cig.update_cig_string()
        >>> assert(test_cig._cigar_string == "10M1I89M")
        """
        if self._needs_update:
            symbols = ''
            current_sym = self._cigar[0]
            current_count = 1
            if 'D' in current_sym:
                current_sym = current_sym[-1]
            for k in range(1, len(self._cigar)):
                next_sym = self._cigar[k]
                if len(next_sym) == 1 and next_sym == current_sym:
                    current_count += 1
                else:
                    symbols += str(current_count) + current_sym
                    if 'D' in next_sym:
                        symbols += str(next_sym.count('D')) + "D"
                        current_sym = next_sym[-1]
                    else:
                        current_sym = next_sym
                    current_count = 1
            symbols += str(current_count) + current_sym
            self._cigar_string = symbols
            self._needs_update = False
        else:
            pass


if __name__ == '__main__':
    cigar = CigarString('20M')
    new_thing = cigar[1:11]
    assert(type(cigar[1:11]) == CigarString)
    print(cigar.cigar)
    cigar.append("I")
    print(cigar.cigar)
    cigar[3] = "I"
    cigar.update_cig_string()
    print(cigar.cigar)
    print(cigar)
