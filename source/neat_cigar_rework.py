import sys
import pdb
from itertools import groupby
from operator import itemgetter
from source.neat_cigar import CigarString as CigarStringOld


class CigarString:
    """
    This code borrowed in whole from https://github.com/brentp/cigar

    cigar is a simple library for dealing with cigar strings. the most useful
    feature now is soft-masking from left or right. This allows one to adjust
    a SAM record only by changing the cigar string to soft-mask a number of bases
    such that the rest of the SAM record (pos, tlen, etc.) remain valid, but
    downstream tools will not consider the soft-masked bases in further analysis.

    CigarString adds additional functionality to the cigar module, to make it compatible with
    the rest of the NEAT codebase. Many functions are lifted directly from brentp's repo
    """

    read_consuming_ops = ("M", "I", "S", "=", "X")
    ref_consuming_ops = ("M", "D", "N", "=", "X")

    def __init__(self, cigar_string):
        self.cigar = cigar_string

    def items(self):
        if self.cigar == "*":
            yield 0, None
            raise StopIteration
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int("".join(n)), "".join(next(cig_iter)[1])

    def __str__(self):
        return self.cigar

    def __repr__(self):
        return "CigarString('%s')" % self

    def __len__(self):
        """
        sum of MIS=X ops shall equal the sequence length.
        """
        return sum(l for l, op in self.items() if op in CigarString.read_consuming_ops)

    @classmethod
    def string_from_elements(cls, elements):
        return "".join("%i%s" % (l, op) for l, op in elements if l != 0)

    # I moved this from the main Cigar class in order to make the code more consistent
    def merge_like_ops(self):
        """
        >>> CigarString("1S20M").merge_like_ops()
        CigarString('1S20M')
        >>> CigarString("1S1S20M").merge_like_ops()
        CigarString('2S20M')
        >>> CigarString("1S1S1S20M").merge_like_ops()
        CigarString('3S20M')
        >>> CigarString("1S1S1S20M1S1S").merge_like_ops()
        CigarString('3S20M2S')
        """

        cigs = []
        for op, grps in groupby(self.items(), itemgetter(1)):
            cigs.append((sum(g[0] for g in grps), op))

        return CigarString(self.string_from_elements(cigs))

    def find_position(self, position: int) -> (int,int):
        """
        This finds the index of the element of list(self.items()) that contains the starting position given.

        :param position: The position to find of the CigarString
        :return: the index of that in terms of the elements of list(self.items())
        >>> cigar = CigarString('2064M1I12979M3D18439M')
        >>> position = 15414
        >>> index, bases = cigar.find_position(position)
        >>> assert(index == 4)
        >>> assert(bases == 370)
        """
        consumed_bases = 0
        bases_left = position
        index = 0
        for item in self.items():
            if item[1] == 'D':
                index += 1
                continue
            consumed_bases += item[0]
            if consumed_bases >= position:
                return index, bases_left
            else:
                bases_left -= item[0]
                index += 1
        return index, bases_left

    def insert_cigar_element(self, pos1: int, insertion_cigar: 'CigarString', pos2: int = None) -> None:
        """
        Inserts a cigar string in either string or list format to the existing cigar string at position pos.

        :param insertion_cigar: A cigar to insert into current cigar string
        :param pos1: integer position where to insert the input cigar string
        :param pos2: integer position where to pick up the remainder
        :return: None

        >>> str1 = CigarString('50M10D7I23M')
        >>> str2 = CigarString('10I25M')
        >>> iPos = 55
        >>> str1.insert_cigar_element(iPos, str2)
        >>> assert(str1.cigar == "50M10D15I25M2I23M")

        >>> str1 = CigarString('50M10D7I23M')
        >>> iPos = 20
        >>> str1.insert_cigar_element(iPos, str2)
        >>> assert(str1.cigar == "20M10I55M10D7I23M")

        >>> str1 = CigarString('11100M')
        >>> str2 = CigarString('2I')
        >>> iPos = 5785 + 1
        >>> str1.insert_cigar_element(iPos, str2)
        >>> assert (str1.cigar == "5786M2I5314M")

        >>> str1 = CigarString('11100M')
        >>> str2 = CigarString('1D')
        >>> iPos = 6610 + 1
        >>> str1.insert_cigar_element(iPos, str2)
        >>> assert(len(str1) == 11100)
        >>> assert (str1.cigar == "6611M1D4489M")

        >>> str1 = CigarString('10M1I2D10M')
        >>> str2 = CigarString('1D')
        >>> pos1 = 10
        >>> pos2 = pos1 + 2
        >>> str1.insert_cigar_element(pos1, str2, pos2)
        >>> cigar_to_insert = '1D1M'
        >>> str1_list = str1.string_to_list()
        >>> str1_list[pos1] = cigar_to_insert
        >>> test_string = CigarStringOld.list_to_string(str1_list)
        >>> assert(str1.cigar == test_string)

        >>> original_string = '100M1I10M'
        >>> original_cigar = CigarString(original_string)
        >>> v_pos = 9
        >>> v_pos2 = 9 + 6
        >>> cigar_to_insert = CigarString('6I')
        >>> original_cigar.insert_cigar_element(v_pos, cigar_to_insert, v_pos2)
        >>> cigar_list = CigarStringOld.string_to_list('100M1I10M')
        >>> indel_length = 6
        >>> test_string = cigar_list[:v_pos] + ['I'] * indel_length + cigar_list[v_pos2:]
        >>> test_string = CigarStringOld.list_to_string(test_string)
        >>> assert (original_cigar.cigar == test_string)
        """

        if insertion_cigar is None:
            print('\nError: Invalid insertion operation in CigarString\n')
            sys.exit(1)

        if pos1 < 0 or pos1 >= len(self):
            print('\nError: Invalid insertion position in CigarString\n')
            sys.exit(1)

        if pos2 and (pos2 < pos1 or pos2 >= len(self)):
            print('\nError: Invalid second position\n')
            sys.exit(1)

        try:
            found = False
            start_index, bases_remain = self.find_position(pos1)
            extra_bases_to_discard = 0
            if pos2:
                extra_bases_to_discard = pos2 - pos1
            # TODO search for all instances of insert and get fragment and make sure they all work
            new_element = list(self.items())[:start_index]
            items_to_iterate = list(self.items())[start_index:]
            new_item = ((bases_remain, items_to_iterate[0][1]), (items_to_iterate[0][0]
                                                                 - bases_remain, items_to_iterate[0][1]))
            new_element.append(new_item[0])
            items_to_iterate[0:1] = new_item[1:]
            new_element += [item for item in insertion_cigar.items()]
            for item in items_to_iterate:
                if item[0] == 0:
                    continue
                if extra_bases_to_discard <= 0 and item[1] == 'D':
                    new_element.append(item)
                elif extra_bases_to_discard > 0 and item[1] == 'D':
                    continue
                elif found:
                    new_element.append(item)
                elif extra_bases_to_discard > item[0]:
                    extra_bases_to_discard -= item[0]
                elif extra_bases_to_discard == item[0]:
                    extra_bases_to_discard = 0
                    found = True
                elif extra_bases_to_discard < item[0]:
                    new_element.append((item[0] - extra_bases_to_discard, item[1]))
                    extra_bases_to_discard = 0
                    found = True

            # The point of this block of code is just to trim off any starting or ending "D"s, so we don't have
            # something like "100M1D" as the return.
            ret = new_element
            for i in range(len(new_element)):
                if new_element[0][1] == 'D':
                    ret = new_element[1:]
                if new_element[-1][1] == 'D':
                    ret = new_element[:-1]
            new_string = self.string_from_elements(ret)
            self.cigar = CigarString(new_string).merge_like_ops().cigar
        except ValueError:
            print('\nBug: Problem with insertion.\n')
            sys.exit(1)

    def get_cigar_fragment(self, start: int, end: int) -> 'CigarString':
        """
        Take a slice of a cigar string. E.g., if we have a cigar string that is "1000M" and we want a slice of the
        first 101 characters, then the return should be '101M'. If we had "20M100I20M" and want the
        first 101 characters, it would return '20M81I'.

        :param start: start point of sequence to retrieve
        :param end: end point of sequence to retrieve
        :return: The sequence that spans start to end
        >>> str1 = CigarString('10M5D10M')
        >>> start = 0
        >>> end = 14
        >>> frag = str1.get_cigar_fragment(start, end)
        >>> assert(len(frag) == 14)
        >>> assert(frag.cigar == "10M5D4M")

        >>> str1 = CigarString('90M1D10M')
        >>> start = 0
        >>> end = 100
        >>> frag = str1.get_cigar_fragment(start, end)
        >>> assert(frag.cigar == "90M1D10M")

        >>> str1 = CigarString('10M2D10M')
        >>> start = 10
        >>> end = 12
        >>> frag = str1.get_cigar_fragment(start, end)
        >>> assert(frag.cigar == "2M")

        >>> str1 = CigarString('10M1D10M')
        >>> start = 10
        >>> end = 12
        >>> frag = str1.get_cigar_fragment(start, end)
        >>> assert(frag.cigar == "2M")

        >>> str1 = CigarString('102M2I10000M')
        >>> start1 = 1
        >>> end1 = 102
        >>> frag1 = str1.get_cigar_fragment(start1, end1)
        >>> start2 = 102
        >>> end2 = 203
        >>> frag2 = str1.get_cigar_fragment(start2, end2)
        >>> assert(frag1.cigar == "101M")
        >>> assert(frag2.cigar == "2I99M")

        >>> temp_symbol_string = CigarString('25179M1D8304M')
        >>> start = 25080
        >>> end = 25180
        >>> frag = temp_symbol_string.get_cigar_fragment(start, end)
        >>> assert(frag.cigar == '99M1D1M')

        >>> temp_symbol_string = CigarString('25179M1D8304M')
        >>> start = 25079
        >>> end = start + 100
        >>> frag = temp_symbol_string.get_cigar_fragment(start, end)
        >>> assert(frag.cigar == "100M")
        """
        # Minus 1 because python slices don't include the end coordinate
        window_size = end - start
        remaining_window = window_size
        if window_size < 0:
            print(f'\nError: start and end coordinates for get_cigar_fragment '
                  f'are wrong: start: {start}, end: {end}')
        if start is None or end is None:
            print('\nError: Invalid coordinates in get_cigar_fragment\n')
            sys.exit(1)

        if start < 0 or start >= len(self):
            print('\nError: Invalid start point in get_cigar_fragment\n')
            sys.exit(1)

        if end > len(self) or end < 0:
            print('\nError: Invalid end point in get_cigar_fragment')
            sys.exit(1)

        try:
            start_index, bases_remain = self.find_position(start)
            end_index = self.find_position(end)[0]
            new_element = []
            # +1 here to make it inclusive
            list_of_items = list(self.items())[start_index:end_index+1]
            if len(list_of_items) == 1:
                new_element.append((window_size, list_of_items[0][1]))
            else:
                for item in list_of_items:
                    if sum([item[0] for item in new_element if item[1] != 'D']) == window_size:
                        break
                    if item[1] == 'D':
                        new_element.append(item)
                    elif bases_remain > 0:
                        bases_remain = item[0] - bases_remain
                        if bases_remain < 0:
                            print("Something went wrong retrieving fragment")
                            # pdb.set_trace()
                            sys.exit(1)
                        elif bases_remain == 0:
                            continue
                        new_element.append((bases_remain, item[1]))
                        remaining_window -= bases_remain
                        bases_remain = 0
                    elif remaining_window > item[0]:
                        new_element.append(item)
                        remaining_window -= item[0]
                    elif remaining_window <= item[0]:
                        new_element.append((remaining_window, item[1]))
                        break

            # The point of this block of code is just to trim off any starting or ending "D"s, so we don't have
            # something like "100M1D" as the return.
            ret = new_element
            for i in range(len(new_element)):
                if new_element[0][1] == 'D':
                    ret = new_element[1:]
                if new_element[-1][1] == 'D':
                    ret = new_element[:-1]
            new_string = self.string_from_elements(ret)
            return CigarString(new_string).merge_like_ops()

        except ValueError:
            print('\nBug: Problem retrieving fragment.\n')
            # pdb.set_trace()
            sys.exit(1)

    # TODO use this method or delete it.
    def count_elements(self, element: str) -> int:
        count = 0
        for item in self.items():
            if item[1] == element:
                count += item[0]
        return count

    def string_to_list(self) -> list:
        """
        Converts the cigar string to a full list.

        :return: cigar string in list format.
        """
        cigar_dat = []
        d_reserve = 0
        for item in self.items():
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


if __name__ == "__main__":
    str1 = CigarString('10M1I2D10M')
    str2 = CigarString('1D')
    pos1 = 10
    pos2 = pos1 + 2
    str1.insert_cigar_element(pos1, str2, pos2)
    cigar_to_insert = '1D1M'
    str1_list = str1.string_to_list()
    str1_list[pos1] = cigar_to_insert
    test_string = CigarStringOld.list_to_string(str1_list)
    assert (str1.cigar == test_string)



