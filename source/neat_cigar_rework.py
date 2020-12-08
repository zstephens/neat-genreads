import sys
from itertools import groupby
from operator import itemgetter

class Cigar(object):
    """
    This code borrowed in whole from https://github.com/brentp/cigar

    cigar is a simple library for dealing with cigar strings. the most useful
    feature now is soft-masking from left or right. This allows one to adjust
    a SAM record only by changing the cigar string to soft-mask a number of bases
    such that the rest of the SAM record (pos, tlen, etc.) remain valid, but
    downstream tools will not consider the soft-masked bases in further analysis.
    >>> c = Cigar('100M')
    >>> len(c)
    100
    >>> str(c)
    '100M'
    >>> list(c.items())
    [(100, 'M')]
    >>> c = Cigar('20H20M20S')
    >>> len(c)
    40
    >>> str(c)
    '20H20M20S'
    >>> list(c.items())
    [(20, 'H'), (20, 'M'), (20, 'S')]
    >>> c.mask_left(29).cigar, c.cigar
    ('20H9S11M20S', '20H20M20S')
    >>> c = Cigar('10M20S10M')
    >>> c.mask_left(10).cigar
    '30S10M'
    >>> c.mask_left(9).cigar
    '9S1M20S10M'
    >>> Cigar('10S').mask_left(10).cigar
    '10S'
    >>> Cigar('10H').mask_left(10).cigar
    '10H'
    >>> Cigar('10H').mask_left(11).cigar
    '10H'
    >>> Cigar('10H').mask_left(9).cigar
    '10H'
    >>> Cigar('1M10H').mask_left(9).cigar
    '1S10H'
    >>> Cigar('5M10H').mask_left(9).cigar
    '5S10H'
    >>> c = Cigar('1S1H1S5H1S5M10H')
    >>> c.mask_left(9).cigar == c.cigar
    True
    >>> c = Cigar('1S1H1S5H1S5M10H')
    >>> c.mask_right(9).cigar == c.cigar
    True
    >>> c.mask_right(11).cigar
    '1S1H1S5H1S4M1S10H'
    """
    read_consuming_ops = ("M", "I", "S", "=", "X")
    ref_consuming_ops = ("M", "D", "N", "=", "X")

    def __init__(self, cigar_string):
        self.cigar = cigar_string

    def items(self):
        if self.cigar == "*":
            yield (0, None)
            raise StopIteration
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int("".join(n)), "".join(next(cig_iter)[1])

    def __str__(self):
        return self.cigar

    def __repr__(self):
        return "Cigar('%s')" % self

    def __len__(self):
        """
        sum of MIS=X ops shall equal the sequence length.
        """
        return sum(l for l, op,in self.items() \
                               if op in Cigar.read_consuming_ops)

    def reference_length(self):
        return sum(l for l, op in self.items() \
                               if op in Cigar.ref_consuming_ops)

    def mask_left(self, n_seq_bases, mask="S"):
        """
        Return a new cigar with cigar string where the first `n_seq_bases` are
        soft-masked unless they are already hard-masked.
        """
        cigs = list(self.items())
        new_cigs = []

        c, cum_len  = self.cigar, 0
        for i, (l, op) in enumerate(cigs):
            if op in Cigar.read_consuming_ops:
                cum_len += l
            if op == "H":
                cum_len += l
                new_cigs.append(cigs[i])
            elif cum_len < n_seq_bases:
                new_cigs.append(cigs[i])
            else:
                # the current cigar element is split by the masking.
                right_extra = cum_len - n_seq_bases
                new_cigs.append((l - right_extra, 'S'))
                if right_extra != 0:
                    new_cigs.append((right_extra, cigs[i][1]))
            if cum_len >= n_seq_bases: break
        else:
            pass

        new_cigs[:i] = [(l, op if op in "HS" else "S") for l, op in
                new_cigs[:i]]
        new_cigs.extend(cigs[i + 1:])
        return Cigar(Cigar.string_from_elements(new_cigs)).merge_like_ops()


    @classmethod
    def string_from_elements(self, elements):
        return "".join("%i%s" % (l, op) for l, op in elements if l !=0)


    def mask_right(self, n_seq_bases, mask="S"):
        """
        Return a new cigar with cigar string where the last `n_seq_bases` are
        soft-masked unless they are already hard-masked.
        """
        return Cigar(Cigar(self._reverse_cigar()).mask_left(n_seq_bases, mask)._reverse_cigar())


    def _reverse_cigar(self):
        return Cigar.string_from_elements(list(self.items())[::-1])

    def merge_like_ops(self):
        """
        >>> Cigar("1S20M").merge_like_ops()
        Cigar('1S20M')
        >>> Cigar("1S1S20M").merge_like_ops()
        Cigar('2S20M')
        >>> Cigar("1S1S1S20M").merge_like_ops()
        Cigar('3S20M')
        >>> Cigar("1S1S1S20M1S1S").merge_like_ops()
        Cigar('3S20M2S')
        """

        cigs = []
        for op, grps in groupby(self.items(), itemgetter(1)):
            cigs.append((sum(g[0] for g in grps), op))

        return Cigar(self.string_from_elements(cigs))


class CigarString(Cigar):
    """
    This adds odditional functionality to the cigar module above, to make it compatible with
    the rest of the NEAT codebase
    """

    def insert_cigar_element(self, pos: int,
                             insertion_cigar: 'CigarString',
                             length: int = 0) -> None:
        """
        Inserts a cigar string in either string or list format to the existing cigar string at position pos.

        :param length: If consuming spaces, length of cigar to be added
        :param insertion_cigar: A cigar to insert into current cigar string
        :param pos: integer position where to insert the input cigar string
        :return:
        """

        if insertion_cigar is None:
            print('\nError: Invalid insertion operation in CigarString\n')
            sys.exit(1)

        if pos < 0 or pos >= len(self):
            print('\nError: Invalid insertion position in CigarString\n')
            sys.exit(1)

        try:
            current_pos = 0
            previous_pos = 0
            new_element = []
            found = False
            for item in self.items():
                current_pos += item[0]
                if current_pos > pos and not found:
                    current_block = item[0]
                    current_index = pos - previous_pos
                    new_element.append((current_index, item[1]))
                    for insert in insertion_cigar.items():
                        new_element.append(insert)
                    if length != 0:
                        new_element.append((current_block-current_index-abs(length), item[1]))
                    else:
                        new_element.append((current_block-current_index, item[1]))
                    found = True
                else:
                    new_element.append(item)
                previous_pos = current_pos
            new_string = self.string_from_elements(new_element)
            self.cigar = CigarString(new_string).merge_like_ops().cigar
        except ValueError:
            print('\nBug: Problem with insertion.\n')
            sys.exit(1)

    def get_cigar_fragment(self, start: int, end: int) -> str:
        """
        Take a slice of a cigar string. E.g., if we have a cigar string that is "1000M" and we want a slice of the
        first 101 characters, then the return should be '101M'. If we had "20M100I20M" and want the
        first 101 characters, it would return '20M81I'.

        :param start: start point of sequence to retrieve
        :param end: end point of sequence to retrieve
        :return: The sequence that spans start to end
        """
        # Minus 1 because python slices don't include the end coordinate
        window_size = end - start
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
            # Minus 1 because of indexing
            current_pos = -1
            previous_pos = -1
            current_index = 0
            start_found = False
            ret = []
            for item in self.items():
                current_pos += item[0]
                current_block_size = item[0]
                current_symbol = item[1]
                if current_pos >= start and not start_found:
                    start_found = True
                    current_index = previous_pos + 1
                    # start - current index = The index where we start relative to the start of the block
                    # if there are not enough items in the current block, we will consume what is left
                    # and move to the next block
                    if current_block_size - (start - current_index) >= window_size:
                        # If from where we start to the end of the current block is sufficient for our window,
                        # then the return is just everything of that symbol
                        ret.append((window_size, current_symbol))
                        break
                    else:
                        # if there are not enough items in the current block, we will consume what is left
                        # and move to the next block
                        ret.append((current_block_size - (start - current_index), current_symbol))
                        # shrink the window to account for the items we've found
                        window_size -= current_block_size - (start - current_index)
                elif start_found:
                    if current_block_size >= window_size:
                        # If this block can cover the remaining slice, then append the end and break
                        ret.append((window_size, current_symbol))
                        break
                    else:
                        # If we still need more, then consume this entire block
                        ret.append((current_block_size, current_symbol))
                        # shrink the window to account for the items we've found
                        window_size -= current_block_size
                else:
                    pass

                previous_pos = current_pos

        except ValueError:
            print('\nBug: Problem retrieving fragment.\n')
            sys.exit(1)
        return Cigar.string_from_elements(ret)

    def count_elements(self, element: str) -> int:
        count = 0
        for item in self.items():
            if item[1] == element:
                count += item[0]
        return count


if __name__ == '__main__':
    print('testing CigarString class...')

    str1 = CigarString('50M10D7I23M')
    str2 = CigarString('10I25M')
    iPos = 55

    str1.insert_cigar_element(iPos, str2)
    print(str1.cigar)
    assert(str1.cigar == "50M5D10I25M5D7I23M")
    print("passed")

    str1 = CigarString('50M10D7I23M')
    iPos = 20
    str1.insert_cigar_element(iPos, str2)
    print(str1.cigar)
    assert(str1.cigar == "20M10I55M10D7I23M")
    print("passed")

    str1 = CigarString('11100M')
    str2 = CigarString('2I')
    iPos = 5785 + 1
    str1.insert_cigar_element(iPos, str2, 1)
    print(str1.cigar)
    assert (str1.cigar == "5786M2I5313M")
    print("passed")

    str1 = CigarString('11100M')
    str2 = CigarString('1D1M')
    iPos = 6610 + 1
    str1.insert_cigar_element(iPos, str2, -1)
    print(str1.cigar)
    assert (str1.cigar == "6611M1D4489M")
    print("passed")

    str1 = CigarString('10M5D10M')
    # [('10', 'M'), ('5', 'D'), ('10', 'M')]
    start = 0
    end = 14
    frag = str1.get_cigar_fragment(start, end)
    print(frag)
    assert(frag == "10M4D")

    str1 = CigarString('10M2D10M')
    # [(10, 'M'), (2, 'D'), (10, 'M')]
    start = 10
    end = 12
    frag = str1.get_cigar_fragment(start, end)
    print(frag)
    assert(frag == "2D")
    print("passed")

    str1 = CigarString('10M1D10M')
    # [(10, 'M'), (1, 'D'), (10, 'M')]
    start = 10
    end = 12
    frag = str1.get_cigar_fragment(start, end)
    print(frag)
    assert(frag == "1D1M")
    print("passed")

    str1 = CigarString('102M2I10000M')
    start1 = 1
    end1 = 102
    frag1 = str1.get_cigar_fragment(start1, end1)
    start2 = 102
    end2 = 203
    frag2 = str1.get_cigar_fragment(start2, end2)
    print(frag1)
    print(frag2)
    assert(frag1 == "101M")
    assert(frag2 == "2I99M")



