from source.neat_cigar_rework import CigarString as CigarStringNew
from source.neat_cigar import CigarString
from time import time

ploidy = 2
start_time = time()
temp_symbol_string = CigarStringNew('10777M1I322M')
cigar_string_new_time = time() - start_time
print(f'CigarStringNew time = {cigar_string_new_time}')

start_time = time()
temp_symbol_string2 = CigarString.string_to_list('10777M1I322M')
old_method_time = time() - start_time
print(f'String-to-list time = {old_method_time}')
fm_pos1 = [[] for _ in range(ploidy)]
fm_span1 = [[] for _ in range(ploidy)]
fm_pos2 = [[] for _ in range(ploidy)]
fm_span2 = [[] for _ in range(ploidy)]
span_dif = 0
read_len = 101

for i in range(ploidy):
    md_so_far = 0
    start_time = time()
    for item in temp_symbol_string.items():
        if item[1] in ['M', 'D']:
            for j in range(item[0]):
                fm_pos1[i].append(md_so_far)
                span_dif = CigarStringNew(temp_symbol_string.get_cigar_fragment(j, j + read_len)).count_elements('M')
                fm_span1[i].append(fm_pos1[i][-1] + span_dif)
                md_so_far += 1
        else:
            pass
    print(f'For new method, ploid {i}, time = {time()-start_time + cigar_string_new_time}')

    start_time = time()
    for j in range(len(temp_symbol_string2)):
        fm_pos2[i].append(md_so_far)
        # fix an edge case with deletions
        if temp_symbol_string2[j] == 'D':
            fm_pos2[i][-1] += temp_symbol_string2[j].count('D')
        # compute number of ref matches for each read
        span_dif = len([n for n in temp_symbol_string2[j:j + read_len] if n == 'M'])
        fm_span2[i].append(fm_pos2[i][-1] + span_dif)
        md_so_far += temp_symbol_string2[j].count('M') + temp_symbol_string2[j].count('D')
    print(f'For old method, ploid {i}, time = {time() - start_time + old_method_time}')


assert(fm_pos1 == fm_pos2)
assert(fm_span1 == fm_span2)
