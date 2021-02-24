import random


def generate_random_dna(lnth: int, seed: int = None) -> str:
    """
    Takes a parameter length and returns a randomly generated DNA string of that length
    :param lnth: how long of a string to generate
    :param seed: Optional seed to produce reproducibly random results
    :return: randomly generated string
    """
    set = ["A", "G", "C", "T"]
    if seed:
        random.seed(seed)
    else:
        random.seed()
    ret = ""
    for i in range(lnth):
        ret += random.choice(set)
    return ret


if __name__ == '__main__':
    print(generate_random_dna(10))
    print(generate_random_dna(10, 1))
    print(generate_random_dna(10, 1))
