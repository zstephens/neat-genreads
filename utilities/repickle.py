import sys
import os
import pathlib
import glob
import pickle


def main():
    dir_to_repickle = pathlib.Path(sys.argv[1])

    if not dir_to_repickle.is_dir():
        print("Input is not a directory.")
        sys.exit(1)

    os.chdir(dir_to_repickle)
    for file in glob.glob("*.p"):
        data = pickle.load(open(file, 'rb'), encoding="bytes")
        pickle.dump(data, open(file, "wb"), pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    main()
