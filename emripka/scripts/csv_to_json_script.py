from converters import *
import sys

# to-do: clean up pathing so can call this from any directory

def main():
    csv_to_json(
        MP_output_dir=r"../data/training/materialsproject_output/",
        MP_json_dir=r"../data/training/materialsproject_json/" ,
        fname=sys.argv[1]
    )

if __name__ == '__main__':
    main()
