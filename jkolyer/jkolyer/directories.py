"""DirectoryTree: utility for generating sample files for uploading.
   Uses `mkfile` command line utility.
"""
import os
import subprocess
import pathlib
import argparse
import pathlib
import sys
import random

class DirectoryTree:
    def __init__(self, root_dir, tree_depth):
        """Instance constructor.  
        :param root_dir: parent directory of files
        :param tree_depth: number of tree levels
        :return: type describe
        """
        self.root_dir = root_dir
        self.tree_depth = tree_depth

    def _generate_file(self, file_path):
        """Utility for generating a file at the given path.
           Selects a random file size under 100 of units
           bytes, kilobytes, megabytes (chosen at random).
        :param file_path: file location relative to root
        """
        sizetypes = ['b','k','m']
        sizetype = sizetypes[random.randint(0,len(sizetypes)-1)]
        size = random.randint(8,100)
        args = f"{size}{sizetype}"
        cmd = ['mkfile', "-n", args, file_path]
        print(cmd)
        process = subprocess.run(['mkfile', "-n", args, file_path])

    def _generate(self, path, dir_level):
        """Iterates over `tree_depth` levels.  Calls itself
           recursively if `dir_level` < `tree_depth`.
           Creates directory at given path.
        :param path: directory to insert child files/dirs
        :param dir_level: current level
        :return: type describe
        """
        print(f"makedir: {path}")
        os.makedirs(path)

        for file_num in range(self.tree_depth):
            _path = f"{path}{os.sep}{file_num}"
            self._generate_file(_path)
            
        if dir_level < self.tree_depth:
            _path = f"{path}{os.sep}d{dir_level}"
            self._generate(_path, dir_level + 1)
            
    def generate(self):
        for level in range(self.tree_depth):
            path = f"{self.root_dir}{os.sep}d{level}"
            self._generate(path, 1)
            
def parse_cmd_line_arguments():
    """Parses command line for `root_dir` and `tree_depth`
       Expected usage: `python directories.py --tree_depth N ROOT_DIR
    :return: parsed arguments
    """
    parser = argparse.ArgumentParser(
        prog="tree",
        description="RP Tree, a directory tree generator",
        epilog="Thanks for using RP Tree!",
    )
    parser.add_argument(
        "root_dir",
        metavar="ROOT_DIR",
        nargs="?",
        default=".",
        help="Generate a full directory tree starting at ROOT_DIR",
    )
    parser.add_argument(
        "--tree_depth",
        metavar="TREE_DEPTH",
        nargs=1,
        default="3",
        type=int,
        required=True,
        help="How many directory levels to create",
    )
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_cmd_line_arguments()
    root_dir = pathlib.Path(args.root_dir)
    tree_depth = args.tree_depth
    if not root_dir.is_dir():
        print("The specified root directory doesn't exist")
        sys.exit()
    tree = DirectoryTree(root_dir, tree_depth[0])
    tree.generate()
