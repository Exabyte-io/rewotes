import os
import subprocess
import pathlib
import argparse
import pathlib
import sys
import random


PIPE = "│"
ELBOW = "└──"
TEE = "├──"
PIPE_PREFIX = "│   "
SPACE_PREFIX = "    "

class _TreeGenerator:
    def __init__(self, root_dir, tree_depth):
        self._root_dir = pathlib.Path(root_dir)
        self.tree_depth = tree_depth
        self._tree = []

    def build_tree(self):
        self._tree_head()
        self._tree_body(self._root_dir)
        return self._tree

    def _tree_head(self):
        self._tree.append(f"{self._root_dir}{os.sep}")
        self._tree.append(PIPE)

    def _tree_body(self, directory, prefix=""):
        entries = directory.iterdir()
        entries = sorted(entries, key=lambda entry: entry.is_file())
        entries_count = len(entries)
        for index, entry in enumerate(entries):
            connector = ELBOW if index == entries_count - 1 else TEE
            if entry.is_dir():
                self._add_directory(
                    entry, index, entries_count, prefix, connector
                )
            else:
                self._add_file(entry, prefix, connector)
                
    def _add_directory(
        self, directory, index, entries_count, prefix, connector
    ):
        self._tree.append(f"{prefix}{connector} {directory.name}{os.sep}")
        if index != entries_count - 1:
            prefix += PIPE_PREFIX
        else:
            prefix += SPACE_PREFIX
        self._tree_body(
            directory=directory,
            prefix=prefix,
        )
        self._tree.append(prefix.rstrip())

    def _add_file(self, file, prefix, connector):
        self._tree.append(f"{prefix}{connector} {file.name}")

class DirectoryTree:
    def __init__(self, root_dir, tree_depth):
        self.root_dir = root_dir
        self.tree_depth = tree_depth
        self._generator = _TreeGenerator(root_dir, tree_depth)

    # def generate(self):
    #     tree = self._generator.build_tree()
    #     for entry in tree:
    #         print(entry)

    def _generate_file(self, filename):
        sizetypes = ['b','k','m']
        sizetype = sizetypes[random.randint(0,len(sizetypes)-1)]
        size = random.randint(8,100)
        args = f"{size}{sizetype}"
        cmd = ['mkfile', "-n", args, filename]
        print(cmd)
        process = subprocess.run(['mkfile', "-n", args, filename])

    def _generate(self, path, dir_level):
        print(f"makedir: {path}")
        os.makedirs(path)

        for file_num in range(self.tree_depth):
            _path = f"{path}{os.sep}{file_num}"
            self._generate_file(_path)
            
        if dir_level < self.tree_depth:
            _path = f"{path}{os.sep}d{dir_level}"
            self._generate(_path, dir_level + 1)
            
    def generate(self):
        # os.makedirs(self.root_dir)
        for level in range(self.tree_depth):
            path = f"{self.root_dir}{os.sep}d{level}"
            self._generate(path, 1)
            
def parse_cmd_line_arguments():
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
