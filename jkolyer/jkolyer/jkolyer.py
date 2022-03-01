import pathlib
import argparse
import pathlib
import sys

"""
Commands
flush database
flush samples
perform upload
"""


def parse_cmd_line_arguments():
    """Describe
    :param name: describe
    :param name: describe
    :return: type describe
    """
    parser = argparse.ArgumentParser(
        prog="tree",
        description="PFU: Parallel File Upload",
        epilog="Thanks for using the service!",
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


