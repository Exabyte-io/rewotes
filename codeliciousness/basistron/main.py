# -*- coding: utf-8 -*-


if __name__ == "__main__":
    from basistron import cli
    parser = cli.get_parser()
    driver = cli.process_args(parser.parse_args())
