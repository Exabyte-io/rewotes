import os
import random

RESOURCES_DIR = "./resources"

FILE_CONTENT = "hi1234"


def create_files(target_dir, num_files, max_depth):
    depth = 0
    current_dir = target_dir
    for num in range(0, num_files):
        random_num = random.randint(0, 100)
        if random_num < 20 and depth < max_depth:
            current_dir += "/child-dir-" + str(depth)
            if not os.path.exists(current_dir):
                os.mkdir(current_dir)
            depth += 1
        target_file_name = current_dir + "/" + str(num) + ".txt"
        with open(target_file_name, 'w') as f:
            f.write(FILE_CONTENT)


# Run script to create files for testing purposes
if __name__ == "__main__":
    create_files(RESOURCES_DIR, num_files=1000, max_depth=5)
