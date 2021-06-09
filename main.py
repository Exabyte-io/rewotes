import os
import schedule
import sys
import time

import scanner
import client

def is_testing():
    return len(sys.argv) > 1 and '--test' in sys.argv[1:]

def sync():
    c = client.Client(os.environ.get('AWS_ID', 0),
                      os.environ.get('AWS_KEY', 0),
                      os.environ.get('AWS_BUCKET', "main"),
                      is_testing())
    s = scanner.Scanner(c)
    s.sync()

def main():
    # Schedule a sync every hour.
    schedule.every().hour.do(sync)

    # For testing, attempt a sync right away.
    if is_testing():
        sync()
        print("Done syncing.")

    while True:
        schedule.run_pending()
        time.sleep(60)


if __name__ == '__main__':
    main()