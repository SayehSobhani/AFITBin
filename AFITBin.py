from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('filename')           # positional argument
parser.add_argument('-c', '--count')      # option that takes a value
parser.add_argument('-v', '--verbose',
                    action='store_true')
args = parser.parse_args()
print(args.filename, args.count, args.verbose)