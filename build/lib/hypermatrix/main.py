# hypermatrix/main.py

import argparse

def main():
    parser = argparse.ArgumentParser(description='Hypermatrix command-line tool')
    parser.add_argument('command', type=str, help='The command to run, e.g., ABcluster')
    parser.add_argument('--variable1', type=str, help='Description for variable1')
    parser.add_argument('--variable2', type=int, help='Description for variable2')
    # Add more arguments as needed

    args = parser.parse_args()

    # Dispatch the command
    if args.command == 'ABcluster':
        abcluster(args)
    else:
        print(f"Unknown command: {args.command}")

def abcluster(args):
    # Implement your ABcluster logic here, using args.variable1, args.variable2, etc.
    print(f"Running ABcluster with variable1={args.variable1} and variable2={args.variable2}")

if __name__ == "__main__":
    main()

