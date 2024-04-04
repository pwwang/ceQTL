import sys
import argx

commands = {
    "call": "Call the ceQTLs",
    "preproc": "Preprocess the data for the ceQTL pipeline",
    "simulate": "Simulate data for testing the ceQTL pipeline",
    "genoqc": "Quality control for genotype data",
    "exprqc": "Quality control for expression data",
    "eqtl": "Run eQTL analysis",
    "roc": "Run ROC analysis",
}

parser = argx.ArgumentParser(description="ceQTL pipeline and tools.")
for name, help in commands.items():
    command = parser.add_command(name, help=help, prefix_chars="#")
    command.add_argument(
        "arguments",
        nargs="*",
        help="Arguments to pass to the sub-pipeline",
    )


def main():
    args = parser.parse_args()

    if args.COMMAND not in commands:
        print("Unknown command: {}".format(args.COMMAND))
        parser.print_help()
        exit(1)

    sys.argv = [f"{sys.argv[0]} {args.COMMAND}"] + args.arguments
    module = f"{args.COMMAND}.main"
    try:
        imported = __import__(module, fromlist=["main"])
    except ModuleNotFoundError:
        imported = __import__(args.COMMAND)

    imported.main()
