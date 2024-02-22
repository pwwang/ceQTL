from pipen import Pipen

from .processes import DataPreparation


class Pipeline(Pipen):
    name = "ceqtl_call"
    desc = "Calling co-expression QTLs"
    starts = DataPreparation


def main():
    Pipeline().run()
