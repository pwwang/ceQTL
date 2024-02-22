from pipen import Pipen

from .processes import ExampleProcess


class Pipeline(Pipen):
    name = "ceqtl_preproc"
    desc = "Preprocess and prepre data for ceqtl pipeline"
    starts = ExampleProcess


def main():
    Pipeline().run()
