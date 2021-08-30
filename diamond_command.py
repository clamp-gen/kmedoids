from Bio.Application import _Option, AbstractCommandline, _Switch

class DiamondCommandline(AbstractCommandline):
    """Base Commandline object for (new) Diamond wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the diamond makedb and blastp tools.
    """

    def __init__(self, cmd="diamond", **kwargs):
        assert cmd is not None
        self.parameters= [
            # Core:
            _Switch(
                ["help", "h"],
                "Print USAGE, DESCRIPTION and ARGUMENTS description; "
                "ignore other arguments.",
            ),
            _Switch(
                ["blastp", "bp"],
                "use blastp.",
            ),
            _Switch(
                ["makedb", "createdb"],
                "use makedb .",
            ),
            _Switch(
                ["--sensitive", "sensitive"],
                "Use sensitive mode.",
            ),
            _Switch(
                ["--more-sensitive", "more_sensitive"],
                " Use more sensitive mode.",
            ),
            # Output configuration options
            _Option(
                ["--query", "query"],
                "The sequence to search with.",
                filename=True,
                equate=False,
            ),  # Should this be required?
            _Option(
                ["--in", "infile"],
                "input fasta for db.",
                filename=True,
                equate=False,
            ),  # Should this be required?
            _Option(
                ["--out", "out"],
                "Output file for alignment.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-d", "db"],
                "query db",
                filename=True,
                equate=False,
            ),
            # Formatting options:
            _Option(
                ["-f", "outfmt"],
                "Alignment view.  Typically an integer 0-14 but for some "
                "formats can be named columns like 'BLAST tabular.  "
                "Use 5 for XML output. ",
                filename=True,  # to ensure spaced inputs are quoted
                equate=False,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

    def _validate_incompatibilities(self, incompatibles):
        """Validate parameters for incompatibilities (PRIVATE).

        Used by the _validate method.
        """
        for a in incompatibles:
            if self._get_parameter(a):
                for b in incompatibles[a]:
                    if self._get_parameter(b):
                        raise ValueError("Options %s and %s are incompatible." % (a, b))

# for test
def main():
    output = DiamondblastpCommandline(bp=True, query="seq1.fasta", db="seq1", outfmt=5)
    print(output)
    output()


if __name__ == "__main__":
    main()