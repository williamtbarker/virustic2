# Security and data handling

Virustic2 accepts local sequence files and does not make network requests.
Treat all sequence inputs as potentially sensitive and review data ownership
before committing examples or results to a public repository.

Output JSON contains the input paths supplied on the command line. Remove or redact that field before
sharing reports if paths contain sample identifiers or protected project information.

Please report software vulnerabilities through a private GitHub security
advisory rather than a public issue.
