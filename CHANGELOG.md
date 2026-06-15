# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- **`wgs`: replaced `--exclude-unpaired-reads` with `--include-unpaired-reads`.**
  The old flag was a boolean that defaulted to `true` with no way to negate it
  from the command line, so unpaired reads (and reads with an unmapped mate)
  were always excluded and the option could never actually be toggled. The new
  `--include-unpaired-reads` flag defaults to `false` (behavior unchanged) and,
  when passed, counts unpaired reads toward coverage. This is a breaking change
  to the CLI, but the removed flag had no working behavior to preserve.
  ([#34](https://github.com/fulcrumgenomics/riker/issues/34),
  [#35](https://github.com/fulcrumgenomics/riker/issues/35))

## [0.2.0]

Initial documented release. See the
[GitHub release notes](https://github.com/fulcrumgenomics/riker/releases) for
the history prior to this changelog.
