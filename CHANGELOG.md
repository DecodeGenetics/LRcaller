# Changelog

## v1.0

### Results

  * **Compared with v0.2, results have changed in some cases.**
  * Various bugs have been fixed that slighty affected results, improving the results.
  * Alignment is now performed with a band by default; in rare cases this might degrade results. This may be changed/deactivated via the `--band` option.

### Performance

  * Significant performance improvements; **speed-ups of 3x to more than 10x** (depends on data and options).
  * To achieve full speed, you need to build with machine-specific optimisations (see build instructions).

### Features

  * Dynamically adapt the window size to the length of the current read (via `--dyn-w-size`).
  * Automatically cache all input files in a local folder to reduce load on remote ressources (via `--cache-data-in-tmp`).

### Command-line interface

  * The option `--bam2` has been removed; instead just give a list of the combined files as single argument.
  * The options that influence scoring in alignment have been deprecated. They are no longer visible on the command-line and will be removed in a future update.
  * Please change all scripts to use the long form of the options (those that begin with `--`), because a future update will remove/change the short-forms.
  * The option `--dyn-w-size` was added (see above).
  * The option `--cache-data-in-tmp` was added (see above).

### Under-the-hood

  * Extensive refactoring, modernisation and clean-up of code.
  * Continuous Integration via Github Actions.
  * Code-formatting via clang-format.
