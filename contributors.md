# Guide for Developers and Contributors

## Getting Started

1. If you don't already have one, create and set up a free [GitHub account](https://github.com).
2. Install `git` on your local system, if it isn't already there.
3. Fork `RNMC`. Go to [the RNMC GitHub page](https://github.com/BlauGroup/RNMC) and select "Fork". This will create a copy of the repository that you own and control.
4. Copy the source code to your local system:

```
$ git clone https://github.com/<USERNAME>/RNMC
```

This will create a directory called "RNMC".

## Adding New Features

1. Make a new branch:

```
$ git pull  # Get most up-to-date code
$ git checkout -b <BRANCH_NAME>  # Create new branch and switch to it
```

2. Add the main repository as a remote:

```
$ git remote add upstream https://github.com/BlauGroup/RNMC
```

3. Code! Make commits frequently to save your progress, and push to sync your local changes to GitHub. The first time, you'll need to set an origin:

```
$ git push --set-upstream origin <BRANCH_NAME> 
```

4. When you're ready, create a [pull request](https://github.com/BlauGroup/RNMC/pulls) on GitHub. “Work-in-progress” pull requests are encouraged, especially if this is your first time contributing to `RNMC`. If your feature is not finished, put "\[WIP\]" in the pull request title.

## Raising Issues

If you have a problem installing or using `RNMC`, you can raise an issue on our [GitHub page](https://github.com/BlauGroup/RNMC/issues). The issues page is also an appropriate place to ask questions or to suggest a new feature. If you create an issue, you can label it to help us triage appropriately.

## Your First Contribution

If you have a feature that you want to add but don't know where to start: reach out to [Sam Blau](mailto:smblau@lbl.gov) and describe what you have in mind. He, or one of the developers, may be able to help you get started.

<!-- TODO: link to the page for expanding RNMC -->
If you're specifically interested in adding a new type of simulator to RNMC: check back soon! We're working on a guide to expanding `RNMC`.

If you want to help but don't have a particular idea in mind, check out the issues page and see if there's something that you would be interested in contributing.

## Coding Guidelines and Style

### Versions

`RNMC` is written in C++17. Before opening a pull request, please ensure that all modifications and additions are consistent with this version.

### Dependencies

`RNMC` is a slim codebase, with `sqlite3` and `GSL` being the only dependencies. In general, we encourage you to only use the C/C++ standard library and these existing dependencies and not to add any external dependencies. If you believe that an external dependency is necessary, please discuss this with the `RNMC` developers and maintainers.

### Documentation

TODO

### Testing

We test `RNMC` in two ways: unit tests and end-to-end tests. Unit tests evaluate discrete parts of the code base with self-contained examples. They're especially helpful for debugging. End-to-end tests perform full simulations (typically on small and simple systems) that ensure that the output of a simulation is correct, typically referencing some human-verified baseline results.

All new features must have unit tests. New simulators are encouraged to have end-to-end tests. These are not required if the simulators have strong non-deterministic behavior that make end-to-end correctness or consistency checking challenging (as is the case for certain simulations using `LGMC`).

See the [testing page](./Testing.md) for more detail about our testing suite.