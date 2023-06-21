# Unit Testing

The `testing` folder contains unit tests using [GoogleTest](https://google.github.io/googletest/primer.html). Unit tests are available for the following classes: `reaction_network`, `nano_particle`, `lattice_reaction_network`, and `lattice`. The unit tests can be used to help determine if changes to the open-source code introduce bugs but are not completely comprehensive.

## Makefile 
- Create executables for each test 
    ```
    $ make all
    ```
- Create single executables for class foo
    ```
    $ make foo_test
    ```
- Clean
    ```
    $ make clean
    ```

## Output

Each unit test will either fail or pass. GoogleTest can use colors in its terminal output to make it easier to spot the important information:

<pre>
... tests that pass ...

<font color="green">[----------]</font> 1 test from FooTest
<font color="green">[ RUN      ]</font> FooTest.DoesAbc
<font color="green">[       OK ]</font> FooTest.DoesAbc
<font color="green">[----------]</font> 2 tests from BarTest
<font color="green">[ RUN      ]</font> BarTest.HasXyzProperty
<font color="green">[       OK ]</font> BarTest.HasXyzProperty
<font color="green">[ RUN      ]</font> BarTest.ReturnsTrueOnSuccess

... tests that fail ...

<font color="red">[   FAILED ]</font> BarTest.ReturnsTrueOnSuccess
...
<font color="green">[==========]</font> 30 tests from 14 test suites ran.
<font color="green">[   PASSED ]</font> 28 tests.
<font color="red">[   FAILED ]</font> 2 tests, listed below:
<font color="red">[   FAILED ]</font> BarTest.ReturnsTrueOnSuccess
<font color="red">[   FAILED ]</font> AnotherTest.DoesXyz

 2 FAILED TESTS
</pre>