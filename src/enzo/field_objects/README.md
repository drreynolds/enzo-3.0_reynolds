# What Is This

This repository includes basic sketches of how FieldCollection objects and
FieldDescriptor objects might look inside Enzo.

# Roadmap

Once the tests pass and we have a reasonable amount of code coverage for tests
as well as functionality, this can be brought in to Enzo.

## What Fields Should Have

 * Centering
 * Dimensions
 * Rank
 * Left Edge (integers)
 * Interpolation Type
 * Units
 * Name

## What Fields Should Do

 * IO
 * Interpolate
 * Copy Overlapping
 * Binary math operations (add, sub, mul)
 * Unary math operations (min, max)
 * Copy
 * Allocate
 * Free
 * Convert to CGS

## This Repository

This repository will include .h and .C files that represent in-progress working
sketches of the field and field descriptor objects.  In the directory tests/
will be test files, using GoogleTest.  For the latest GoogleTest, visit

https://googletest.googlecode.com/

This might be it:

https://googletest.googlecode.com/files/gtest-1.6.0.zip

For ease of use, we suggest compiling and then symbolically linking the compile
directory here as gtest -- i.e.,

    ln -s /path/to/gtest-1.6.0 gtest

You can run the tests with:

    make runtests
