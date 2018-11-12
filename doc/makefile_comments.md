Notes on the Makefile:

* The TODO at the top reminds me that I am using a different version of a
library in development and it must be removed before deployment.
* The `TARGET` is the main executable of the project, in this case `bin/runner`.
Type make and this is what gets built.
* I’m using *g++* because it’s the same on Mac OS X and on the production
Linux boxes.
* If I uncomment the clang line, I get a failed link as the libraries are
incompatible (or comment out the last line under `$(TARGET):`). But then I get
the benefit of a clang static analyzer run help me make my code better, well
worth it.
* I use the fewest number of compiler `CFLAGS` when developing as possible,
optimization happens later.
* The `SOURCES` list is dynamic, I don’t want to manually have to maintain this
list as I program. Anything in the `src` folder will be included in the compile
as long as it has a `SRCEXT` extension.
* The `OBJECTS` list is also dynamic and uses a Makefile trick to build the
list based on available sources.
* The `LIB` in this case uses a local library for MongoDB as I am testing it,
but uses the default homebrew or yum installed libraries for boost. I normally
do not use boost, but Mongo needs it.
* The `INC` ensures all headers in the include folder are accessible.
* I like to see the commands that run, hence the multitude of `@echo`'s.
* Since there are so few of them, I manually add spikes and test builds as a
new Makefile target, see the ticket: target for example.
* The `.PHONY` clean is brilliant, it nukes the build folder and the main
executable. It does not clean spike or test executables though.

Posted by *Hilton Lipschitz*